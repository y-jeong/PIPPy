      program pip
      use iso_c_binding
#ifdef _USE_GPULIB
      use libgpu_fc_interface
#endif
c
      implicit double precision (a-h,p-z)
c MPI
      include 'mpif.h'
      integer my_id,nprocs,ierr
      integer status(MPI_STATUS_SIZE)

      include 'param.inc'
      dimension coef(maxcoef,maxsurf),sig(maxsurf,maxdata)
      dimension coef_1d(maxcoef*maxsurf)
      dimension basis(maxterm)
      dimension ind(maxterm,maxpair)

      dimension vvs(maxsurf2),vv(maxsurf,maxdata)
      dimension vve(maxsurf3,maxdata) ! energies only
      dimension rrr(maxdata,maxpair)
      dimension cut(100),nc(100),na(100),wc(100),wa(100),
     &          erra(100),errc(100)
      dimension x(maxatom),y(maxatom),z(maxatom)
      dimension sigtmp(maxsurf,maxdata),vvtmp(maxsurf,maxdata)
      dimension rrrtmp(maxdata,maxpair)
      !dimension basistmp(maxterm,maxdata)

      character*16 datafit,datatest
      character*2 dum
      !CHARACTER(LEN=16) :: datafit,datatest
      !CHARACTER(LEN=2) :: dum

      integer natom1,a,b,itmp(100)
      integer nsurf3,nsurf2,nsurf,isurf(maxsurf)
c      nsurf3: real # of surfaces (nsurf2 = nsurf3*nsurf3)
c      nsurf2: #surface total
c      nsurf:  #surfaces of interest
c      isurf:  indices of interest
      integer k,l 
c      k,l: real surface indices in (1, nsurf3)
      logical linteronly,lwrite(-1:100),lcust
!      logical fit_ex, test_ex
      logical ltraindata(maxsurf,maxdata),ltestdata(maxsurf,maxdata) ! YJ

#ifdef _USE_GPULIB
      dimension ind_1d(maxterm*maxpair)
      integer :: num_gpus, gpu_id
      type (c_ptr) :: gpu

      double precision, target :: v_
      double precision, pointer :: ptr_v
      type (c_ptr) :: addr_v
 
      dimension ibasis(maxterm)
#endif
      common/foox/rrr,nncoef,natom1,linter,lwrite,nterm,ibasis,npairs,
     &  ind
      linteronly = linter
      open(99,file='debug')

#ifdef _USE_GPULIB
      ptr_v => v_
      addr_v = c_loc(ptr_v)
      !write(99,*)"v_==",v_
      !write(99,*)"ptr_v==",ptr_v

      ! -- GPU setup (called once and at very beginning of program)
      gpu=libgpu_fc_create_device()
      num_gpus=libgpu_fc_get_num_devices(gpu)
      gpu_id=0
      call libgpu_fc_set_device(gpu, gpu_id)
      write(*,*) '# of devices detected= ',num_gpus,
     & ' running on device id= ',gpu_id
#endif

ccccc MPI
      call MPI_INIT(ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,my_id,ierr)

      IF (my_id.eq.0) THEN
        timestart=MPI_WTIME()
c        write(6,*)"Began at ",timestart
      ENDIF
!      write(6,"(a, i3)") " MPI current proc: ", my_id

      lcust = .false.
!      lcust = .true. ! turn on nonstandard options

      IF (my_id.eq.0) THEN
!      write(6,"(a, i3)") " MPI num procs: ", nprocs
      write(6,'(100("*"))')
      write(6,*)
!      write(6,*)"PIP: A Fortran code for fitting permutationally",
!     &  " invariant polynomials"
!      write(6,*)"     Daniel R. Moberg and Ahren W. Jasper, Argonne,",
!     & " 2021"
      write(6,*)"PIPPy: A massively parallel code for",
     & " permutationally invariant polynomial expansion",
     & " construction, evaluation, and optimization"
      write(6,*)"     Yeonjun Jeong, Christopher Knight, ",
     & "Michael J. Davis, and Ahren W. Jasper, Argonne National",
     & " Laboratory, 2024"
      write(6,*)

      !datafit='train.dat'
      !datatest='test.dat'
      open(unit=7,file='input')
      read(7,*)datafit,datatest
      read(7,*)nsurf3,nsurf2,nsurf,(isurf(j),j=1,nsurf)
      read(7,*)nwrite,(itmp(j),j=1,nwrite)
      read(7,*)sflag,wflag,epsilon,vvref
      read(7,*)ncut,(cut(j),j=1,ncut)
      ENDIF
      !write(99,*)nwrite
      !write(99,*)itmp

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call MPI_BCAST(datafit, 1, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(datatest,1, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

!      inquire(file=datafit, exist=fit_ex)
!      if (.not.fit_ex) then
!        write(6,*)"No file named ",datafit," found, exiting"
!        stop
!      endif

      call MPI_BCAST(epsilon, 1, MPI_DOUBLE_PRECISION, 0,
     &               MPI_COMM_WORLD, ierr)
      call MPI_BCAST(vvref, 1, MPI_DOUBLE_PRECISION, 0,
     &               MPI_COMM_WORLD, ierr)
      call MPI_BCAST(ncut, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(cut, ncut, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(nwrite, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
!      call MPI_BCAST(lwrite,nwrite+2,MPI_LOGICAL,0,MPI_COMM_WORLD, ierr)

      if (nwrite.eq.-1) then
        do i=1,100
          lwrite(i)=.true.
        enddo
      else if (nwrite.eq.0) then
        do i=1,100
          lwrite(i)=.false.
        enddo
      else
        do i=1,100
          lwrite(i)=.false.
        enddo
        do i=1,nwrite
          lwrite(itmp(i))=.true.
        enddo
      endif
      !write(99,*)lwrite
      !write(99,*)lwrite(20)

      call prepot ! generate or read basis

#ifdef _USE_GPULIB
      ! -- repack data...
      do i=1,nsurf
        do j=1,ncoef
          coef_1d((i-1)*ncoef+j)=coef(j,i)
        enddo
      enddo

      do i=1,nterm
        do j=1,npairs
          ind_1d((i-1)*npairs+j)=ind(i,j)
        enddo
      enddo
#endif

c      IF (.false.) THEN ! To not fit function
      IF (my_id.eq.0) THEN
      ncoef=nncoef
 
      open(7,file=datafit)
      read(7,*)ndat2
      if (ndat2.gt.maxdata) then
        write(6,*)"ndat = ",ndat2," while maxdata = ",maxdata
        stop
      endif
 
      !ndat=0
      do i=1,ndat2
        read(7,*)natom
        read(7,*)iz,dum,(vvs(j),j=1,nsurf2)
        do j=1,natom
          read(7,*)dum,x(j),y(j),z(j)
        enddo

c       Handle the surfaces of interest
        do j=1,nsurf
        !vvx=vvs(isurf(j))
        !vv(j,i)=vvx
        vv(j,i)=vvs(isurf(j))
        if (wflag.eq.2) then ! save the energy surfaces
          l=mod(isurf(j)-1,nsurf3) + 1 
          ! index of the second surface (energy only)
          k=(isurf(j)-l)/nsurf3 + 1 
          ! index of the first surface (energy only)
          m=(k-1)*nsurf3 + k ! index of the first energy surface
          n=(l-1)*nsurf3 + l ! index of the second energy surface
          vve(k,i)=vvs(m)
          vve(l,i)=vvs(n)
        endif
        enddo

c       Calculate the interatomic distances
        ii=0
        a=natom
        if (linteronly) a=natom1
        do j=1,a
          b=j
          if (linteronly) b=natom1
          do k=b+1,natom
            ii=ii+1
            rrr(i,ii)=dsqrt((x(j)-x(k))**2+(y(j)-y(k))**2
     &                                         +(z(j)-z(k))**2)
          enddo  
        enddo  
      enddo
      close(7)

c     Apply surface flag YJ
      if (sflag.eq.1) then ! fit |V_ij|'s
        do i=1,ndat2
        do j=1,nsurf
          vv(j,i)=dabs(vv(j,i))
        enddo
        enddo
      endif

c     Apply cuts
      do i=1,ndat2
      do j=1,nsurf
        vvx=vv(j,i)
        if ((vvx.lt.cut(1).and.vvx.gt.cut(ncut)).and.
     &   (.not.lcust.or.vvx.ne.0.d0)) then
          ltraindata(j,i)=.true.
        else
          ltraindata(j,i)=.false.
        endif
      enddo
      enddo

c     Assign weights YJ
      do i=1,ndat2
      do j=1,nsurf
      if (ltraindata(j,i)) then
        if (wflag.eq.0) then ! energy fitting
          sig(j,i)=1.d0/(epsilon/ 
     &         (dabs(vv(j,i)-vvref)+epsilon))
        else if (wflag.eq.1) then ! diabatic coupling fitting 1
          sig(j,i)=1.d0/(epsilon/ 
     &         (vvref+epsilon-dabs(vv(j,i))))
        else if (wflag.eq.2) then ! diabatic coupling fitting 2
          l=mod(isurf(j)-1,nsurf3) + 1 
          ! index of the second surface (energy only)
          k=(isurf(j)-l)/nsurf3 + 1 
          ! index of the first surface (energy only)
          sig(j,i)=1.d0/(epsilon/ 
     &         (vvref+epsilon-dabs(vv(j,i))**2/dabs(vve(k,i)-vve(l,i))))
        else
          stop "Invalid wflag"
        endif
      else
        sig(j,i)=huge(0.d0) ! Unused
      endif
      enddo
      enddo

      write(6,'(34("*  "))')
      write(6,*)
      write(6,*)"Fitting to the training data in ",datafit
      if (datatest.eq."none".or.datatest.eq."skip") then 
        write(6,*)"Skipping out of sampling testing"
      else
        write(6,*)"Out of sample testing using data in ",datatest
      endif
      write(6,*)"Weight function parameters: epsilon = ",
     & epsilon," and vref = ",vvref
      write(6,*)

      open(56,file="coef.dat")
      !if (lwrite(20)) write(99,*)"Opening xmat.dat"
      if (lwrite(20)) open(20,file="xmat.dat")
      if (lwrite(11)) open(11,file="vtrain.dat")

c     Eliminate data out of the predefined energy range before fitting
c     YJ
      vvtmp=vv
      sigtmp=sig
      rrrtmp=rrr
      do i=1,nsurf
      ndat=0
      do j=1,ndat2
        if (ltraindata(i,j)) then
          ndat=ndat+1
          vv(i,ndat)=vvtmp(i,j)
          sig(i,ndat)=sigtmp(i,j)
          rrr(ndat,:)=rrrtmp(j,:)
        endif
      enddo 
      !write(99,*)"Surface ",isurf(i)
      !write(99,*)"ndat == ",ndat
      !write(99,*)vv(i,1:ndat2)
      !write(99,*)sig(i,1:ndat2)
 
      write(6,*)"Surface ",isurf(i)
      write(6,*)"Using ",ndat," of the ",ndat2," provided data"
      !write(99,*)"svdfit has started"
      !write(99,*)"coef before fitting:"
      !write(99,*)coef
      call svdfit(vv(i,:),sig(i,:),ndat,coef(:,i),ncoef,ndat,ncoef)
      !write(99,*)"coef after fitting:"
      !write(99,*)coef
      !write(99,*)"svdfit is done"
      do j=1,ncut
        erra(j)=0.d0 ! RMSE for E < cut(j)
        errc(j)=0.d0 ! RMSE for cut(j+1) < E < cut(j)
        wc(j)=0.d0
        wa(j)=0.d0
        nc(j)=0
        na(j)=0
      enddo

      vvimin=1d50
      vvxmin=1d50
      if (lwrite(20)) write(20,*)" INDEX SIGMA VAI BASIS(1...NCOEF)"
      if (lwrite(11)) write(11,'(a10,i5)')"Surface ",isurf(i)
      if (lwrite(11)) write(11,'(a10,3a18)')"INDEX","SIGMA","VAI","VFIT"

      !if (i.eq.1) then
      !  do j=1,ndat
      !    call funcs1(j,basis,ncoef) ! Calculate basis for geometry j
      !    basistmp(:,j)=basis
      !  enddo
      !endif

#ifdef _USE_GPULIB
      ! -- GPU memory setup (called once and after data initialized)
!      write(99,*)coef(:,i)
!      call libgpu_fc_setup(gpu, rrr, ind_1d, coef(:,i), ibasis,
!     &   nterm, nncoef, npairs)
      call libgpu_fc_setup(gpu, ind_1d, coef_1d, ibasis,
     &   nterm, nncoef, npairs, nsurf)
#endif

      call cpu_time(t0)
      do j=1,ndat
        !basis=basistmp(:,j)
#if _USE_GPULIB
        ! -- GPU compute (called many times with different r)
        !call libgpu_fc_compute_v(gpu, rrr(j,:), addr_v)
        call libgpu_fc_compute_v(gpu, rrr(j,:), addr_v, i)
        vvx=v_*1.d0
#else
        call funcs1(j,basis,ncoef) ! Calculate basis for geometry j
        vvx=0.d0
        do k=1,ncoef
          vvx=vvx+coef(k,i)*basis(k)
          if (j.eq.1) write(56,'(i5,e20.10)')k,coef(k,i)
        enddo
#endif
        if (i.eq.1) then
          if (lwrite(20)) write(20,'(i10,*(e18.6))')j,1./sig(i,j),
     &   vv(i,j),(basis(k),k=1,ncoef)
        endif
        !write(99,*)"i",i
        !write(99,*)"1/sig(i,j): ",1/sig(i,j)
        !write(99,*)"vv(i,j): ",vv(i,j)
        !write(99,*)"vvx: ",vvx
        if (lwrite(11)) write(11,'(i2,i8,3e18.6)')i,j,1./sig(i,j),
     &   vv(i,j),vvx
        do k=1,ncut
          if (vv(i,j).lt.cut(k)) then
            erra(k)=erra(k)+(vvx-vv(i,j))**2/sig(i,j)**2
            wa(k)=wa(k)+1.d0/sig(i,j)**2
            na(k)=na(k)+1
          endif
        enddo
        do k=1,ncut-1
          if (vv(i,j).lt.cut(k).and.vv(i,j).ge.cut(k+1)) then
            errc(k)=errc(k)+(vvx-vv(i,j))**2/sig(i,j)**2
            wc(k)=wc(k)+1.d0/sig(i,j)**2
            nc(k)=nc(k)+1
          endif
        enddo
        if (vvx.lt.vvxmin) then
          vvxmin=vvx
          vvxa=vv(i,j)
          vvxb=vvx
        endif
        if (vv(i,j).lt.vvimin) then
          vvimin=vv(i,j)
          vvia=vv(i,j)
          vvib=vvx
        endif
      enddo ! ndat
      call cpu_time(t1)
      print '(a,f12.5,a)', "train set evaluation time: ",
     & (t1-t0)*1000, " ms"

      errc(ncut)=erra(ncut)
      wc(ncut)=wa(ncut)
      nc(ncut)=na(ncut)
      write(6,*)'Weighted errors between (below) user-provided energies'
      write(6,*)'          E       number    %weight       error',
     &                  '   (     number        error   )'
      do k=1,ncut
        if (wa(k).ne.0.d0) erra(k)=dsqrt(erra(k)/wa(k))
        if (wc(k).ne.0.d0) errc(k)=dsqrt(errc(k)/wc(k))
        write(6,'(f15.6,i10,f10.1,f15.5," ( ",i10,f15.5," ) ")')
     &    cut(k),nc(k),wc(k)/wa(1)*100.,errc(k),na(k),erra(k)
      enddo
      write(6,*) 
      write(6,*) "Comparison of low energy points found while fitting"
      write(6,*) "                        data           fit        ",
     & "  difference "
      write(6,'(a,3f15.5)')"   data minimum ",vvia,vvib,vvia-vvib
      write(6,'(a,3f15.5)')" fitted minumum ",vvxa,vvxb,vvxa-vvxb
      write(6,*)
      enddo ! nsurf
      close(56)
 
c     test set
      if (datatest.eq."none".or.datatest.eq."skip") then ! skip
      else
!      inquire(file=datatest, exist=test_ex)
!      if (.not.test_ex) then
!        write(6,*)"No file named ",datatest," found"
!        write(6,*)"Skipping test set."
!        go to 290
!      endif
      open(7,file=datatest)
      rewind(7)
      read(7,*)ndat2
      if (ndat2.gt.maxdata) then
        write(6,*)"ndat = ",ndat2," while maxdata = ",maxdata
        stop
      endif
  
      ndat=0
      do i=1,ndat2
        read(7,*)
        read(7,*)iz,dum,(vvs(j),j=1,nsurf2)
        do j=1,natom
          read(7,*)dum,x(j),y(j),z(j)
        enddo
 
        do j=1,nsurf
        !vvx=vvs(isurf(j))
        vv(j,i)=vvs(isurf(j))
        if (wflag.eq.2) then ! save the energy surfaces
          l=mod(isurf(j)-1,nsurf3) + 1 
          ! index of the second surface (energy only)
          k=(isurf(j)-l)/nsurf3 + 1 
          ! index of the first surface (energy only)
          m=(k-1)*nsurf3 + k ! index of the first energy surface
          n=(l-1)*nsurf3 + l ! index of the second energy surface
          vve(k,i)=vvs(m)
          vve(l,i)=vvs(n)
        endif
        enddo

        ii=0
        a=natom
        if (linteronly) a=natom1
        do j=1,a
          b=j
          if (linteronly) b=natom1
          do k=b+1,natom
            ii=ii+1
      rrr(i,ii)=dsqrt((x(j)-x(k))**2+(y(j)-y(k))**2+(z(j)-z(k))**2)
          enddo  
        enddo
      enddo  

c     Apply surface flag YJ
      if (sflag.eq.1) then ! examine fit to |V_ij|'s
        do i=1,ndat2
        do j=1,nsurf
          vv(j,i)=dabs(vv(j,i))
        enddo
        enddo
      endif

c     Apply cuts
      xcut=cut(1) 
      if (lcust.and.ncut.gt.1) xcut=cut(2)
      do i=1,ndat2
      do j=1,nsurf
        vvx=vv(j,i)
        if ((vvx.lt.xcut.and.vvx.gt.cut(ncut)).and.
     &   (.not.lcust.or.vvx.ne.0.d0)) then
          ltestdata(j,i)=.true.
        else
          ltestdata(j,i)=.false.
        endif
      enddo
      enddo

c     Assign weights YJ
      do i=1,ndat2
      do j=1,nsurf
      if (ltestdata(j,i)) then
        if (wflag.eq.0) then ! energy fitting
          sig(j,i)=1.d0/(epsilon/ 
     &         (dabs(vv(j,i)-vvref)+epsilon))
        else if (wflag.eq.1) then ! diabatic coupling fitting 1
          sig(j,i)=1.d0/(epsilon/ 
     &         (vvref+epsilon-dabs(vv(j,i))))
        else if (wflag.eq.2) then ! diabatic coupling fitting 2
          l=mod(isurf(j)-1,nsurf3) + 1 
          ! index of the second surface (energy only)
          k=(isurf(j)-l)/nsurf3 + 1 
          ! index of the first surface (energy only)
          sig(j,i)=1.d0/(epsilon/ 
     &         (vvref+epsilon-dabs(vv(j,i))**2/dabs(vve(k,i)-vve(l,i))))
        else
          stop "Invalid wflag"
        endif
      else
        sig(j,i)=huge(0.d0) ! Unused
      endif
      enddo
      enddo

c     Test set: Eliminate data out of the predefined energy
c     YJ
      vvtmp=vv
      sigtmp=sig
      rrrtmp=rrr
      do i=1,nsurf
      ndat=0
        do j=1,ndat2
          if (ltestdata(i,j)) then
            ndat=ndat+1
            vv(i,ndat)=vvtmp(i,j)
            sig(i,ndat)=sigtmp(i,j)
            rrr(ndat,:)=rrrtmp(j,:)
          endif
        enddo 
        write(6,*)"Test set surface ",isurf(i),": Using ",ndat,
     & " of the ",ndat2," provided data"
 
      err3=0.d0
      !ndat3=0
      wn=0.d0
      if (lwrite(12)) open(12,file="vtest.dat")
      if (lwrite(12)) write(12,'(a10,3a18)')"INDEX","SIGMA","VAI","VFIT"

      !if (i.eq.1) then
      !  do j=1,ndat
      !    call funcs1(j,basis,ncoef)
      !    basistmp(:,j)=basis
      !  enddo
      !endif

#ifdef _USE_GPULIB
      ! -- GPU memory setup (called once and after data initialized)
!      write(99,*)coef(:,i)
!      call libgpu_fc_setup(gpu, rrr, ind_1d, coef(:,i), ibasis,
!     &   nterm, nncoef, npairs)
      call libgpu_fc_setup(gpu, ind_1d, coef_1d, ibasis,
     &   nterm, nncoef, npairs, nsurf)
#endif

      call cpu_time(t2)
      do j=1,ndat
        !basis=basistmp(:,j)
#ifdef _USE_GPULIB
        ! -- GPU compute (called many times with different r)
        !call libgpu_fc_compute_v(gpu, rrr(j,:), addr_v)
        call libgpu_fc_compute_v(gpu, rrr(j,:), addr_v, i)
        vvx=v_*1.d0
#else
        call funcs1(j,basis,ncoef)
        vvx=0.d0
        do k=1,ncoef
          vvx=vvx+coef(k,i)*basis(k)
        enddo
#endif
        err3=err3+(vvx-vv(i,j))**2/sig(i,j)**2
        wn=wn+1.d0/sig(i,j)**2
        if (lwrite(12)) write(12,'(i2,i8,3e18.6)')i,j,1./sig(i,j),
     &   vv(i,j),vvx
      enddo
      enddo ! nsurf

#ifdef _USE_GPULIB
      call libgpu_fc_destroy_device(gpu)
#endif

      err3=dsqrt(err3/wn)
      call cpu_time(t3)
      print '(a,f12.5,a)', "test set evaluation time: ",
     & (t3-t2)*1000, " ms"
 
      write(6,*)
      write(6,*)"Out of sample test set error: ",err3 
      endif ! test set
! 290  continue

      write(6,*)
      write(6,*)"Optimized coefficients written to coef.dat"
      write(6,*)
      write(6,'(100("*"))')

      if (lcust) then
        ix=1
        if (ncut.gt.1) ix=2 
        write(6,*) 
        write(6,*) "summary",ncoef,erra(ix),err3
        write(6,*) 
      endif
 
      ENDIF ! serial section for fitting procedure

      IF (my_id.eq.0) THEN
        timeend=MPI_WTIME()
        write(6,'(/,"Total time: ",f20.10," s",/)'),(timeend-timestart)
      ENDIF

      call MPI_FINALIZE(ierr)

      close(99)
      end
