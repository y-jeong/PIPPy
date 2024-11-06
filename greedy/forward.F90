program greedy

#if _USE_MAGMA
use iso_c_binding
use magma2
#endif

implicit double precision(a-h,o-z)
character*2 dum

! PIPPy is used to write the B matrix as defined in eq 2
! This code assumes the B matrix and the training data and weights are in xmat.dat
! This code also reads basis.dat to get the size of the full basis, ncoef
! The termination condition is hard coded such that it runs until it uses 1/2 the basis.

! set these
logical :: ldofull = .false.   ! Compute error for the full B matrix?
logical :: lreadin = .false.    ! Start with basis in greedy.dat?
integer :: ntarget = 10       ! End when you reach ntarget basis functions

character(len=10) :: arg
double precision, dimension(:), allocatable :: ww,vv  ! w(i) & v(i), weights and training energies
double precision, dimension(:,:), allocatable :: bb   ! b-matrix(i,k)
! i is used to index training data points w/ max ndat determined by reading xmat.dat
! k is used to index terms and permutation groups w/ max ncoef read from basis.dat

! least-squares
double precision, dimension(:), allocatable :: xty  ! X^T*y(k)
double precision, dimension(:), allocatable :: cc  ! c(k), coefficients
double precision, dimension(:,:), allocatable :: cov  ! cov. matrix (X^TX)^-1 (k1,k2)

integer, dimension(:), allocatable :: ipiv ! lapack
double precision, dimension(:), allocatable :: work ! lapack

integer, dimension(:), allocatable :: used,list ! flag and list of groups added so far
double precision, dimension(:), allocatable :: u1,u2,u3 ! efficient updating scheme variables (k)
double precision, dimension(:,:), allocatable :: covu ! temp cov. matrix(k1,k2)

#if defined(_USE_LAPACK) || defined(_USE_LAPACK_TEST) || defined(_USE_MAGMA)
double precision, dimension(:,:), allocatable :: bba  ! b-matrix, compressed
double precision, dimension(:), allocatable :: vfits, uvfits
double precision, dimension(:), allocatable :: uvv
double precision, dimension(:), allocatable :: errs
#endif

#if _USE_MAGMA
integer :: dev
type(c_ptr) :: queue !! magma_queue_t
type(c_ptr) :: dcov
type(c_ptr) :: dbba
type(c_ptr) :: dcc
type(c_ptr) :: dwork
type(c_ptr) :: dvfits

write(6,*)"--------------- init"
call magma_init()
    
write(6,*)"--------------- create queue"
dev = 0
call magma_queue_create( dev, queue )

#endif

! command line arguments
nargs = COMMAND_ARGUMENT_COUNT()
if (nargs.eq.3) then
call GET_COMMAND_ARGUMENT (1, arg)
read(arg,*)ldofull
call GET_COMMAND_ARGUMENT (2, arg)
read(arg,*)ntarget
call GET_COMMAND_ARGUMENT (3, arg)
read(arg,*)lreadin
else
print *,"Usage: forward.x LDOFULL NTARGET LREADIN"
stop
endif

! debug outputs
#if _USE_LAPACK
open(99,file='debuglapack')
#elif _USE_LAPACK_TEST
open(99,file='debuglapacktest')
#elif _USE_MAGMA
open(99,file='debuggpu')
#else
open(99,file='debugcpu')
#endif

call cpu_time(t1)
! read the basis
open(10,file="basis.dat")
read(10,*)natom,npair,ncoef,nterm
print *,"From basis.dat: natom,npair,ncoef,nterm = ",natom,npair,ncoef,nterm
close(10)
call cpu_time(t2)
write(99,*)"Time spent in reading the basis: ",t2-t1," s"
if (ntarget.eq.-1) ntarget=ncoef

! read the B matrix (maybe change to binary, as this matrix can be large)
call cpu_time(t3)
open(20,file="xmat.dat")
read(20,*)
ndat=0
do
  read(20,*,end = 20)
  ndat=ndat+1
end do
20 continue
print *,"Reading xmat.dat(ndat,ncoef) = B(",ndat,ncoef,")"
print *
rewind(20)
allocate(ww(ndat))
allocate(vv(ndat))
allocate(bb(ndat,ncoef))
#if defined(_USE_LAPACK) || defined(_USE_LAPACK_TEST) || defined(_USE_MAGMA)
allocate(bba(ndat,ncoef))
allocate(vfits(ndat))
allocate(uvfits(ndat))
allocate(uvv(ndat))
allocate(errs(ndat))
#endif
read(20,*)
i=0
do
  i=i+1
  read(20,*,end = 21)ii,ww(i),vv(i),(bb(i,k),k=1,ncoef)
end do
21 continue
close(20)
call cpu_time(t4)
write(99,*)"Time spent in reading xmat.dat: ",t4-t3," s"

call cpu_time(t5)
! use weights
do i=1,ndat
vv(i)=vv(i)*ww(i)  ! y = W*V (eq 3)
do k=1,ncoef
bb(i,k)=bb(i,k)*ww(i) ! X = W*B (eq 4)
enddo
enddo
call cpu_time(t6)
write(99,*)"Time spent in y=W*V and X=W*B (1): ",t6-t5," s"

allocate(cov(ncoef,ncoef))
allocate(ipiv(ncoef))
allocate(work(ncoef))
allocate(xty(ncoef))
allocate(cc(ncoef))
#if _USE_MAGMA
info = magma_dmalloc( dcov, int(ncoef, c_size_t)*ncoef )
info = magma_dmalloc( dbba, int(ndat , c_size_t)*ncoef )
info = magma_dmalloc( dvfits, int(ndat , c_size_t) )
#endif

call cpu_time(t7)
do k=1,ncoef
xty(k)=0.d0
do i=1,ndat
xty(k)=xty(k)+bb(i,k)*vv(i) ! X^T*y in eq 5
enddo
enddo
call cpu_time(t8)
write(99,*)"Time spent in X^T*y (1): ",t8-t7," s"

! solve for the full B matrix
if (ldofull) then
! compute (X^T*X)^-1 = cov(k1,k2)
do k1=1,ncoef
do k2=1,ncoef
cov(k1,k2)=0.d0
do i=1,ndat
cov(k1,k2)=cov(k1,k2)+bb(i,k1)*bb(i,k2)  ! X^T*X
enddo
enddo
enddo
lwork=ncoef
call dgetrf(ncoef,ncoef,cov,ncoef,ipiv,info)       ! together, these two LAPACK routines
call dgetri(ncoef,cov,ncoef,ipiv,work,lwork,info)  ! return the inverse of X^T*X
! now, cov = (X^T*X)^-1, i.e., the first part of eq 5

! next we compute xty = X^T*y, i.e., the rest of eq 5
!do k=1,ncoef
!xty(k)=0.d0
!do i=1,ndat
!xty(k)=xty(k)+bb(i,k)*vv(i) ! X^T*y in eq 5
!enddo
!enddo

! compute coefficients
do k1=1,ncoef
cc(k1)=0.d0
do k2=1,ncoef
cc(k1)=cc(k1)+cov(k1,k2)*xty(k2)   ! c in eq 5
enddo
enddo

! test fit for the full B matrix
vf2=1000.
vv1=1000.
wrr=0.d0
err=0.d0
do i=1,ndat
vfit=0.d0
do k=1,ncoef
vfit=vfit+bb(i,k)/ww(i)*cc(k) ! eq 2 (remember to un-weight bb)
enddo
wrr=wrr+ww(i)**2
err=err+(ww(i)*(vfit-vv(i)/ww(i)))**2 ! least squares error
if (vfit.lt.vf2) then ! keep track of minimum fitted energy
vf1=vv(i)/ww(i)
vf2=vfit
endif
if (vv(i)/ww(i).lt.vv1) then ! keep track of minimum training energy
vv1=vv(i)/ww(i)
vv2=vfit
endif
enddo
err=dsqrt(err/wrr)
write(6,123)ncoef,err,0,vv1,vv2,vf1,vf2
endif











! INITIATE THE RUN
! used is a flag that will keep track of whether or not we have already added each group
! list will keep a growing list of which groups we have added
! nopt is the number of groups we added so far to our "optimized" expansion
nopt=0
allocate(used(ncoef))
allocate(list(ncoef))
do k=1,ncoef
used(k)=0
list(k)=0
enddo
nopt=1
used(1)=1
list(1)=1
call cpu_time(t15)
if (lreadin) then
! read in groups
! we will then add these 1 at a time to use our "forward" machinery
open(11,file='greedy.dat')
do while (.true.)
read(11,*,end=111)i,xdum,ii
if (ii.eq.0.or.ii.eq.1) then ! we always include ii=1 as our first basis func, we skip rows with ii=0
else if (ii.gt.0) then
nopt=nopt+1
list(nopt)=ii
else if (ii.lt.0) then
iii=0
do i=1,nopt
if (list(i).ne.abs(ii)) then
iii=iii+1
list(iii)=list(i)
endif
enddo
list(nopt)=0
nopt=nopt-1
endif
enddo
111 continue
close(11)
endif
call cpu_time(t16)
write(99,*)"Time spent in reading greedy.dat: ",t16-t15," s"
!nopt=1
! if nopt = 1 we add each one at a time
! if nopt is not reset, then it will start with the full read-in matrix

! evaluate the fit as above...
call cpu_time(t17)
do k1=1,nopt
used(list(k1))=1
do k2=1,nopt
cov(k1,k2)=0.d0
do i=1,ndat
cov(k1,k2)=cov(k1,k2)+bb(i,list(k1))*bb(i,list(k2)) ! X^T*X
enddo
enddo
enddo
lwork=nopt
call dgetrf(nopt,nopt,cov,ncoef,ipiv,info)           ! invert it using LAPACK
call dgetri(nopt,cov,ncoef,ipiv,work,lwork,info)
do k1=1,nopt
cc(k1)=0.d0
do k2=1,nopt
cc(k1)=cc(k1)+cov(k1,k2)*xty(list(k2)) ! c in eq 5
enddo
enddo

vf2=1000.
vv1=1000.
wrr=0.d0
err=0.d0
do i=1,ndat
vfit=0.d0
do k=1,nopt
vfit=vfit+bb(i,list(k))/ww(i)*cc(k)  ! eq 2 (remember to un-weight bb)
enddo
wrr=wrr+ww(i)**2
err=err+(ww(i)*(vfit-vv(i)/ww(i)))**2
if (vfit.lt.vf2) then
vf1=vv(i)/ww(i)
vf2=vfit
endif
if (vv(i)/ww(i).lt.vv1) then
vv1=vv(i)/ww(i)
vv2=vfit
endif
enddo
err=dsqrt(err/wrr)
if (nopt.eq.1) write(6,123)nopt,err,1,vv1,vv2,vf1,vf2
if (nopt.ne.1) write(6,123)nopt,err,-1,vv1,vv2,vf1,vf2

if (nopt.ge.ntarget) then
print *,"NOPT >= NTARGET (",nopt,">=",ntarget,") so there's nothing to do!"
stop
endif
call cpu_time(t18)
write(99,*)"Time spent in evaluating the fit from greedy.dat: ",t18-t17," s"


write(99,*)
write(99,*)" forward greedy algorithm"
! forward greedy algorithm
allocate(u1(ncoef))
allocate(u2(ncoef))
allocate(u3(ncoef))
allocate(covu(ncoef,ncoef))

100 continue
call cpu_time(tt1)
errmin=10000.
write(99,*)"nopt == ",nopt
tta=0.d0 ! average time for evaluating covu
ttb=0.d0 ! average time for evaluating the candidate fit
ttba=0.d0 ! a breakdown of ttb
ttbb=0.d0 ! a breakdown of ttb
do ii=1,ncoef ! main loop
call cpu_time(t19)
if (nopt.eq.ntarget) go to 999
if (list(nopt+1).ne.0.and.ii.ne.list(nopt+1)) cycle ! we know which basis to add because we read in the list
if (used(ii).eq.1) cycle ! we already used this group, so skip to the next group

! update cov = (X^T*X)^-1 using the efficient algorithm in emtiyaz.github.io/Writings/OneColInv.pdf
! notation follows those notes closely, where italic-B in those notes = our cov
! we generate the updated cov (called here covu) for each new group ii, then
!    compute the error for this expanded basis, and keep track of the one with
!    the lowest error. We permanently add that one and repeat the cycle. Groups
!    can't get added twice thanks to the used flag.
do k=1,nopt
u1(k)=0.d0
do i=1,ndat
u1(k)=u1(k)+bb(i,list(k))*bb(i,ii)
enddo
enddo

do k1=1,nopt
u2(k1)=0.d0
do k2=1,nopt
u2(k1)=u2(k1)+cov(k1,k2)*u1(k2)
enddo
enddo

d1=0.d0
do i=1,ndat
d1=d1+bb(i,ii)**2
enddo
d2=0.d0
do k=1,nopt
d2=d2+u1(k)*u2(k)
enddo
d0=1.d0/(d1-d2)
do k=1,nopt
u3(k)=u2(k)*d0
enddo

do k1=1,nopt
do k2=1,nopt
covu(k1,k2)=cov(k1,k2)+d0*u2(k1)*u2(k2)
enddo
covu(k1,nopt+1)=-u3(k1)   ! updated cov matrix 
covu(nopt+1,k1)=-u3(k1)
covu(nopt+1,nopt+1)=d0
enddo

do k1=1,nopt+1
cc(k1)=0.d0
do k2=1,nopt
cc(k1)=cc(k1)+covu(k1,k2)*xty(list(k2))
enddo
cc(k1)=cc(k1)+covu(k1,nopt+1)*xty(ii)
enddo
call cpu_time(t20)
!write(99,*)"Time spent in updating cov: ",t20-t19," s"
tta=tta+t20-t19
!write(99,*)"current taa == ",taa

call cpu_time(t21)
! test fit
 ! Is iterating over the first index faster?
 ! Might distribute (ncoef-nopt)*ndat vector-vector multiplications over gpus? This will require a significant rewriting of the code, though
#if defined(_USE_LAPACK) || defined(_USE_LAPACK_TEST) || defined(_USE_MAGMA)
! compress bb into bba
do k1=1,nopt
  do i=1,ndat
    bba(i,k1)=bb(i,list(k1))
  enddo
enddo
do i=1,ndat
  bba(i,nopt+1)=bb(i,ii)
enddo
#endif

#if defined(_USE_LAPACK) || defined(_USE_MAGMA)
vf2=1000.
vv1=1000.
wrr=0.d0
err=0.d0
call cpu_time(t211)
do i=1,ndat
  vfits(i)=0.d0
enddo
#if _USE_LAPACK
call dgemv('N',ndat,nopt+1,1.d0,bba,ndat,cc,1,0.d0,vfits,1)
#endif
#if _USE_MAGMA
call magma_dsetmatrix(ndat,nopt,bba,ndat,dbba,ndat,queue)
call magma_dsetvector(nopt,cc,1,dcc,1,queue)
call magma_dgemv(magmaNoTrans,ndat,nopt+1,1.d0,dbba,ndat,dcc,1,0.d0,dvfits,1,queue)
!call magma_dgetmatrix(ndat,nopt,dbba,ndat,bba,ndat,queue)
call magma_dgetvector(nopt,dvfits,1,vfits,1,queue)
#endif
call cpu_time(t212)
t21a=t212-t211
call cpu_time(t213)
do i=1,ndat
  errs(i)=vfits(i)-vv(i)
enddo
err=ddot(ndat,errs,1,errs,1)
wrr=ddot(ndat,ww,1,ww,1)
do i=1,ndat
  uvfits(i)=vfits(i)/ww(i)
  uvv(i)=vv(i)/ww(i)
enddo
iminvfit=minloc(uvfits,DIM=1)
iminvv=minloc(uvv,DIM=1)
vf1=uvv(iminvfit)
vf2=uvfits(iminvfit)
vv1=uvv(iminvv)
vv2=uvfits(iminvv)
call cpu_time(t214)
t21b=t214-t213

err=dsqrt(err/wrr)
if (err.lt.errmin) then
iisav=ii
errmin=err
endif

#else
vf2=1000.
vv1=1000.
wrr=0.d0
err=0.d0
t21a=0.d0 ! timing
do i=1,ndat
call cpu_time(t215)
vfit=0.d0
do k=1,nopt
vfit=vfit+bb(i,list(k))*cc(k) ! vfit is now weighted
enddo
vfit=vfit+bb(i,ii)*cc(nopt+1)
call cpu_time(t216)
t21a=t21a+t216-t215
call cpu_time(t217)
wrr=wrr+ww(i)**2
err=err+(vfit-vv(i))**2
if (vfit/ww(i).lt.vf2) then
vf1=vv(i)/ww(i)
vf2=vfit/ww(i)
endif
if (vv(i)/ww(i).lt.vv1) then
vv1=vv(i)/ww(i)
vv2=vfit/ww(i)
endif
enddo
err=dsqrt(err/wrr)
if (err.lt.errmin) then
iisav=ii
errmin=err
endif
call cpu_time(t218)
t21b=t218-t217

#endif
call cpu_time(t22)
!write(99,*)"Time spent in evaluating the candidate fit: ",t22-t21," s"
ttb=ttb+t22-t21
ttba=ttba+t21a
ttbb=ttbb+t21b
!write(99,*)"current tab == ",tab
enddo ! main loop end
!taa=tta/ncoef
!tab=ttb/ncoef

! permanently update the basis with the best one
call cpu_time(t23)
nopt=nopt+1 ! size of new basis
write(6,123)nopt,errmin,iisav,vv1,vv2,vf1,vf2
used(iisav)=1 ! flag the newly added one as used
list(nopt)=iisav ! add it to the list
call cpu_time(t231)
#if defined(_USE_LAPACK) || defined(_USE_LAPACK_TEST) || defined(_USE_MAGMA)
do k1=1,nopt
  do k2=1,nopt
    cov(k1,k2)=0.d0
  enddo
enddo
#endif

#if _USE_LAPACK
call dgemm('T','N',nopt,nopt,ndat,1.d0,bba,ndat,bba,ndat,1.d0,cov,ncoef) ! dgemm appears to be faster than dsyrk for larger matrices > 20000x100
#elif _USE_LAPACK_TEST
call dsyrk('u','T',nopt,ndat,1.d0,bba,ndat,0.d0,cov,ncoef)
do k1=2,nopt
  do k2=1,k1-1
    cov(k1,k2)=cov(k2,k1)
  enddo
enddo
#elif _USE_MAGMA
call magma_dsetmatrix(ndat,nopt,bba,ndat,dbba,ndat,queue)
call magma_dsetmatrix(nopt,nopt,cov,ncoef,dcov,ncoef,queue)
call magma_dgemm(MagmaTrans,MagmaNoTrans,nopt,nopt,ndat,1.d0,dbba,ndat,dbba,ndat,1.d0,dcov,ncoef,queue)
!call magma_dgetmatrix(ndat,nopt,dbba,ndat,bba,ndat,queue)
call magma_dgetmatrix(nopt,nopt,dcov,ncoef,cov,ncoef,queue)
#else
do k1=1,nopt ! cov needs to be re-evaluated due to the numerical instability of covu
do k2=1,nopt
cov(k1,k2)=0.d0
do i=1,ndat
cov(k1,k2)=cov(k1,k2)+bb(i,list(k1))*bb(i,list(k2))  ! X^T*X
enddo
enddo
enddo
#endif
call cpu_time(t232)
lwork=nopt
call cpu_time(t24)
#if _USE_MAGMA
call cpu_time(t241)
call magma_dgetrf(nopt,nopt,cov,ncoef,ipiv,info)
call cpu_time(t242)
if (info.ne.0) then
    write(6,*)"magma_dgetrf error: info ==",info
endif
!call magma_dsetmatrix(nopt,nopt,cov,ncoef,dcov,ncoef,queue)
!call magma_dgetri_gpu(nopt,dcov,ncoef,ipiv,dwork,lwork,info)
!call magma_dgetmatrix(nopt,nopt,dcov,ncoef,cov,ncoef,queue)
!if (info.ne.0) then
!    write(6,*)"magma_dgetri error: info ==",info
!endif
call dgetri(nopt,cov,ncoef,ipiv,work,lwork,info)
#else
call cpu_time(t241)
call dgetrf(nopt,nopt,cov,ncoef,ipiv,info)           ! invert it using LAPACK
call cpu_time(t242)
call dgetri(nopt,cov,ncoef,ipiv,work,lwork,info)
#endif
call cpu_time(t25)

call cpu_time(tt2)
!write(99,*)"Avg. time updating cov (evaluating covu): ",taa," s"
!write(99,*)"Avg. time evaluating the candidate fit: ",tab," s"
write(99,*)"Time spent in updating cov (ncoef-nopt times): ",tta," s"
write(99,*)"Time spent in evaluating the candidate fit (ncoef-nopt times): ",ttb," s"
write(99,*)"    Time spent in evaluating X*c (ncoef-nopt times): ",ttba," s"
write(99,*)"    Time spent in analysis of the candidate fit (ncoef-nopt times): ",ttbb," s"
write(99,*)"Time spent in updating the basis (including a re-evaluation of cov): ",t25-t23," s"
write(99,*)"Time spent in X^T * X: ",t232-t231," s"
write(99,*)"Total time adding a coefficient (taa+tab)*ncoef + updating basis: ",tt2-tt1," s"
write(99,*)"Time spent in dgetrf and dgetrt",t25-t24," s"
#if _USE_MAGMA
write(99,*)"Time spent in magma_dgetrf",t242-t241," s"
#else
write(99,*)"Time spent in dgetrf",t242-t241," s"
#endif
write(99,*)""
! decide when to end
if (nopt.eq.ntarget) go to 999
go to 100 ! repeat the cycle

123 format(i10,f12.4,i10,2(2x,2f12.4))

999 continue

#if _USE_MAGMA
write(6,*)"--------------- destroy queue"
call magma_queue_destroy(queue)

write(6,*)"--------------- finalize"
call magma_finalize()
write(6,*)"done"
#endif

close(99)
end

!subroutine print_matrix( A )
!    double precision :: A(:,:)
    
!    integer :: i
!    do i = 1, ubound(A,1)
!        write(6,*)(A(i,j),j=1,ubound(A,2))
!    enddo
!    write(6,*)
!end subroutine print_matrix
