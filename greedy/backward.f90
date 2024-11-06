program greedy

implicit double precision(a-h,o-z)
character*2 dum

! PIPPy is used to write the B matrix as defined in eq 2
! This code assumes the B matrix and the training data and weights are in xmat.dat
! This code also reads basis.dat to get the size of the full basis, ncoef
! The termination condition is hard coded such that it runs until it uses 1/2 the basis.

! set these
logical :: ldofull = .true.   ! Compute error for the full B matrix?
logical :: lreadin = .true.    ! Start with basis in greedy.dat?
integer :: ntarget = 50       ! End when you reach ntarget basis functions

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
double precision, dimension(:,:), allocatable :: covu ! temp cov. matrices(k1,k2)

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
print *,"Usage: backward.x LDOFULL NTARGET LREADIN"
stop
endif

! read the basis
open(10,file="basis.dat")
read(10,*)natom,npair,ncoef,nterm
print *,"From basis.dat: natom,npair,ncoef,nterm = ",natom,npair,ncoef,nterm
close(10)
if (ntarget.eq.-1) ntarget=ncoef

! read the B matrix (maybe change to binary, as this matrix can be large)
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
read(20,*)
i=0
do
  i=i+1
  read(20,*,end = 21)ii,ww(i),vv(i),(bb(i,k),k=1,ncoef)
end do
21 continue
close(20)

! use weights
do i=1,ndat
vv(i)=vv(i)*ww(i)  ! y = W*V (eq 3)
do k=1,ncoef
bb(i,k)=bb(i,k)*ww(i) ! X = W*B (eq 4)
enddo
enddo

allocate(cov(ncoef,ncoef))
allocate(ipiv(ncoef))
allocate(work(ncoef))
allocate(xty(ncoef))
allocate(cc(ncoef))

! solve for the full B matrix
if (ldofull) then
! compute (X^T*X)^-1 = cov(k1,k2)
lwork=ncoef
do k1=1,ncoef
do k2=1,ncoef
cov(k1,k2)=0.d0
do i=1,ndat
cov(k1,k2)=cov(k1,k2)+bb(i,k1)*bb(i,k2)  ! X^T*X
enddo
enddo
enddo
call dgetrf(ncoef,ncoef,cov,ncoef,ipiv,info)       ! together, these two LAPACK routines
call dgetri(ncoef,cov,ncoef,ipiv,work,lwork,info)  ! return the inverse of X^T*X
! now, cov = (X^T*X)^-1, i.e., the first part of eq 5
! next we compute xty = X^T*y, i.e., the rest of eq 5
do k=1,ncoef
xty(k)=0.d0
do i=1,ndat
xty(k)=xty(k)+bb(i,k)*vv(i) ! X^T*y in eq 5
enddo
enddo

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
else
do i=2,ncoef ! DOEST WORK
nopt=nopt+1
list(nopt)=nopt
enddo
endif
nread=nopt
!nopt=1

if (nopt.le.ntarget) then
print *,"NOPT <= NTARGET (",nopt,"<=",ntarget,") so there's nothing to do!"
stop
endif


! evaluate the fit as above...
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





! forward greedy algorithm used to read in basis
! this isn't used anymore in backward.f90 (nopt=nread so it goes to 999 right away)!!
! i'm leaving it in case we want to have a forward-backward code someday
allocate(u1(ncoef))
allocate(u2(ncoef))
allocate(u3(ncoef))
allocate(covu(ncoef,ncoef))

100 continue
errmin=10000.
do ii=1,ncoef ! main loop
if (nopt.eq.nread) go to 999
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

! test fit
vf2=1000.
vv1=1000.
wrr=0.d0
err=0.d0
do i=1,ndat
vfit=0.d0
do k=1,nopt
vfit=vfit+bb(i,list(k))/ww(i)*cc(k)
enddo
vfit=vfit+bb(i,ii)/ww(i)*cc(nopt+1)
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
if (err.lt.errmin) then
iisav=ii
errmin=err
endif
enddo

! permanenly update the basis wit the best one
nopt=nopt+1 ! size of new basis
write(6,123)nopt,errmin,iisav,vv1,vv2,vf1,vf2
used(iisav)=1 ! flag the newly added one as used
list(nopt)=iisav ! add it to the list
do k1=1,nopt
do k2=1,nopt
cov(k1,k2)=0.d0
do i=1,ndat
cov(k1,k2)=cov(k1,k2)+bb(i,list(k1))*bb(i,list(k2))  ! X^T*X
enddo
enddo
enddo
lwork=nopt
call dgetrf(nopt,nopt,cov,ncoef,ipiv,info)           ! invert it using LAPACK
call dgetri(nopt,cov,ncoef,ipiv,work,lwork,info)

! finished reading the basis
if (nopt.eq.nread) go to 999
go to 100 ! repeat the cycle

999 continue






! backward greedy algorithm
200 continue
errmin=10000.
do ii=2,nopt ! main loop

kk1=0 ! copy cov matrix but skip iith basis
do k1=1,nopt
if (k1.ne.ii) then
kk1=kk1+1
kk2=0
do k2=1,nopt
if (k2.ne.ii) then
kk2=kk2+1
covu(kk1,kk2)=cov(k1,k2)
endif
enddo
endif
enddo

d0 = cov(ii,ii)

kk=0
do k=1,nopt
if (k.ne.ii) then
kk=kk+1
u3(kk)=-cov(k,ii)
u2(kk)=-cov(k,ii)/d0
endif
enddo

do k1=1,nopt-1
do k2=1,nopt-1
covu(k1,k2)=covu(k1,k2)-d0*u2(k1)*u2(k2)
enddo
enddo

do k1=1,nopt-1
cc(k1)=0.d0
kk2=0
do k2=1,nopt-1
kk2=kk2+1
if (k2.eq.ii) kk2=kk2+1
cc(k1)=cc(k1)+covu(k1,k2)*xty(list(kk2))
enddo
enddo

! test fit
vf2=1000.
vv1=1000.
wrr=0.d0
err=0.d0
do i=1,ndat
vfit=0.d0
kk=0
do k=1,nopt
if (k.ne.ii) then
kk=kk+1
vfit=vfit+bb(i,list(k))/ww(i)*cc(kk)
endif
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
if (err.lt.errmin) then
iisav=ii
errmin=err
endif
enddo

! permanenly update the basis wit the best one
nopt=nopt-1 ! size of new basis
write(6,123)nopt,errmin,-list(iisav),vv1,vv2,vf1,vf2
kk=0
do k=1,nopt
kk=kk+1
if (k.eq.iisav) kk=kk+1
list(k)=list(kk)
enddo
do k1=1,nopt
do k2=1,nopt
cov(k1,k2)=0.d0
do i=1,ndat
cov(k1,k2)=cov(k1,k2)+bb(i,list(k1))*bb(i,list(k2))  ! X^T*X
enddo
enddo
enddo
lwork=nopt
call dgetrf(nopt,nopt,cov,ncoef,ipiv,info)           ! invert it using LAPACK
call dgetri(nopt,cov,ncoef,ipiv,work,lwork,info)

! finished reading the basis
if (nopt.eq.ntarget) go to 990
go to 200 ! repeat the cycle


990 continue










123 format(i10,f12.4,i10,2(2x,2f12.4))
end


