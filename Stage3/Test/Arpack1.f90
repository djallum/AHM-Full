subroutine Arpack1(dim,matrix,eigvalues,eigvectors)

  implicit none

  integer, intent(in) :: dim
  real, intent(in) :: matrix(dim,dim)
  real, intent(out) :: eigvalues(dim)
  real, intent(out) :: eigvectors(dim,dim)

  integer :: IDO, iparam(11), ipntr(11), n, nev
  integer :: ncv, lworkl, info, ierr, j, nx
  integer :: ishfts, maxitr, model, nconv
  double precision :: tol, sigma
  double precision, allocatable :: v(:,:), workl(:), ax(:)
  double precision, allocatable :: workd(:), d(:), resid(:), select(:)
  character :: bmat*1, which*2

  bmat = 'I'
  which = 'SA'
  nx = dim
  n = nx * nx
  nev = nx
  ncv = 2*nev
  lworkl = ncv*(ncv+8)

  allocate(v(nx,ncv),workd(3*nx),workl(lworkl),resid(nx),sel(ncv))
  
  tol = 0d.0
  info = 0
  IDO = 0

  ishfts = 1
  maxitr = 300
  model = 1

  iparam(1) = ishfts
  iparam(3) = maxitr
  iparam(7) = model

  call dsaupd( ido, bmat, nx, which, nev, tol, resid, ncv, v, nx, &
       iparam, ipntr, workd, workl, lworkl, info )
  
  do while ( (ido .eq. -1) .or. (ido .eq. 1) ) then

     call MaxtrixProduct!I believe I need a routine to calculate the matrix product of the array 'matrix' (input) with workd(ipntr(1)) and output that to workd(ipntr(2))

     call dsaupd( ido, bmat, nx, which, nev, tol, resid, ncv, v, nx, &
       iparam, ipntr, workd, workl, lworkl, info ) 
  end do

  if ( info .lt. 0 ) then
     print*, ""
     print*, "Error in ARPack, info: ", info
     !may need more errors
     stop
  end if

  allocate(d(nev))
  call dseupd( .True. , 'A', sel, d, v, nx, sigma, bmat, nx, which, nev, &
       tol, resid, ncv, v, nx, iparam, ipntr, workd, workl, lworkl, info )  ! This routine takes what has been calculated in dsaupd and gets the eigenvalues and vectors

  if ( info .lt. 0 ) then
     print*, ""
     print*, "Error in ARPack, info: ", info
     !may need to print more
     stop
  end if
  

  ! May need to add input of a logical to determine if eigenvectors will be computed, will be first arguement of dseupd

end subroutine Arpack1

subroutine MatrixProduct( dim, A, u, w) !Computes Au = w where A is a dim x dim matrix and u is a vector of length dim (similarly for w)
  implicit none
  integer, intent(in) :: dim
  real, intent(in) :: A(dim,dim), u(dim)
  real, intent(out) :: w(dim)
  integer :: i, j
  real :: sum

  
  do i=1,dim
     sum = 0
     
     do j=1,dim
        if ( A(i,j) .eq. 0 ) CYCLE
        sum = sum + A(i,j)*u(j)
     end do
     
     w(i) = sum

  end do

end subroutine MatrixProduct

  
  

  
