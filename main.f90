!new scheme
!reference for program structure  
!Thor Gjesdal. "Programming Multigrid in F90", 1993

!>Module levels contains grid data type (data structure) definition,
!!the number of levels, and the array of grids with the highest level being the
!!finest grid, and the level 0 the coarsest
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                       !
!                                                       !
! Here, the algorithm ends and the data output begins   !
!                                                       !
!                                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
include "Subr_MG.f90"

program main
use levels
real, allocatable, dimension (:,:,:) :: inF, inA
real, allocatable, dimension (:,:) ::   inDivF, outSol
integer :: i, n, nx, j, k, l, param1, param2
real :: dx, mean, var_x, var_y, length
real(8),  parameter :: PI_8  = 4 * atan (1.0_8)

nx=8
length=pi_8
dx=length/nx

allocate(inF(1:nx+1,1:nx+1,0:1),inA(1:nx+1,1:nx+1,0:1),outSol(1:nx,1:nx))
outSol=0
inF=0
inA=1
var_x=1.
var_y=var_x
mean=pi_8/2
!...instead of having if statements inside, put them between the first loop and the second
!...like in the forcing interpolation
do k=0,1
do j=1,nx+1
do i=1,nx+1
	if (k==0) then
		x=dx*(real(i-1))
		y=dx*(real(j-1)+0.5)
	else 
		x=dx*(real(i-1)+0.5)
		y=dx*(real(j-1))
	endif
inA(i,j,k)=1!-exp(-0.5*( ((x-mean)/var_x)**2+ ((y-mean)/var_y)**2 ) ) / (var_x*var_y*2*pi_8) 
inF(i,j,k)=sin(x)+sin(y)
enddo
enddo
enddo

!!!!first argument is the number of iterations per level in the V-cycle
call MG(100,nx,1,1,inF,inA,outSol)
do i=level(number_of_levels -1)%nx,1,-1
print *, outSol(i,:)
enddo
end program
