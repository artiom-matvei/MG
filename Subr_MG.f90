!new scheme
!reference for program structure  
!Thor Gjesdal. "Programming Multigrid in F90", 1993

!>Module levels contains grid data type (data structure) definition,
!!the number of levels, and the array of grids with the highest level being the
!!finest grid, and the level 0 the coarsest
module Levels

	type :: grid
                !>Number of points
		integer :: nx
                !>Distance between two points in the x (or y) direction
		real :: dx
                !>The coefficient in front of v
		real, pointer, dimension (:,:,:) :: a
                !>The given function f, since the RHS is div (f)
                real, pointer, dimension (:,:,:) :: F
                !>The approximation of the unknown function
		real, pointer, dimension (:,:) :: v
                !>The RHS of the equation
                real, pointer, dimension (:,:) :: divF
	end type grid

        !>The number of levels in the multi-grid algorithm
	integer :: number_of_levels

        !>Array of grids named level
	type (grid), dimension(:), allocatable :: level

end module Levels

!>Allocate memory for the grids by first determining the total number of gr and  1-D 
subroutine Allocate_grids ( nx_fine, length ) !Set_up in Thor's
	use levels
	integer :: nx_fine !nx should be a power of two!
	real :: length
	number_of_levels = No_levels(nx_fine)
	allocate ( level(0:number_of_levels-1) )
	nx=nx_fine*2
	do i=number_of_levels-1,0,-1
		nx=nx/2
		level(i)%nx=nx
		level(i)%dx=length/nx
		allocate (level(i)%a(1:nx+1,1:nx+1,0:1))
		allocate (level(i)%F(1:nx+1,1:nx+1,0:1))
		allocate (level(i)%divF(1:nx,1:nx))
		allocate (level(i)%v(0:nx+1,0:nx+1))
                level(i)%v=0
		!allocate (level(i)%v_0(-1:nx,-1:nx))
		!allocate (level(i)%v_1(-1:nx,-1:nx))
	enddo	
	contains
		function No_levels(n) result(ans)
		
		integer :: n, ans, i, nn
                nn=n; i=1
                do 
                  if (mod(nn,2)/= 0)  exit
                  if (nn/2<4)         exit
                  nn=nn/2; i=i+1
                enddo

		ans=i
		!print *, ans
		end function No_levels
		
end subroutine Allocate_grids

subroutine Initialize_forcings (pi_8,nx,param1,param2,inF,inA)
	use Levels
	
	integer :: i, n, nx, j, k, l, param1, param2
	real :: dx, mean, var_x, var_y
        real(8), intent(in) :: pi_8
        real, dimension (1:nx+1,1:nx+1,0:1) :: inF, inA

        !...if condition for given forcings begins here
	!...initialize finest grid forcings
	n = number_of_levels-1; nx=level(n)%nx	
	dx=level(n)%dx

        level(n)%F=inF
        level(n)%a=inA

	!...TODO
	!...add divF computation!
        do j=1,nx
        do i=1,nx
        level(n)%divF(i,j)=(level(n)%F(i+1,j,0)-level(n)%F(i,j,0)+&
          level(n)%F(i,j+1,1)-level(n)%F(i,j,1))/dx
        enddo
        enddo

	!...interpolate coarser grid forcings
	do l=n-1,0,-1
        dx=level(l)%dx
	nx=level(l)%nx
	!do k=0, left and right boundaries
	do j=1,nx
        !j=nx+1 is not interpolated because it is outside the grid and not used by the algorithm, 
        !the same is true for i=nx+1 for k=1
        i=1
        level(l)%F(i,j,0)=(level(l+1)%F(2*i-1,2*j,1)+&
          level(l+1)%F(2*i-1,2*j,0)+&
          level(l+1)%F(2*i-1,2*j-1,0))/3
	level(l)%a(i,j,0)=(level(l+1)%a(2*i-1,2*j,1)+&
          level(l+1)%a(2*i-1,2*j,0)+&
          level(l+1)%a(2*i-1,2*j-1,0))/3
	
	do i=2,nx !interpolate left and right boundaries separately
	level(l)%F(i,j,0)=(level(l+1)%F(2*i-1,2*j,1)+&
          level(l+1)%F(2*i-2,2*j,1)+&
          level(l+1)%F(2*i-1,2*j,0)+&
          level(l+1)%F(2*i-1,2*j-1,0))/4
	level(l)%a(i,j,0)=(level(l+1)%a(2*i-1,2*j,1)+&
          level(l+1)%a(2*i-2,2*j,1)+&
          level(l+1)%a(2*i-1,2*j,0)+&
          level(l+1)%a(2*i-1,2*j-1,0))/4
	enddo !i=2,nx
        i=nx+1
	level(l)%F(i,j,0)=(level(l+1)%F(2*i-2,2*j,1)+&
          level(l+1)%F(2*i-1,2*j,0)+&
          level(l+1)%F(2*i-1,2*j-1,0))/3
	level(l)%a(i,j,0)=(level(l+1)%a(2*i-2,2*j,1)+&
          level(l+1)%a(2*i-1,2*j,0)+&
          level(l+1)%a(2*i-1,2*j-1,0))/3
	enddo!j=1,nx

        !do k=1, top and bottom boundaries
       	do i=1,nx
        !i=nx+1 is not interpolated because it is outside the grid and not used by the algorithm, 
        j=1
        level(l)%F(i,j,1)=(level(l+1)%F(2*i,2*j-1,1)+&
          level(l+1)%F(2*i-1,2*j-1,1)+&
          level(l+1)%F(2*i,2*j-1,0))/3
	level(l)%a(i,j,1)=(level(l+1)%a(2*i,2*j-1,1)+&
          level(l+1)%a(2*i-1,2*j-1,1)+&
          level(l+1)%a(2*i,2*j-1,0))/3
        do j=2,nx
	level(l)%F(i,j,1)=(level(l+1)%F(2*i,2*j-1,1)+&
          level(l+1)%F(2*i-1,2*j-1,1)+&
          level(l+1)%F(2*i,2*j-1,0)+&
          level(l+1)%F(2*i,2*j-2,0))/4
	level(l)%a(i,j,1)=(level(l+1)%a(2*i,2*j-1,1)+&
          level(l+1)%a(2*i-1,2*j-1,1)+&
          level(l+1)%a(2*i,2*j-1,0)+&
          level(l+1)%a(2*i,2*j-2,0))/4
        enddo
        j=nx+1
        level(l)%F(i,j,1)=(level(l+1)%F(2*i,2*j-1,1)+&
          level(l+1)%F(2*i-1,2*j-1,1)+&
          level(l+1)%F(2*i,2*j-2,0))/3
	level(l)%a(i,j,1)=(level(l+1)%a(2*i,2*j-1,1)+&
          level(l+1)%a(2*i-1,2*j-1,1)+&
          level(l+1)%a(2*i,2*j-2,0))/3
	enddo !i=1,nx

      !TODO
      !probably the divF for low levels isn't even needed...
      !delete or comment out this code 
        do j=1,nx
        do i=1,nx
        level(l)%divF(i,j)=(level(l)%F(i+1,j,0)-level(l)%F(i,j,0)+&
          level(l)%F(i,j+1,1)-level(l)%F(i,j,1))/dx
        enddo
        enddo
	enddo! l=nb_of_levels-2,0,-1
end subroutine Initialize_forcings

subroutine Smoothing(num_its, k)
use levels
integer nx, num_its, l, i, j, k 
real, dimension (:,:), pointer :: divF, v_0, v_1
real, dimension (:,:,:), pointer :: a 
real :: dx, mean
a=>level(k)%a
divF=>level(k)%divF
v_0=>level(k)%v
dx=level(k)%dx
nx=level(k)%nx


do l=0,num_its 
	
        mean=0

	v_0(0,1:nx)=v_0(1,1:nx)
	v_0(nx+1,1:nx)=v_0(nx,1:nx)
	v_0(1:nx,0)=v_0(1:nx,1)
	v_0(1:nx,nx+1)=v_0(1:nx,nx)
	do j=1,nx
	do i=1,nx

	v_0(i,j)=(v_0(i+1,j)*a(i+1,j,0)+v_0(i-1,j)*a(i,j,0)+v_0(i,j+1)*a(i,j+1,1)+v_0(i,j-1)*a(i,j,1)&
			&-dx**2*divF(i,j))/(a(i+1,j,0)+a(i,j,0)+a(i,j+1,1)+a(i,j,1))
        mean=mean+v_0(i,j)/nx**2
	enddo
	enddo

	do j=1,nx
	do i=1,nx
	v_0(i,j)=v_0(i,j)-mean
	enddo
	enddo

enddo

end subroutine Smoothing

subroutine Res_comp_and_coars(k)
use levels
integer nx, num_its, i, j, k 
real, dimension (:,:), pointer :: divF, v_0, v_1
real, dimension (:,:,:), pointer :: a 
real :: dx
a=>level(k)%a
divF=>level(k)%divF
v_0=>level(k)%v
dx=level(k)%dx
nx=level(k)%nx

!here v_1 is used to compute divF - div (a*grad u)
allocate(v_1(1:nx,1:nx))

do j=1,nx
do i=1,nx
	v_1(i,j)=divF(i,j)-(v_0(i+1,j)*a(i+1,j,0)+v_0(i-1,j)*a(i,j,0)+v_0(i,j+1)*a(i,j+1,1)+v_0(i,j-1)*a(i,j,1)&
			&-v_0(i,j)*(a(i+1,j,0)+a(i,j,0)+a(i,j+1,1)+a(i,j,1)))/dx**2
enddo
enddo

!...continue here by passing the residual to the next lvl and coarsening it
do j=1,level(k-1)%nx
do i=1,level(k-1)%nx
level(k-1)%divF(i,j)=(v_1(2*i,2*j)+v_1(2*i-1,2*j)+v_1(2*i,2*j-1)+v_1(2*i-1,2*j-1))/4
enddo
enddo

end subroutine Res_comp_and_coars

subroutine Interp_er_and_cor(k)
!k is the level where the current error is situated        
use levels
integer nx, num_its, i, j, k 
real, dimension (:,:), pointer :: divF, v, v_1
real, dimension (:,:,:), pointer :: a 
real :: dx
a=>level(k)%a
divF=>level(k)%divF
v=>level(k)%v
dx=level(k)%dx
nx=level(k)%nx

!...TODO
!...VERY IMPORTANT
!...verify somehow if I need to add or substract the error...!!!
!...change the ghost points adjustment to take place before every smoothing iteration, not after, to ensure everything is set up

!interpolate the boundaries separately, go counterclockwise from bottom left "coin"
!or go in the order as if inside the loop
i=1; j=1
level(k+1)%v(i,j)=level(k+1)%v(i,j)+v(i,j)
do j=1,nx-1
level(k+1)%v(i,2*j)=level(k+1)%v(i,2*j)        +(3*v(i,j)+1*v(i,j+1))/4
level(k+1)%v(i,2*j+1)=level(k+1)%v(i,2*j+1)    +(1*v(i,j)+3*v(i,j+1))/4
enddo
i=1; j=nx
level(k+1)%v(i,2*j)=level(k+1)%v(i,2*j)+v(i,j)
do i=1,nx-1
level(k+1)%v(2*i,2*j)=level(k+1)%v(2*i,2*j)        +(3*v(i,j)+1*v(i+1,j))/4
level(k+1)%v(2*i+1,2*j)=level(k+1)%v(2*i+1,2*j)    +(1*v(i,j)+3*v(i+1,j))/4
enddo
i=nx; j=nx
level(k+1)%v(2*i,2*j)=level(k+1)%v(2*i,2*j)        +v(i,j)
do j=1,nx-1
level(k+1)%v(2*i,2*j)=level(k+1)%v(2*i,2*j)        +(3*v(i,j)+v(i,j+1))/4
level(k+1)%v(2*i,2*j+1)=level(k+1)%v(2*i,2*j+1)    +(v(i,j)+3*v(i,j+1))/4
enddo
i=nx; j=1
level(k+1)%v(2*i,j)=level(k+1)%v(2*i,j)+v(i,j)
do i=1, nx-1
level(k+1)%v(2*i,j)=level(k+1)%v(2*i,j)        +(3*v(i,j)+1*v(i+1,j))/4
level(k+1)%v(2*i+1,j)=level(k+1)%v(2*i+1,j)    +(1*v(i,j)+3*v(i+1,j))/4
enddo

do j=1,nx-1
do i=1,nx-1
level(k+1)%v(2*i,2*j)=level(k+1)%v(2*i,2*j)        +(9*v(i,j)+3*v(i+1,j)+3*v(i,j+1)+1*v(i+1,j+1))/16
level(k+1)%v(2*i,2*j+1)=level(k+1)%v(2*i,2*j+1)    +(3*v(i,j)+1*v(i+1,j)+9*v(i,j+1)+3*v(i+1,j+1))/16
level(k+1)%v(2*i+1,2*j)=level(k+1)%v(2*i+1,2*j)    +(3*v(i,j)+9*v(i+1,j)+1*v(i,j+1)+3*v(i+1,j+1))/16
level(k+1)%v(2*i+1,2*j+1)=level(k+1)%v(2*i+1,2*j+1)+(1*v(i,j)+3*v(i+1,j)+3*v(i,j+1)+9*v(i+1,j+1))/16
enddo
enddo

end subroutine Interp_er_and_cor

subroutine MG(num_its,nx,param1,param2,inF, inA, outSol)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !This program needs two parameters when called! ! !The first is the number of iterations,         !  the second is the number of grid points        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use levels
integer nx, num_its, l, i, j, k, param1, param2

real, dimension (1:nx+1,1:nx+1,0:1) :: inF, inA
real, dimension (1:nx,1:nx) :: outSol

real, dimension (:,:), pointer :: divF, v_0, v_1
real, dimension (:,:,:), pointer :: a, F 
real length, dx, mean, var_x, var_y, x, y
real(8),  parameter :: PI_8  = 4 * atan (1.0_8)
character(100) :: num1char, num2char!, num3char
!allocate(v_1(1:nx,1:nx)) !call get_command_argument(1,num1char) !call get_command_argument(1,num1char) !call get_command_argument(2,num2char) !read(num1char,*)num_its !read(num2char,*)nx
length=pi_8
!allocate (v_0(-1:nx,-1:nx), v_1(-1:nx,-1:nx), divF(0:nx-1,0:nx-1), a(0:nx,0:nx,0:1), F(0:nx,0:nx,0:1))

call Allocate_grids(nx, length)
call Initialize_forcings(pi_8,nx,param1,param2,inF, inA)

do k=number_of_levels-1,1,-1
call Smoothing (num_its, k)
call Res_comp_and_coars(k)
enddo

!print *, level(number_of_levels-1)%v

do k=0,number_of_levels-2
call Smoothing (num_its,k)
call Interp_er_and_cor(k)
enddo

k=number_of_levels-1
call Smoothing (num_its,k)

nx=level(number_of_levels -1)%nx
outSol=level(k)%v(1:nx,1:nx)
end subroutine MG
