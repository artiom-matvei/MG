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
program main
include Subr_MG.f90
use levels
F=>level(k)%F
a=>level(k)%a
divF=>level(k)%divF
v_0=>level(k)%v
dx=level(k)%dx
nx=level(k)%nx
allocate(v_1(1:nx,1:nx))
do j=1,nx
do i=1,nx
	v_1(i,j)=divF(i,j)-(v_0(i+1,j)*a(i+1,j,0)+v_0(i-1,j)*a(i,j,0)+v_0(i,j+1)*a(i,j+1,1)+v_0(i,j-1)*a(i,j,1)&
			&-v_0(i,j)*(a(i+1,j,0)+a(i,j,0)+a(i,j+1,1)+a(i,j,1)))/dx**2
enddo
enddo
do i=nx,1,-1
print *, v_0(1:nx,i)
enddo

!print *, "The above table is the solution, the below is the residual." 

do j=nx,1,-1
!print *, v_1(:,j)
enddo
!print *, "Below is the divF"
do j=nx,1,-1
        print *, v_1(:,j)
enddo
do j=level(k-1)%nx+1,1,-1
        !print *, level(k-1)%a(:,j,0)
	enddo

open(20, file = 'error_gnu', status = 'unknown')
write(20,*) 'set cbrange[-1:1] #use to change coloarbar' 
write(20,*) 'set zrange[-1:1] #use to whiteout extreme values'
write(20,*) '# indicates comment'
write(20,*) '# the first two lines are optional'
write(20,*) 'set pm3d map'
write(20,*) 'splot "-"'

do j = 1,nx
	y =  float(j)*dx-dx/2 !zero to one
	do i = 1,nx
		x =  float(i)*dx-dx/2
		write(20,*) x, y, v_1(i,j)
	enddo 
	write(20,*) ! write blank line between rows 
enddo  
close(20)

open(20, file = 'h_gnu', status = 'unknown')!a(i,j,0)
write(20,*) 'set cbrange[-1:1] #use to change coloarbar' 
write(20,*) 'set zrange[-1:1] #use to whiteout extreme values'
write(20,*) '# indicates comment'
write(20,*) '# the first two lines are optional'
write(20,*) 'set pm3d map'
write(20,*) 'splot "-"'

do j = 1,nx
	y =  float(j)*dx  !zero to one
	do i = 1,nx+1
		x =  float(i)*dx
		write(20,*) x, y, a(i,j,0)
	enddo 
	write(20,*) ! write blank line between rows 
enddo  
close(20)

open(20, file = 'rhs_gnu', status = 'unknown')
write(20,*) 'set cbrange[-1:1] #use to change coloarbar' 
write(20,*) 'set zrange[-1:1] #use to whiteout extreme values'
write(20,*) '# indicates comment'
write(20,*) '# the first two lines are optional'
write(20,*) 'set pm3d map'
write(20,*) 'splot "-"'

do j = 1,nx
	y =  float(j)*dx-dx/2 !zero to one
	do i = 1,nx
		x =  float(i)*dx-dx/2
		write(20,*) x, y, divF(i,j)
	enddo 
	write(20,*) ! write blank line between rows 
enddo  
close(20)



open(20, file = 'phi_gnu', status = 'unknown')
write(20,*) 'set cbrange[-1:1] #use to change coloarbar' 
write(20,*) 'set zrange[-1:1] #use to whiteout extreme values'
write(20,*) '# indicates comment'
write(20,*) '# the first two lines are optional'
write(20,*) 'set pm3d map'
write(20,*) 'splot "-"'

do j = 1,nx
	y =  float(j)*dx-dx/2 !zero to one
	do i = 1,nx
		x =  float(i)*dx-dx/2
		write(20,*) x, y, v_0(i,j)
	enddo 
	write(20,*) ! write blank line between rows 
enddo  
close(20)


open(20, file = 'fx_gnu', status = 'unknown')!F(i,j,0)
write(20,*) 'set cbrange[-1:1] #use to change coloarbar' 
write(20,*) 'set zrange[-1:1] #use to whiteout extreme values'
write(20,*) '# indicates comment'
write(20,*) '# the first two lines are optional'
write(20,*) 'set pm3d map'
write(20,*) 'splot "-"'

do j = 1,nx
	y =  float(j)*dx  !zero to one
	do i = 1,nx+1
		x =  float(i)*dx
		write(20,*) x, y, F(i,j,0)
	enddo 
	write(20,*) ! write blank line between rows 
enddo  
close(20)

end program
!tests to do
!0. verify against an analytic function
!1. loop at the coarsest level until convergence (stationary point)
!2. test against other solvers like:
!3. First grid configuration 
!4. new grid configuration
