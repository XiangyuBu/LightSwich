subroutine initialize()
USE global_parameters
USE utility_routines 
IMPLICIT none
INCLUDE 'mpif.h'

Integer*8 npts,mynpts
Integer :: flag_c, index_con, n_graft_point, x_flag, x_stat
Integer :: i, j, k, j_temp
integer :: n, m , l, length
integer ::  change,change_1  ! change is the flag of MC, 1 for accepte the move.
integer :: ddt(8) 

DOUBLE PRECISION :: axis(3)
DOUBLE PRECISION :: r, cos_t, sin_t, phi
DOUBLE PRECISION :: unew(3), uold(3) 
DOUBLE PRECISION :: alpha, beta, angle, dotp
DOUBLE PRECISION :: ratio_sigma, a_plane,loa_azo 
DOUBLE PRECISION, PARAMETER :: TOL = 1.0D-5 
DOUBLE PRECISION :: r_radius,r_radius_1, zz, rr, temp, x_r, y_r, temp_1


type(node) :: graft_point(0:401)  
character*7 res,resres,res_s
character res0,res1,res2
logical alive,check

call date_and_time(values=ddt)
seed=(ddt(8)-500)*654321*(myid+1) + 88888*myid
!print*, "myid",myid,"seed",seed
!seed=(ddt(8)-500)*54321 + 11223344
open(unit=15,file='log.txt')
open(unit=10,file='input.txt')

read(10,*) azo_position    ! azo position / L from the free end
read(10,*) ratio_sigma     ! grafting density ratio between substrate and sphere sigma_p/sigma_s
read(10,*) csoL            ! distance between surfaces of substrate and center of sphere / L
read(10,*) tau             ! end interaction parameter
read(10,*) nu              ! interaction parameter
read(10,*) Loa             ! L/a
read(10,*) roL             ! radius of sphere/L
read(10,*) N_chain         ! number of grafted chains
read(10,*) lambda          ! mixing parameter
read(10,*) Nm              ! number of bonds
read(10,*) Nm_pol          ! number of polymer's bonds
read(10,*) Nm_sub          ! number of substrate except the azo
read(10,*) index_con       ! number of conformations  
read(10,*) Nr              ! number of points in r direction
read(10,*) Nz              ! number of points in z direction
read(10,*) Npre            ! number of prerotated   
read(10,*) Nmove           ! number of move in a MC step 
read(10,*) Max_iter        ! Max of iterations
read(10,*) num             ! num==1 ,pivot move with no small pivot move 
read(10,*) rotate          ! the limitation of the random move angle 
read(10,*) rotate_s        ! the limitation of the rotate of sphere angle
read(10,*) move_max        ! max of sphere trial move
read(10,*) hahah           ! the real loa      !!!!!It will be change!!!!!
read(10,*) loa_azo         ! loa_azo
close(10)

allocate(polymer(1:N_chain,0:Nm_pol))
allocate(w(1:Nr+1,1:Nz),w_new(1:Nr+1,1:Nz),eta_new(1:Nr+1,1:Nz),eta(1:Nr+1,1:Nz))
allocate(eta_azo_new(1:Nr+1,1:Nz),eta_azo(1:Nr+1,1:Nz))
allocate(eta_sub_new(1:Nr+1,1:Nz),eta_sub(1:Nr+1,1:Nz))
allocate(ir(1:N_chain,0:Nm_pol),iz(1:N_chain,0:Nm_pol),bond_vector(1:N_chain))
allocate(r_a(1:Nr+1),rr_r(1:Nr+1),zz_r(1:Nz))
allocate(trial_move_rate(1:4))

trial_move_rate(1) = 0.6d0 ! rate of azo move to sphere move
trial_move_rate(2) = 0.4d0 ! rate of polymer move in sphere move
!trial_move_rate(3) = 1.0d0 ! rate interval of sphere rotate 
trial_move_rate(3) = 0.7d0 ! rate interval of sphere rotate
trial_move_rate(4) = 0.5d0 
!trial_move_rate(3) = 1.d0 for the position of sphere is fixed.


NMCs = 5*10**index_con



!!!!!!!!!!!!!!!!!!
! calculate rho_0
!!!!!!!!!!!!!!!!!!!!!
!rho_0 = 1.0d0*Nm * (4*n_azo + N_chain) / (4*Nz*Nr*Nr )
rho_0 = 1.0d0*(4*n_azo*nm + 4*n_sub*nm_sub + N_chain*Nm_pol) / (4*Nz*Nr*Nr )
write(15,*) "rho_0", rho_0

!!!!!!!!!!!!!!!!!!!!!
! the i_azo
!!!!!!!!!!!!!!!!!!!!!

! depend on Loa , hahah depend on the real loa
! compute bending coefficent of chain
epsilon_azo = 1.0d0*Nm/(4.0d0*loa_azo)       
epsilon = 1.0d0*Nm/(4.0d0*hahah)             
deltaS = 1.0d0*Loa/Nm                        


! depend on roL
! r used in MC
r_sphere = roL*Nm                                   
r_sphere_2 = r_sphere*r_sphere

Lr = Nm*(roL+1)
!Lz = Nm*(roL+1+csoL)  
Lz = Nm*(roL+8.25+csoL) 

dr = 1.0d0*Lr/Nr
dz = 1.0d0*Lz/Nz
Lbox = Lr

r_a(1)= 0.5d0*dr
do i= 2, nr + 1
    r_a(i) = r_a(i-1) + dr    
end do

! r used in physical world
r_dr = 1.0d0*Loa*Lr/Nr/Nm
r_dz = 1.0d0*Loa*Lz/Nz/Nm

rr_r(1)= 0.5d0*r_dr
do i= 2, nr + 1
    rr_r(i) = rr_r(i-1) + r_dr    
end do
zz_r(1)= r_dz
do i= 2, nz
    zz_r(i) = zz_r(i-1) + r_dz    
end do

!!!!!!!!!!!!
! initialize the omega
!!!!!!!!!!!!                 
w = 0
w_new = 0
eta = 0
eta_azo = 0
eta_sub = 0

inquire(file='wfield.txt',exist=alive) 
if(alive) then                          
    open(unit=42,file='wfield.txt',status='old')
    read(42,*)
    print*, "read omega"
    do j=1,Nz
        do i=1,Nr
            read(42,*) w(i,j), eta(i,j)
        enddo
    enddo
    close(42)
endif

!!!! initialize the chains on the substrate!!!!!!have some questions!!!!!!!!!
p_sphere_2 = csoL * Nm 
a_plane = dsqrt( 2.0d0*pi*r_sphere_2/(dsqrt(3.0d0)*ratio_sigma*N_chain) )             
write(15,*) "lattice constant on graftted surface",a_plane
graft_point(1)%x = a_plane*0.5d0
graft_point(1)%y = a_plane*dsqrt(3.0d0)*0.5d0
graft_point(1)%z = p_sphere_2
x_stat = 1
x_flag = 1
n_graft_point = 1
do while(x_flag == 1)
    n_graft_point = n_graft_point + 1
    do while(x_flag == 1)
        temp = graft_point(n_graft_point-1)%x + a_plane
        if(temp>Lbox)then
            x_stat = x_stat + 1
            exit
        end if
        graft_point(n_graft_point)%x = temp
        graft_point(n_graft_point)%y = graft_point(n_graft_point-1)%y
        graft_point(n_graft_point)%z = p_sphere_2
        temp_1 = (graft_point(n_graft_point)%x)**2 + (graft_point(n_graft_point)%y)**2 
        do while(temp_1 <= (2*r_sphere)**2)            
            graft_point(n_graft_point)%x = graft_point(n_graft_point)%x + a_plane                 
            temp_1 = (graft_point(n_graft_point)%x)**2 + (graft_point(n_graft_point)%y)**2
        end do        
        n_graft_point = n_graft_point + 1
    end do    
    temp = graft_point(n_graft_point-1)%y + a_plane*dsqrt(3.0d0)*0.5d0
    if(temp>Lbox)then
        exit
    end if  
    graft_point(n_graft_point)%x = a_plane*0.5d0 * mod(x_stat,2)
    graft_point(n_graft_point)%y = temp
    graft_point(n_graft_point)%z = p_sphere_2
    temp_1 = (graft_point(n_graft_point)%x)**2 + (graft_point(n_graft_point)%y)**2 
    do while(temp_1 <= (2*r_sphere)**2)
        graft_point(n_graft_point)%x = graft_point(n_graft_point)%x + a_plane       
        temp_1 = (graft_point(n_graft_point)%x)**2 + (graft_point(n_graft_point)%y)**2        
    end do    
end do
n_azo = n_graft_point - 2
!print*,  n_azo,"number of substrate-polymers and azo_polymers"
write(15,*) "number of substrate-polymers and azo_polymers",n_azo
 
open(unit=45,file='graft_point.txt')
do i = 1, n_azo   
    write(45,*) graft_point(i+1)%x, graft_point(i+1)%y, graft_point(i+1)%z
end do
close(45)
!print*, i_azo,"i_azo"

n_azo = floor(n_azo/2)
n_sub = n_azo
allocate(azo(1:n_azo,0:Nm), ir_azo(1:n_azo,0:Nm), iz_azo(1:n_azo,0:Nm), i_azo(1:n_azo))
allocate(sub(1:n_sub,0:Nm_sub), ir_sub(1:n_sub,0:Nm_sub), iz_sub(1:n_sub,0:Nm_sub))

i_azo_temp = floor(azo_position*Nm) !this is used in the file of utility ,line 551.

i_azo = Nm + 1
do i = 1, n_azo 
    i_azo(i) = floor(azo_position*Nm)     
end do

do j=1,n_azo 
    azo(j,0)%x = graft_point(2*j)%x 
    azo(j,0)%y = graft_point(2*j)%y
    azo(j,0)%z = graft_point(2*j)%z
	  do i=1,Nm
        change = 0
        if ( i==i_azo(j) ) then      
            unew(1) = azo(j,i-1)%x + 1.0d0
            unew(2) = azo(j,i-1)%y
            unew(3) = azo(j,i-1)%z
            r_radius = unew(1)*unew(1) + unew(2)*unew(2) + unew(3)*unew(3)
            if (r_sphere_2 > r_radius .or. unew(3) > p_sphere_2 ) then     ! jiao die le
                change = 0 
                stop "initial wrong"
            else
            	azo(j,i)%x = unew(1)
                azo(j,i)%y = unew(2)
                azo(j,i)%z = unew(3) 
            end if
        else if ( i==i_azo(j)+1 ) then  
!            print*,"i_azo+1"
            unew(1) = azo(j,i-1)%x - 0.5d0
            unew(2) = azo(j,i-1)%y + 0.5d0*dsqrt(3.0d0)
            unew(3) = azo(j,i-1)%z 
            r_radius = unew(1)*unew(1) + unew(2)*unew(2) + unew(3)*unew(3)
            if (r_sphere_2 > r_radius .or. unew(3) > p_sphere_2 ) then     ! jiao die le
                change = 0 
                stop "initial wrong"
            else
            	azo(j,i)%x = unew(1)
                azo(j,i)%y = unew(2)
                azo(j,i)%z = unew(3) 
            end if 
        else        
            do while(change == 0)                 ! make sure wall is impentrate
                change = 1        
                cos_t=(2*ran2(seed)-1)*0.99999999d0
                sin_t=dsqrt(1.0d0-cos_t**2) 

                phi = (2*ran2(seed) - 1)*pi      !generate angle of phi \in (-pi,+pi)
                axis(1) = sin_t*dcos(phi)
                axis(2) = sin_t*dsin(phi)
                axis(3) = cos_t
        		    unew(1) = azo(j,i-1)%x + axis(1)
                unew(2) = azo(j,i-1)%y + axis(2)
    	      	  unew(3) = azo(j,i-1)%z + axis(3) 
                r_radius = unew(1)*unew(1) + unew(2)*unew(2) + unew(3)*unew(3)
                if (r_sphere_2 > r_radius .or. unew(3) > p_sphere_2 ) then     ! jiao die le
                    change = 0 
                else
            	   azo(j,i)%x = unew(1)
                    azo(j,i)%y = unew(2)
                    azo(j,i)%z = unew(3) 
                    exit
                end if                
            end do !! do while
        end if    
    end do
end do

do j=1,n_sub 
    sub(j,0)%x = graft_point(2*j + 1)%x 
    sub(j,0)%y = graft_point(2*j + 1)%y
    sub(j,0)%z = graft_point(2*j + 1)%z
	do i=1,Nm_sub
        change = 0
        do while(change == 0)                 ! make sure wall is impentrate
            change = 1        
            cos_t=(2*ran2(seed)-1)*0.99999999d0
            sin_t=dsqrt(1.0d0-cos_t**2) 

            phi = (2*ran2(seed) - 1)*pi      !generate angle of phi \in (-pi,+pi)
            axis(1) = sin_t*dcos(phi)
            axis(2) = sin_t*dsin(phi)
            axis(3) = cos_t
        		unew(1) = sub(j,i-1)%x + axis(1)
            unew(2) = sub(j,i-1)%y + axis(2)
    	      unew(3) = sub(j,i-1)%z + axis(3) 
            r_radius = unew(1)*unew(1) + unew(2)*unew(2) + unew(3)*unew(3)
            if (r_sphere_2 > r_radius .or. unew(3) > p_sphere_2 ) then     ! jiao die le
                change = 0 
            else
                sub(j,i)%x = unew(1)
                sub(j,i)%y = unew(2)
                sub(j,i)%z = unew(3) 
                exit
            end if                
        end do !! do while    
	end do
end do

!read grafting points
open(unit=43,file='grafting_points.txt')
do j=1,N_chain
	read(43,*) bond_vector(j)%x,bond_vector(j)%y, bond_vector(j)%z
	r_radius = dsqrt( bond_vector(j)%x**2 + bond_vector(j)%y**2 + bond_vector(j)%z**2 )
    bond_vector(j)%x = bond_vector(j)%x / r_radius
    bond_vector(j)%y = bond_vector(j)%y / r_radius
    bond_vector(j)%Z = bond_vector(j)%z / r_radius
    polymer(j,0)%x = r_sphere * bond_vector(j)%x 
    polymer(j,0)%y = r_sphere * bond_vector(j)%y
    polymer(j,0)%z = r_sphere * bond_vector(j)%z
end do
close(43)

do j=1,N_chain
    do i=1,Nm_pol
        change = 0
		    do while(change == 0)                 ! make sure wall is impentrate
            change = 1        
            cos_t=(2*ran2(seed)-1)*0.99999999d0
            sin_t=dsqrt(1.0d0-cos_t**2) 

            phi = (2*ran2(seed) - 1)*pi      !generate angle of phi \in (-pi,+pi)
            axis(1) = sin_t*dcos(phi)
            axis(2) = sin_t*dsin(phi)
            axis(3) = cos_t
    		unew(1) = polymer(j,i-1)%x + axis(1)
			unew(2) = polymer(j,i-1)%y + axis(2)
    	   	unew(3) = polymer(j,i-1)%z + axis(3) 
            r_radius = unew(1)*unew(1) + unew(2)*unew(2) + unew(3)*unew(3)
            if (r_sphere_2 > r_radius .or. unew(3) > p_sphere_2 ) then     ! jiao die le
                change = 0 
            else
        		polymer(j,i)%x = unew(1)
                polymer(j,i)%y = unew(2)
        		polymer(j,i)%z = unew(3) 
                exit
            end if        
        
        end do !! do while    
    end do
end do


call checkpolymer (flag_c)



open(22,file='iriz.dat')
do j=1,N_chain
	do i=1,Nm_pol
	    r_radius = dsqrt( polymer(j,i)%x*polymer(j,i)%x + polymer(j,i)%y*polymer(j,i)%y )
    	ir(j,i) = floor( r_radius / dr ) + 1    
    	iz(j,i) = floor( ( p_sphere_2 - polymer(j,i)%z ) / dz ) + 1
		write(22,*) j,i, ir(j,i), iz(j,i)
	end do   
end do
close(22)
open(22,file='iriz_azo.dat')
do j=1,N_azo
	do i=1,Nm
        if ( azo(j,i)%x<=Lbox .and. azo(j,i)%y<=Lbox ) then
            x_r = azo(j,i)%x
            y_r = azo(j,i)%y             
        else if ( azo(j,i)%x>Lbox .and. azo(j,i)%y>Lbox ) then
            x_r = azo(j,i)%x - 2*Lbox
            y_r = azo(j,i)%y - 2*Lbox
        else if ( azo(j,i)%x>Lbox ) then
            x_r = azo(j,i)%x - 2*Lbox
            y_r = azo(j,i)%y
        else
            x_r = azo(j,i)%x
            y_r = azo(j,i)%y - 2*Lbox
        end if    
        r_radius = dsqrt( x_r*x_r + y_r*y_r )
        if(r_radius<Lr)then        
			     ir_azo(j,i) = floor( r_radius / dr ) + 1
        else
            ir_azo(j,i) = Nr + 1
		end if
    	iz_azo(j,i) = floor( ( p_sphere_2 - azo(j,i)%z ) / dz ) + 1		
		write(22,*) j,i, ir_azo(j,i), iz_azo(j,i)
	end do       
end do
close(22)
open(22,file='iriz_sub.dat')
do j=1,N_sub
	do i=1,Nm_sub
        if ( sub(j,i)%x<=Lbox .and. sub(j,i)%y<=Lbox ) then
            x_r = sub(j,i)%x
            y_r = sub(j,i)%y             
        else if ( sub(j,i)%x>Lbox .and. sub(j,i)%y>Lbox ) then
            x_r = sub(j,i)%x - 2*Lbox
            y_r = sub(j,i)%y - 2*Lbox
        else if ( sub(j,i)%x>Lbox ) then
            x_r = sub(j,i)%x - 2*Lbox
            y_r = sub(j,i)%y
        else
            x_r = sub(j,i)%x
            y_r = sub(j,i)%y - 2*Lbox
        end if    
        r_radius = dsqrt( x_r*x_r + y_r*y_r )
        if(r_radius<Lr)then        
			     ir_sub(j,i) = floor( r_radius / dr ) + 1
        else
            ir_sub(j,i) = Nr + 1
		end if
    	iz_sub(j,i) = floor( ( p_sphere_2 - sub(j,i)%z ) / dz ) + 1		
		write(22,*) j,i, ir_sub(j,i), iz_sub(j,i)
	end do       
end do
close(22)

!print*,"initial OK"

length = 1

j=0
do i=1,length
 
    call move_sphere(change)
    if (change==1) then
        j = j + 1
        call checkpolymer (flag_c)
        if (mod(j,50)==0) then
        write(15,*) j, p_sphere_2
    	end if
    end if
 
end do
call checkpolymer (flag_c)
!print*, 1.0d0*j/length,"sphere moved"
 


j=0
do i=1,length
 
    call pivot_azo(change)
    if (change==1) then
        j = j + 1
    end if
end do
call checkpolymer (flag_c)
!print*, 1.0d0*j/length,"pivot azo"


j=0
do i=1,length
 
    call pivot_sub(change)
    if (change==1) then
        j = j + 1
    end if
end do
call checkpolymer (flag_c)
!print*, 1.0d0*j/length,"pivot sub"


 
j=0
do i=1,length
 
    call pivot(change)
    if (change==1) then
        j = j + 1
    end if
end do
call checkpolymer (flag_c)
!print*, 1.0d0*j/length,"pivot move" 



j = 0
do i=1,length
 
    call rotate_sphere(change)
    if (change==1)then
    	j = j + 1
    end if 
end do
!print*, 1.0d0 * j/length,"rotate move"

!open(unit=50,file='azo_ini.txt')
!open(unit=51,file='sub_ini.txt')
!open(unit=52,file='pol_ini.txt')
!do j=1,n_azo
!    do i=0,nm
!        write(50,*) azo(j,i)%x, azo(j,i)%y, azo(j,i)%z
!    end do
!end do
!do j=1,n_sub
!    do i=0,nm_sub
!        write(51,*) sub(j,i)%x, sub(j,i)%y, sub(j,i)%z
!    end do
!end do
!do j=1,n_chain
!    do i=0,nm_pol
!        write(52,*) polymer(j,i)%x, polymer(j,i)%y, polymer(j,i)%z
!    end do
!end do
!close(50)
!close(51)
!close(52)

 
call checkpolymer (flag_c)

Deallocate(bond_vector)    

!call comformation_write()
!write(*,*) "CREATED OK"

stop"ini is ok"
end subroutine initialize          
