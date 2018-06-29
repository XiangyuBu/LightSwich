subroutine var_s()
USE global_parameters
USE utility_routines    

implicit none
INCLUDE 'mpif.h' 

Integer :: i_count
Integer :: df,dff,dfff,dffff
Integer :: i, j, k,j_temp, is, rotate_i
integer :: i_move,i_rotate, i_pivot, i_small, i_azo_test
integer :: n_move,n_rotate, n_pivot, n_small  
integer :: change ! change is the flag of MC, 1 for accepte the move.
 
DOUBLE PRECISION :: zz, rr, r_radius 
DOUBLE PRECISION :: r, End1_End2, rhorho,End1End2
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: pp,ppp,pppp,ppppp,psi,rho,subend 
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: density,density_azo, density_end_1, density_end_2,density_end_3
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: pp_temp,ppp_temp,pppp_temp,ppppp_temp,subend_temp
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: density_temp,density_azo_temp
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: density_end_1_temp, density_end_2_temp,density_end_3_temp

character*7 res,resres,res_s,res_o,res_e,res_p,res_a
character res0,res1,res2
allocate( pp(1:Nz),ppp(1:Nz),pppp(1:Nz),ppppp(1:Nz),psi(1:Nz),rho(1:Nz),subend(1:Nz) )
allocate( pp_temp(1:Nz),ppp_temp(1:Nz),pppp_temp(1:Nz),ppppp_temp(1:Nz),subend_temp(1:Nz) )
allocate( density_end_1(1:Nr+1,1:Nz),density_end_2(1:Nr+1,1:Nz),density_end_3(1:Nr+1,1:Nz) )
allocate( density_azo(1:Nr+1,1:Nz), density(1:Nr+1,1:Nz) )
allocate( density_end_1_temp(1:Nr+1,1:Nz),density_end_2_temp(1:Nr+1,1:Nz),density_end_3_temp(1:Nr+1,1:Nz) )
allocate( density_azo_temp(1:Nr+1,1:Nz), density_temp(1:Nr+1,1:Nz) )

i_count = 0
open(unit=61,file='azo.txt')
open(unit=62,file='polymer.txt')
open(unit=63,file='density.txt')
open(unit=44,file='energy.txt')
open(unit=81,file='p_df.txt')        !p_sphere_2 distribution fuction
open(unit=82,file='azo_end.txt')
open(unit=83,file='polymer_end.txt')
open(unit=85,file='azobond.txt')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!open(unit=100,file='w.ome')

    call SCMFT()

!do j = 1, Nz
!    do i = 1, Nr
!		read(100,*)  w(i,j)
!	end do
!end do

!close(100)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
MCS = 0
density = 0
density_end_1 = 0
density_end_2 = 0
density_end_3 = 0
pp = 0
ppp = 0
pppp = 0
ppppp = 0
subend = 0
psi = 0
rho = 0
!!!!!!!!!!!!!!!!!!!!!!!!!!pre-move!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        N_pre = 0
        n_move = 0
        n_rotate = 0 
        n_pivot = 0
        n_small = 0
        i_move = 0
        i_rotate = 0 
        i_pivot = 0
        i_small = 0
       
        do while(N_pre < Npre)              
            r = ran2(seed)
            if ( r>trial_move_rate(1) ) then
            
                r = ran2(seed)
                if ( r<trial_move_rate(2) ) then
                    n_pivot = n_pivot + 1
                    call pivot(change)
                    if (change == 1) then
                        i_pivot = i_pivot + 1
                    end if 
                else if (r>=trial_move_rate(2) .and. r<=trial_move_rate(3) ) then

                    n_rotate = n_rotate + 1
                    call rotate_sphere(change)
                    if (change == 1) then
                        i_rotate = i_rotate + 1
                    end if 

                else if ( r>trial_move_rate(3) ) then

                    n_move = n_move + 1
                    call move_sphere(change)
                    if (change == 1) then
                        i_move = i_move + 1
                    end if 

                end if
                 
            else 
                r = ran2(seed)
                call pivot_azo(change)
                if (change == 1) then
                    i_small = i_small + 1
                end if    
            end if 
            
            if(change == 1)then
                N_pre = 1 + N_pre         
            end if
        end do  
!!!!!!!!!!!!!!!!!!!!!!!!!!end pre-move!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

print*, "var is beginning"
do while(MCS < 5*NMCs)
               
    MCS = MCS + 1
  
    moves = 0
        do while(moves < Nmove)              
        
            r = ran2(seed)
            if ( r>trial_move_rate(1) ) then
            
                r = ran2(seed)
                if ( r<trial_move_rate(2) ) then
                    call pivot(change)
                else if (r>=trial_move_rate(2) .and. r<=trial_move_rate(3) ) then
                    call rotate_sphere(change)
                else if ( r>trial_move_rate(3) ) then
                    call move_sphere(change)
                end if
            else 
                r = ran2(seed)
                call pivot_azo(change)
            end if 

            if(change == 1)then
                moves = 1 + moves         
            end if  
    
        end do   

    dffff = 0

    do j=1,N_azo 
        if(azo_position_1>=1.0d0) then
            i_azo_test = floor((azo_position_1 - 1.0d0) * Nm_azo)    !it is related to azo_position
            dffff = floor( ( Lz_2 + azo(j,i_azo_test)%z ) / dz ) + 1      !!!iz_azo(j,Nm)
            if (dffff < 1 .or.  dffff > 500)   then
                stop "error dffff"
            end if
            ppppp(dffff) = ppppp(dffff) + 1   
        else
            dffff = floor( ( Lz_2 + azo(j,i_azo(j))%z ) / dz ) + 1      !!!iz_azo(j,Nm)
            if (dffff < 1 .or.  dffff > 500)   then
                stop "error dffff"
            end if
            ppppp(dffff) = ppppp(dffff) + 1
        end if
    end do

    dfff = 0
    do j=1,N_chain
        dfff = floor( ( Lz_2 + polymer(j,Nm_pol)%z ) / dz ) + 1 !!!iz(j,Nm)
        if (dfff < 1 .or.  dfff > 500)   then
            stop "error dfff"
        end if
        pppp(dfff) = pppp(dfff) + 1
    end do 

    dff = 0
    do j=1,N_azo
            dff = floor( ( Lz_2 + azo(j,Nm_azo)%z ) / dz ) + 1      !!!iz_azo(j,Nm)
            if (dff < 1 .or.  dff > 500)   then
                stop "error dff"
            end if
            ppp(dff) = ppp(dff) + 1    
    end do 

    df = 0 
    df = floor( (Lz_2 + np_r) / dz ) + 1
    if (df < 1 .or.  df > 500)   then
        stop "error df"
    end if
    pp(df) = pp(df) + 1  

    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!    cal density
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    
!    do j=1,N_chain
!        do i=1,Nm_pol
!            density(ir(j,i),iz(j,i)) = density(ir(j,i),iz(j,i)) + 1                
!        end do
!        density_end_1(ir(j,Nm_pol),iz(j,Nm_pol)) = density_end_1(ir(j,Nm_pol),iz(j,Nm_pol)) + 1
!    end do

    do j=1,N_azo
        do i=1,Nm_azo
            density_azo(ir_azo(j,i),iz_azo(j,i)) = density_azo(ir_azo(j,i),iz_azo(j,i)) + 4                
        end do
            density_end_2(ir_azo(j,Nm_azo),iz_azo(j,Nm_azo)) = density_end_2(ir_azo(j,Nm_azo),iz_azo(j,Nm_azo)) + 4
    end do
  
end do   ! MCS
close(64)
density = deltaS*density/MCS
density_azo = deltaS*density_azo/MCS
!density_end_1 = density_end_1/MCS
density_end_2 = density_end_2/MCS
!density_end_3 = density_end_3/MCS
pp = pp/MCS
ppp = ppp/MCS
pppp = pppp/MCS
ppppp = ppppp/MCS
print*, "var_mcs is ok"
print*, myid, "Ok"
CALL MPI_barrier(MPI_COMM_WORLD, ierr)
print*, myid, "barrier is ok"
call mpi_allreduce(density(1,1), density_temp(1,1), (Nr+1)*Nz, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
call mpi_allreduce(density_azo(1,1), density_azo_temp(1,1), (Nr+1)*Nz, mpi_double_precision, mpi_sum, mpi_comm_world, ierr) 
!call mpi_allreduce(density_end_1(1,1), density_end_1_temp(1,1), (Nr+1)*Nz, mpi_double_precision, mpi_sum, mpi_comm_world, ierr) 

call mpi_allreduce(density_end_2(1,1), density_end_2_temp(1,1), (Nr+1)*Nz, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
!call mpi_allreduce(density_end_3(1,1), density_end_3_temp(1,1), (Nr+1)*Nz, mpi_double_precision, mpi_sum, mpi_comm_world, ierr) 

call mpi_allreduce(pp, pp_temp, Nz, mpi_double_precision, mpi_sum, mpi_comm_world, ierr) 
call mpi_allreduce(ppp, ppp_temp, Nz, mpi_double_precision, mpi_sum, mpi_comm_world, ierr) 
call mpi_allreduce(pppp, pppp_temp, Nz, mpi_double_precision, mpi_sum, mpi_comm_world, ierr) 
call mpi_allreduce(ppppp, ppppp_temp, Nz, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
density=density_temp/numprocs
density_azo=density_azo_temp/numprocs
!density_end_1=density_end_1_temp/numprocs
density_end_2=density_end_2_temp/numprocs
!density_end_3=density_end_3_temp/numprocs
pp=pp_temp/numprocs
ppp=ppp_temp/numprocs
pppp=pppp_temp/numprocs
ppppp=ppppp_temp/numprocs

print*,"mpi is ok"

if (myid==0) then
    zz = Lz/Nz/Nm
    rr = Lr/Nr/Nm

do j=1,nz
    do i=1,nr  
        density (i,j) = density (i,j)/( r_a(i)*r_dr*r_dz )
        density_azo (i,j) = density_azo (i,j)/( r_a(i)*r_dr*r_dz )
!        density_end_1 (i,j) = density_end_1 (i,j)/( r_a(i)*r_dr*r_dz )
        density_end_2 (i,j) = density_end_2 (i,j)/( r_a(i)*r_dr*r_dz )
!        density_end_3 (i,j) = density_end_3 (i,j)/( r_a(i)*r_dr*r_dz )
    end do
end do

do df=1,500
    write(81,*)  zz*df, pp(df)
end do

do dff=1,Nz
    write(82,*)  zz*dff, ppp(dff)
end do

do dfff=1,Nz
    write(83,*)  zz*dfff, pppp(dfff)
end do

do dffff=1,Nz
    write(85,*)  zz*dffff, ppppp(dffff)
end do

!do j=1,Nz
!End1End2 = 0
!    do i=1,Nr  
!        End1End2 = End1End2 + density_end_1 (i,j) * density_end_2 (i,j) * rr_r(i) 
!    end do
!    psi(j) = 4*pi * End1End2 * r_dr * r_dz
!    write(84,*)  zz*j, psi(j)
!end do

!do j=1,Nz
!rhorho = 0
!    do i=1,Nr  
!        rhorho = rhorho + density (i,j) * density (i,j) * rr_r(i) 
!    end do
!    rho(j) = 4*pi* rhorho * r_dr * r_dz
!    write(86,*)  zz*j, rho(j)
!end do

close(81)
close(82)
close(83)
close(85)
!close(84)
!close(86)

!End1_End2 = 0
!do j=1,nz
!    do i=1,nr  
!        End1_End2 = End1_End2 + density_end_1 (i,j) * density_end_2 (i,j) * rr_r(i) 
!    end do
!end do
!End1_End2 = 4*pi* End1_End2 * r_dr * r_dz

!rhorho = 0
!do j=1,nz
!    do i=1,nr  
!        rhorho = rhorho + density (i,j) * density (i,j) * rr_r(i) 
!    end do
!end do
!rhorho = 4*pi* rhorho * r_dr * r_dz

!write(44,"(6E25.13)") hahah, csoL, rhorho, End1_End2

do j=1,N_azo
    do i=0,Nm_azo
        write(61,*)  azo(j,i)%x, azo(j,i)%y, azo(j,i)%z
    end do
end do
close(61)

do j=1,N_chain
    do i=0,Nm_pol
        write(62,*)  polymer(j,i)%x, polymer(j,i)%y, polymer(j,i)%z
    end do
end do
close(62)



do j=1,nz
    do i=1,nr
        write(63,"(7E25.13)") rr*i,zz*j,density_azo (i,j),density_end_2(i,j)
    end do
end do
close(63)
!close(44)

end if
end subroutine var_s
