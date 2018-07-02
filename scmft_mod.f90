subroutine SCMFT()
USE global_parameters
USE utility_routines    

implicit none

Integer :: i, j, k,j_temp, is,ii, rotate_i 
integer :: change ! change is the flag of MC, 1 for accepte the move.
integer :: i_move,i_rotate, i_pivot, i_small
integer :: n_move,n_rotate, n_pivot, n_small 
integer :: flag_c
DOUBLE PRECISION, PARAMETER :: TOL = 1.5D-3
DOUBLE PRECISION :: E2_test
DOUBLE PRECISION ::  w_erro, w_erromax,trial, r, End1_End2, rhorho 
DOUBLE PRECISION :: derro(1:70),erro(1:70),errodown
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: density, density_end_1, density_end_2, density_end_3

character*7 res,resres,res_s,res_o,res_p,res_a,res_sub,resresres
character res0,res1,res2

allocate( density_end_3(1:Nr+1,1:Nz),density_end_1(1:Nr+1,1:Nz),density_end_2(1:Nr+1,1:Nz), density(1:Nr+1,1:Nz) )

!!! iteration
 
do n_iter = 1, Max_iter
    print*, "start", n_iter,"iteration"

    ! pre-move
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
!        write(15,*) i_move, 1.0d0 * i_move/n_move,"sphere move"
!        write(15,*) i_rotate, 1.0d0 * i_rotate/n_rotate,"rotate move"
!        write(15,*) i_pivot, 1.0d0 * i_pivot/n_pivot,"pivot move"    
!        write(15,*) i_small, 1.0d0 * i_small/n_small,"azo move"
!        print*,"pre is ok"
!! find out w_new
    
    MCS = 0
    density = 0
    density_end_1 = 0
    density_end_2 = 0
    density_end_3 = 0

    do while(MCS < NMCs)      
        
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
!        do j=1,N_chain
!            do i=1,Nm_pol
!                density(ir(j,i),iz(j,i)) = density(ir(j,i),iz(j,i)) + 1                
!            end do
!            density_end_1(ir(j,Nm_pol),iz(j,Nm_pol)) = density_end_1(ir(j,Nm_pol),iz(j,Nm_pol)) + 1
!        end do

        do j=1,N_azo
            do i=1,Nm_azo
                density(ir_azo(j,i),iz_azo(j,i)) = density(ir_azo(j,i),iz_azo(j,i)) + 4                
            end do
                density_end_2(ir_azo(j,Nm_azo),iz_azo(j,Nm_azo)) = density_end_2(ir_azo(j,Nm_azo),iz_azo(j,Nm_azo)) + 4
        end do

    end do   ! MCS

    density = deltaS*density / MCS 
!    density_end_1 = density_end_1 / MCS 
    density_end_2 = density_end_2 / MCS 

    do j=1,nz
        do i=1,nr  
            density (i,j) = density (i,j)/( r_a(i)*r_dr*r_dz )
!            density_end_1 (i,j) = density_end_1 (i,j)/( r_a(i)*r_dr*r_dz )
            density_end_2 (i,j) = density_end_2 (i,j)/( r_a(i)*r_dr*r_dz )
        end do
    end do

    w_new = nu * density
!    eta_new = tau * density_end_2
!    eta_azo_new = tau * density_end_1
   ! compute erros
    w_erro = 0

    do j = 1, Nz
        do i = 1, Nr
            w_erro = w_erro + abs(w_new(i,j) - w(i,j)) 
!                            + abs(eta_new(i,j) - eta(i,j)) &
!                            + abs(eta_azo_new(i,j) - eta_azo(i,j)) 
       end do
    end do
 
    w_erro = 0.3d0*w_erro/Nz/Nr    
    
    print*, "SCMFT", n_iter, w_erro



    if ( n_iter == 1 ) then
        open(unit=33,file='density_raw.dat')    
            do j = 1, Nz
                do i = 1, 4*Nr/5
                    write(33,"(7E25.13)") rr_r(i), zz_r(j), density(i,j),density_end_2 (i,j)
                end do
            end do                   
        close(33)
    end if


!stop "first scft"

    if ( mod(n_iter,5) == 0 ) then
        open(unit=61,file='w_temp.ome')
        do j = 1, Nz
            do i = 1, Nr
                write(61,*) n_iter, w(i,j)
            end do
        end do
        close(61)
    end if


    errodown = -1.0d0
    erro(n_iter) = w_erro
    if (w_erro<TOL .and. n_iter>3) then
        open(unit=61,file='w.ome')
        do j = 1, Nz
            do i = 1, Nr
                write(61,*) w(i,j)
            end do
        end do
        close(61)
        exit
    end if
    
    if (n_iter>5) then
        
        derro(n_iter) = erro(n_iter) - erro(n_iter-1)
        if (n_iter>10) then
            errodown = derro(n_iter) + derro(n_iter-1) + derro(n_iter-2) + derro(n_iter-3) + derro(n_iter-4) 
        else 
            errodown = -1.0d0                   
        end if
    end if
    if (errodown>0.0d0) then
        print*,"errodown=",errodown
        open(unit=61,file='w.ome')
        do j = 1, Nz
            do i = 1, Nr
		            write(61,*) w(i,j)
	          end do
        end do
        close(61)
        exit
    end if    
    
    !simple mixing scheme
    w = lambda*w_new + (1-lambda)*w
!    eta = 0.5d0*lambda*eta_new + 0.5d0*(1-lambda)*eta
!    eta_azo = 0.5d0*lambda*eta_azo_new + 0.5d0*(1-lambda)*eta_azo
    ! boundary condition
    do j=1,nz
       w(nr+1,j) = w(nr,j)
!       eta(nr+1,j) = eta(nr,j)
!       eta_azo(nr+1,j) = eta_azo(nr,j)
    end do
        
 
    call checkpolymer (flag_c)  
           
end do  ! enddo n_iter

stop "meanfield is over"
end subroutine SCMFT
