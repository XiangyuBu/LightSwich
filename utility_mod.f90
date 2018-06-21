MODULE utility_routines
IMPLICIT NONE
CONTAINS
subroutine move_sphere (change)
use global_parameters
implicit none
integer i, j, k
integer :: change
Integer :: iz_temp(1:N_chain,0:Nm_pol)
DOUBLE PRECISION :: r_radius, r
DOUBLE PRECISION ::  DE2, DE3 
DOUBLE PRECISION :: z_temp
DOUBLE PRECISION :: z_move, p_sphere_try

         
DE2 = 0
DE3 = 0
z_move = 0
p_sphere_try = 0
do while(z_move==0 .or. p_sphere_try<=r_sphere .or. (p_sphere_try+nm*roL+Nm_pol)>=Lz)
    z_move = 2.0d0*move_max*(2.0d0*ran2(seed) - 1)      
    p_sphere_try = p_sphere_2 + z_move
end do
change = 1

if (z_move>=0) then       ! if move away
    change = 1
    do j = 1, n_chain
        do i= 0, Nm_pol 
    	      if (polymer(j,i)%z > p_sphere_try ) then    
                change = 0
                exit
            end if
        end do
        if (change == 0)then
            exit
        end if
    end do
    if(change==1) then   
    do j = 1, n_chain
        do i= 0, Nm_pol
            iz_temp(j,i) = floor( ( p_sphere_try - polymer(j,i)%z ) / dz ) + 1
            if (iz_temp(j,i) < 1 .or.  iz_temp(j,i) > 500)   then
                stop "111"
            end if
            DE2 = DE2 + w(ir(j,i), iz_temp(j,i)) - w(ir(j,i), iz(j,i))            
        end do
        DE3 = DE3 + eta(ir(j,Nm_pol), iz_temp(j,Nm_pol)) - eta(ir(j,Nm_pol), iz(j,Nm_pol))
    end do
    r = ran2(seed)
    if ( r < dexp ( - deltaS*DE2 + DE3 ) )then
        p_sphere_2 = p_sphere_try        
        do j = 1, n_sub
            do i = 0,Nm_sub          
                sub(j,i)%z = sub(j,i)%z + z_move
            end do
        end do
        do j = 1, n_azo
            do i = 0,Nm          
                azo(j,i)%z = azo(j,i)%z + z_move
            end do
        end do
        do j = 1, n_chain
            do i = 0,Nm_pol          
                iz(j,i) = iz_temp(j,i)
            end do
        end do 
    else
        change = 0
    end if
    end if
else    ! if z_move<0
    change = 1
    do j = 1, n_chain
        do i= 0, Nm_pol 
    	      if (polymer(j,i)%z > p_sphere_try ) then    
                change = 0
                exit
            end if
        end do
        if (change == 0)then
            exit
        end if
    end do   
    if (change == 1) then
        do j = 1, n_azo
            do i= 0, Nm
    	       z_temp = azo(j,i)%z + z_move
    	       r_radius = azo(j,i)%x*azo(j,i)%x + azo(j,i)%y*azo(j,i)%y + z_temp*z_temp
    	       if (r_sphere_2 > r_radius ) then    
                    change = 0
                    exit
                end if
            end do
            if (change == 0)then
                exit
            end if
        end do       
    end if
    if (change == 1) then
        do j = 1, n_sub
            do i= 0, Nm_sub
    	       z_temp = sub(j,i)%z + z_move
    	       r_radius = sub(j,i)%x*sub(j,i)%x + sub(j,i)%y*sub(j,i)%y + z_temp*z_temp
    	       if (r_sphere_2 > r_radius ) then    
                    change = 0
                    exit
                end if
            end do
            if (change == 0)then
                exit
            end if
        end do       
    end if    
!        print*, change,"azo"
    if (change == 1) then    
	     do j = 1, n_chain
    	     do i= 0, Nm_pol
        	     iz_temp(j,i) = floor( ( p_sphere_try - polymer(j,i)%z ) / dz ) + 1
               if (iz_temp(j,i) < 1 .or.  iz_temp(j,i) > 500)   then
                  stop "222"
               end if                
        	     DE2 = DE2 + w(ir(j,i), iz_temp(j,i)) - w(ir(j,i), iz(j,i))           
    	     end do
    	     DE3 = DE3 + eta(ir(j,Nm_pol), iz_temp(j,Nm_pol)) - eta(ir(j,Nm_pol), iz(j,Nm_pol))
        end do
        r = ran2(seed)
	    if ( r < dexp ( - deltaS*DE2 + DE3 ) )then
    	    p_sphere_2 = p_sphere_try        
    	    do j = 1, n_sub
        	    do i = 0,Nm_sub          
                    sub(j,i)%z = sub(j,i)%z + z_move
        	    end do
    	    end do
    	    do j = 1, n_azo
        	    do i = 0,Nm          
                    azo(j,i)%z = azo(j,i)%z + z_move
        	    end do
    	    end do
    	    do j = 1, n_chain
        	    do i = 0,Nm_pol          
            	    iz(j,i) = iz_temp(j,i)
        	    end do
    	    end do 
	    else
    	    change = 0

        end if 
    end if
end if
   
end subroutine move_sphere 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine pivot_azo (change)
    use global_parameters
    implicit none
    integer i, j, length
    integer :: change
    type(node) :: new(0:401)
    Integer :: ir_temp(0:401), iz_temp(0:401)

    DOUBLE PRECISION :: axis(3)
    DOUBLE PRECISION :: r_radius, r, cos_t, sin_t, phi
    DOUBLE PRECISION :: unew(3), uold(3)  
    DOUBLE PRECISION :: alpha, beta, angle, dotp
    DOUBLE PRECISION :: x_r, y_r
    DOUBLE PRECISION :: DE1, DE2, DE3
    change = 0
    change = 1
    !this is the big MC move
    jj = floor( ran2(seed)*0.9999999d0*(N_azo) ) + 1 ! random pickup the chain in [1,N_azo] to be rotated 
    i = floor(ran2(seed)*0.9999999d0*Nm)  ! random pickup the monomer in [0,Nm-1] to be rotated         
    if ( i/= i_azo(jj) )then
        cos_t=(2*ran2(seed)-1)*0.99999999d0
        sin_t=dsqrt(1.0d0-cos_t**2) 
        phi = (2*ran2(seed) - 1)*pi      !generate angle of phi \in (-pi,+pi)
        axis(1) = sin_t*dcos(phi)
        axis(2) = sin_t*dsin(phi)
        axis(3) = cos_t
!       print*,"azo",jj
    else
        axis(1) = azo(jj,i)%x - azo(jj,i-1)%x
        axis(2) = azo(jj,i)%y - azo(jj,i-1)%y
        axis(3) = azo(jj,i)%z - azo(jj,i-1)%z               
    end if

    angle = rotate*(2*ran2(seed) - 1)*pi      !generate angle of phi \in (-pi,+pi)
    alpha = dcos(angle)
    beta = dsin(angle)

    do j = i+1,Nm     
        uold(1) = azo(jj,j)%x - azo(jj,i)%x
        uold(2) = azo(jj,j)%y - azo(jj,i)%y
        uold(3) = azo(jj,j)%z - azo(jj,i)%z
                        
        dotp = axis(1)*uold(1) + axis(2)*uold(2) + axis(3)*uold(3)
        unew(1) = uold(1)*alpha + axis(1)*dotp*(1-alpha) + ( uold(2)*axis(3) - uold(3)*axis(2) )*beta
        unew(2) = uold(2)*alpha + axis(2)*dotp*(1-alpha) + ( uold(3)*axis(1) - uold(1)*axis(3) )*beta
        unew(3) = uold(3)*alpha + axis(3)*dotp*(1-alpha) + ( uold(1)*axis(2) - uold(2)*axis(1) )*beta

        new(j)%x = azo(jj,i)%x + unew(1)
        new(j)%y = azo(jj,i)%y + unew(2)
        new(j)%z = azo(jj,i)%z + unew(3)
        r_radius = new(j)%x*new(j)%x + new(j)%y*new(j)%y + new(j)%z*new(j)%z

        if (r_sphere_2 > r_radius .or.  new(j)%z > p_sphere_2  ) then
            change = 0 
            exit
        end if
    end do 
        
if (change == 1) then

    if (i == 0) then
        DE1 = 0
    else
        DE1 = (new(i+1)%x - 2*azo(jj,i)%x + azo(jj,i-1)%x)**2   &
            + (new(i+1)%y - 2*azo(jj,i)%y + azo(jj,i-1)%y)**2   &  
            + (new(i+1)%z - 2*azo(jj,i)%z + azo(jj,i-1)%z)**2   &
            - (azo(jj,i+1)%x - 2*azo(jj,i)%x + azo(jj,i-1)%x)**2   &
            - (azo(jj,i+1)%y - 2*azo(jj,i)%y + azo(jj,i-1)%y)**2   &  
            - (azo(jj,i+1)%z - 2*azo(jj,i)%z + azo(jj,i-1)%z)**2
           
    end if   !endif i

    DE2 = 0.0d0    
    
    do j = i+1, Nm
        if ( new(j)%x<=Lbox .and. new(j)%y<=Lbox ) then
            x_r = new(j)%x
            y_r = new(j)%y             
        else if ( new(j)%x>Lbox .and. new(j)%y>Lbox ) then
            x_r = new(j)%x - 2*Lbox
            y_r = new(j)%y - 2*Lbox
        else if ( new(j)%x>Lbox ) then
            x_r = new(j)%x - 2*Lbox
            y_r = new(j)%y
        else
            x_r = new(j)%x
            y_r = new(j)%y - 2*Lbox
        end if   
        iz_temp(j) = floor( ( p_sphere_2 - new(j)%z ) / dz ) + 1
        r_radius = dsqrt( x_r*x_r + y_r*y_r )
        if(r_radius<Lr)then
            ir_temp(j) = floor( r_radius / dr ) + 1
    	  else 
            ir_temp(j) = Nr + 1
		    end if        
        DE2 = DE2 + w(ir_temp(j), iz_temp(j)) - w(ir_azo(jj,j), iz_azo(jj,j))
        
    end do
    

    DE3 = eta_azo(ir_temp(Nm),iz_temp(Nm)) - eta_azo(ir_azo(jj,Nm),iz_azo(jj,Nm))

   
    r = ran2(seed)
   
    if ( r < dexp ( -epsilon_azo*DE1 - deltaS*DE2 + DE3 ))then

        do j = i+1,Nm          
            azo(jj,j)%x = new(j)%x
            azo(jj,j)%y = new(j)%y
            azo(jj,j)%z = new(j)%z
            iz_azo(jj,j) = iz_temp(j)
            ir_azo(jj,j) = ir_temp(j)
        end do  

    else
        change = 0
    end if
 
end if       

end subroutine pivot_azo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine pivot_sub (change)
    use global_parameters
    implicit none
    integer i, j, length
    integer :: change
    type(node) :: new(0:401)
    Integer :: ir_temp(0:401), iz_temp(0:401)

    DOUBLE PRECISION :: axis(3)
    DOUBLE PRECISION :: r_radius, r, cos_t, sin_t, phi
    DOUBLE PRECISION :: unew(3), uold(3)  
    DOUBLE PRECISION :: alpha, beta, angle, dotp
    DOUBLE PRECISION :: x_r, y_r
    DOUBLE PRECISION :: DE1, DE2, DE3
    change = 0
    change = 1
    !this is the big MC move
    jj = floor(ran2(seed)*0.9999999d0*(N_sub)) + 1 ! random pickup the chain in [1,N_sub] to be rotated 
    i = floor(ran2(seed)*0.9999999d0*Nm_sub)  ! random pickup the monomer in [0,Nm-1] to be rotated         


    cos_t=(2*ran2(seed)-1)*0.99999999d0
    sin_t=dsqrt(1.0d0-cos_t**2) 
    phi = (2*ran2(seed) - 1)*pi      !generate angle of phi \in (-pi,+pi)
    axis(1) = sin_t*dcos(phi)
    axis(2) = sin_t*dsin(phi)
    axis(3) = cos_t

    angle = rotate*(2*ran2(seed) - 1)*pi      !generate angle of phi \in (-pi,+pi)
    alpha = dcos(angle)
    beta = dsin(angle)

    do j = i+1,Nm_sub     
        uold(1) = sub(jj,j)%x - sub(jj,i)%x
        uold(2) = sub(jj,j)%y - sub(jj,i)%y
        uold(3) = sub(jj,j)%z - sub(jj,i)%z
                        
        dotp = axis(1)*uold(1) + axis(2)*uold(2) + axis(3)*uold(3)
        unew(1) = uold(1)*alpha + axis(1)*dotp*(1-alpha) + ( uold(2)*axis(3) - uold(3)*axis(2) )*beta
        unew(2) = uold(2)*alpha + axis(2)*dotp*(1-alpha) + ( uold(3)*axis(1) - uold(1)*axis(3) )*beta
        unew(3) = uold(3)*alpha + axis(3)*dotp*(1-alpha) + ( uold(1)*axis(2) - uold(2)*axis(1) )*beta

        new(j)%x = sub(jj,i)%x + unew(1)
        new(j)%y = sub(jj,i)%y + unew(2)
        new(j)%z = sub(jj,i)%z + unew(3)
        r_radius = new(j)%x*new(j)%x + new(j)%y*new(j)%y + new(j)%z*new(j)%z

        if (r_sphere_2 > r_radius .or.  new(j)%z > p_sphere_2  ) then
            change = 0 
            exit
        end if
    end do 
        
if (change == 1) then

    if (i == 0) then
        DE1 = 0
    else
        DE1 = (new(i+1)%x - 2*sub(jj,i)%x + sub(jj,i-1)%x)**2   &
            + (new(i+1)%y - 2*sub(jj,i)%y + sub(jj,i-1)%y)**2   &  
            + (new(i+1)%z - 2*sub(jj,i)%z + sub(jj,i-1)%z)**2   &
            - (sub(jj,i+1)%x - 2*sub(jj,i)%x + sub(jj,i-1)%x)**2   &
            - (sub(jj,i+1)%y - 2*sub(jj,i)%y + sub(jj,i-1)%y)**2   &  
            - (sub(jj,i+1)%z - 2*sub(jj,i)%z + sub(jj,i-1)%z)**2
           
    end if   !endif i

    DE2 = 0.0d0    
    
    do j = i+1, Nm_sub
        if ( new(j)%x<=Lbox .and. new(j)%y<=Lbox ) then
            x_r = new(j)%x
            y_r = new(j)%y             
        else if ( new(j)%x>Lbox .and. new(j)%y>Lbox ) then
            x_r = new(j)%x - 2*Lbox
            y_r = new(j)%y - 2*Lbox
        else if ( new(j)%x>Lbox ) then
            x_r = new(j)%x - 2*Lbox
            y_r = new(j)%y
        else
            x_r = new(j)%x
            y_r = new(j)%y - 2*Lbox
        end if   
        iz_temp(j) = floor( ( p_sphere_2 - new(j)%z ) / dz ) + 1
        r_radius = dsqrt( x_r*x_r + y_r*y_r )
        if(r_radius<Lr)then
            ir_temp(j) = floor( r_radius / dr ) + 1
    	  else 
            ir_temp(j) = Nr + 1
		    end if        
        DE2 = DE2 + w(ir_temp(j), iz_temp(j)) - w(ir_sub(jj,j), iz_sub(jj,j))
        
    end do
       
    r = ran2(seed)
   
    if ( r < dexp ( -epsilon_azo*DE1 - deltaS*DE2 ))then

        do j = i+1,Nm_sub          
            sub(jj,j)%x = new(j)%x
            sub(jj,j)%y = new(j)%y
            sub(jj,j)%z = new(j)%z
            iz_sub(jj,j) = iz_temp(j)
            ir_sub(jj,j) = ir_temp(j)
        end do  

    else
        change = 0
    end if
 
end if       

end subroutine pivot_sub
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine pivot (change)
    use global_parameters
    implicit none
    integer i, j, length
    integer :: change
	
    type(node) :: new(0:401)
    Integer :: ir_temp(0:401), iz_temp(0:401)

    DOUBLE PRECISION :: axis(3)
    DOUBLE PRECISION :: r_radius,r_radius_1, r, cos_t, sin_t, phi
    DOUBLE PRECISION :: unew(3), uold(3)  
    DOUBLE PRECISION :: alpha, beta, angle, dotp    
    DOUBLE PRECISION :: DE1, DE2, DE3  
    change = 0
        change = 1
        !this is the big MC move
        jj =  1 + floor( ran2(seed)*0.9999999d0*(N_chain) ) ! random pickup the chain in [1,N_chain] to be rotated        
        i = i_azo(jj)
        do while (i == i_azo(jj))
            i = floor(ran2(seed)*0.9999999d0*Nm_pol)  ! random pickup the monomer in [0,Nm-1] to be rotated 
        end do
        length = Nm_pol - i 

        cos_t=(2*ran2(seed)-1)*0.99999999d0
        sin_t=dsqrt(1.0d0-cos_t**2) 

        phi = (2*ran2(seed) - 1)*pi      !generate angle of phi \in (-pi,+pi)
        axis(1) = sin_t*dcos(phi)
        axis(2) = sin_t*dsin(phi)
        axis(3) = cos_t

        angle = rotate*(2*ran2(seed) - 1)*pi      !generate angle of phi \in (-pi,+pi)
        alpha = dcos(angle)
        beta = dsin(angle)
 
        do j = i+1,Nm_pol     
            uold(1) = polymer(jj,j)%x - polymer(jj,i)%x
            uold(2) = polymer(jj,j)%y - polymer(jj,i)%y
            uold(3) = polymer(jj,j)%z - polymer(jj,i)%z
                        
            dotp = axis(1)*uold(1) + axis(2)*uold(2) + axis(3)*uold(3)
            unew(1) = uold(1)*alpha + axis(1)*dotp*(1-alpha) + ( uold(2)*axis(3) - uold(3)*axis(2) )*beta
            unew(2) = uold(2)*alpha + axis(2)*dotp*(1-alpha) + ( uold(3)*axis(1) - uold(1)*axis(3) )*beta
            unew(3) = uold(3)*alpha + axis(3)*dotp*(1-alpha) + ( uold(1)*axis(2) - uold(2)*axis(1) )*beta

            new(j-i)%x = polymer(jj,i)%x + unew(1)
            new(j-i)%y = polymer(jj,i)%y + unew(2)
            new(j-i)%z = polymer(jj,i)%z + unew(3)
            r_radius = new(j-i)%x*new(j-i)%x + new(j-i)%y*new(j-i)%y + new(j-i)%z*new(j-i)%z

            if (r_sphere_2 > r_radius .or.  new(j-i)%z > p_sphere_2  ) then
                change = 0 
                exit
            end if
        end do 
      
if (change == 1) then

    if (i == 0) then
        DE1 = 0
    else
        DE1 = (new(1)%x - 2*polymer(jj,i)%x + polymer(jj,i-1)%x)**2   &
           + (new(1)%y - 2*polymer(jj,i)%y + polymer(jj,i-1)%y)**2   &  
           + (new(1)%z - 2*polymer(jj,i)%z + polymer(jj,i-1)%z)**2   &
           - (polymer(jj,i+1)%x - 2*polymer(jj,i)%x + polymer(jj,i-1)%x)**2   &
           - (polymer(jj,i+1)%y - 2*polymer(jj,i)%y + polymer(jj,i-1)%y)**2   &  
           - (polymer(jj,i+1)%z - 2*polymer(jj,i)%z + polymer(jj,i-1)%z)**2
           
    end if   !endif i

    r_radius = dsqrt( new(1)%x*new(1)%x + new(1)%y*new(1)%y )

    ir_temp(1) = floor( r_radius / dr ) + 1    
    iz_temp(1) = floor( ( p_sphere_2 - new(1)%z ) / dz ) + 1

    DE2 = w(ir_temp(1),iz_temp(1)) - w(ir(jj,i+1),iz(jj,i+1))    
    
    do j = 2, length
 
        r_radius = dsqrt( new(j)%x*new(j)%x + new(j)%y*new(j)%y )
        ir_temp(j) = floor( r_radius / dr ) + 1
        iz_temp(j) = floor( ( p_sphere_2 - new(j)%z ) / dz ) + 1 

        DE2 = DE2 + w(ir_temp(j), iz_temp(j)) - w(ir(jj,j+i), iz(jj,j+i))

    end do

    DE3 = eta(ir_temp(length),iz_temp(length)) - eta(ir(jj,Nm_pol),iz(jj,Nm_pol))

    r = ran2(seed)
   
    if ( r < dexp ( -epsilon*DE1 - deltaS*DE2 + DE3 ))then

        do j = i+1,Nm_pol          
            polymer(jj,j)%x = new(j-i)%x
            polymer(jj,j)%y = new(j-i)%y
            polymer(jj,j)%z = new(j-i)%z
            iz(jj,j) = iz_temp(j-i)
            ir(jj,j) = ir_temp(j-i)
        end do  

    else
        change = 0
    end if
 
end if       

end subroutine pivot 

subroutine rotate_sphere (change)
    use global_parameters
    implicit none
    integer i, j, k, l, length, trial_number
    integer :: change, change_1
    type(node) :: new(1:N_chain,0:Nm_pol)
    type(node) :: p_new(0:401)
    Integer :: ir_temp(1:N_chain,0:Nm_pol), iz_temp(1:N_chain,0:Nm_pol)

    DOUBLE PRECISION :: axis(3)
    DOUBLE PRECISION :: r_radius,r_radius_1, r, cos_t, sin_t, phi 
    DOUBLE PRECISION :: unew(3),uold(3)  
    DOUBLE PRECISION :: alpha, beta, angle, dotp
    DOUBLE PRECISION :: DE1, DE2, DE3 
    
    cos_t=(2*ran2(seed)-1)*0.99999999d0
    sin_t=dsqrt(1.0d0-cos_t**2) 

    phi = (2*ran2(seed) - 1)*pi      !generate angle of phi \in (-pi,+pi)
    axis(1) = sin_t*dcos(phi)
    axis(2) = sin_t*dsin(phi)
    axis(3) = cos_t

    angle = rotate_s*(2*ran2(seed) - 1)*pi      !generate angle of phi \in (-pi,+pi)
    alpha = dcos(angle)
    beta = dsin(angle)
 
    do j = 1,n_chain
        do i = 0,Nm_pol     
            uold(1) = polymer(j,i)%x 
            uold(2) = polymer(j,i)%y 
            uold(3) = polymer(j,i)%z 
                        
            dotp = axis(1)*uold(1) + axis(2)*uold(2) + axis(3)*uold(3)
            new(j,i)%x = uold(1)*alpha + axis(1)*dotp*(1-alpha) + ( uold(2)*axis(3) - uold(3)*axis(2) )*beta
            new(j,i)%y = uold(2)*alpha + axis(2)*dotp*(1-alpha) + ( uold(3)*axis(1) - uold(1)*axis(3) )*beta
            new(j,i)%z = uold(3)*alpha + axis(3)*dotp*(1-alpha) + ( uold(1)*axis(2) - uold(2)*axis(1) )*beta                
        end do
    end do
    if ( (rol*Nm+Nm_pol) < p_sphere_2 ) then
        change_1 = 1
        DE1 = 0
    else
      	DE1 = 0
        do j = 1,n_chain
            change_1 = 1
            do i=0,Nm_pol 
                if (new(j,i)%z > p_sphere_2) then                   
                    change_1 = 0
                    exit
                end if
            end do
            trial_number = 0 
            do while(change_1 == 0)                 ! make sure wall is impentrate
                change_1 = 1 
                trial_number = trial_number + 1
                i = floor( ran2(seed)*0.9999999d0*Nm_pol )  ! random pickup the monomer in [0,Nm-1] to be rotated 
                    
                cos_t=(2*ran2(seed)-1)*0.99999999d0
                sin_t=dsqrt(1.0d0-cos_t**2) 

                phi = (2*ran2(seed) - 1)*pi      !generate angle of phi \in (-pi,+pi)
                axis(1) = sin_t*dcos(phi)
                axis(2) = sin_t*dsin(phi)
                axis(3) = cos_t

                angle = rotate*(2*ran2(seed) - 1)*pi      !generate angle of phi \in (-pi,+pi)
                alpha = dcos(angle)
                beta = dsin(angle)
 
                do k = i+1,Nm_pol     
                    uold(1) = new(j,k)%x - new(j,i)%x
                    uold(2) = new(j,k)%y - new(j,i)%y
                    uold(3) = new(j,k)%z - new(j,i)%z
                        
                    dotp = axis(1)*uold(1) + axis(2)*uold(2) + axis(3)*uold(3)
                    unew(1) = uold(1)*alpha + axis(1)*dotp*(1-alpha) + ( uold(2)*axis(3) - uold(3)*axis(2) )*beta
                    unew(2) = uold(2)*alpha + axis(2)*dotp*(1-alpha) + ( uold(3)*axis(1) - uold(1)*axis(3) )*beta
                    unew(3) = uold(3)*alpha + axis(3)*dotp*(1-alpha) + ( uold(1)*axis(2) - uold(2)*axis(1) )*beta

                    p_new(k)%x = new(j,i)%x + unew(1)
                    p_new(k)%y = new(j,i)%y + unew(2)
                    p_new(k)%z = new(j,i)%z + unew(3)
                    r_radius   = p_new(k)%x*p_new(k)%x + p_new(k)%y*p_new(k)%y + p_new(k)%z*p_new(k)%z
                            
                    if (r_sphere_2 > r_radius .or. p_new(k)%z > p_sphere_2 ) then
                        change_1 = 0  ! overlap then try again the pivot of the partial chain on j-th chain
                        exit
                    end if                            
                end do ! pivot trial the partial chain
                        
                if (trial_number>10)then
                    change_1 = 0
                    exit
                end if
            end do  ! end do while change_1
                    
            if (change_1 == 1)then                        
                do k = i+1,Nm_pol     
                    new(j,k)%x = p_new(k)%x 
                    new(j,k)%y = p_new(k)%y 
                    new(j,k)%z = p_new(k)%z  
                end do 
    			      if (i /= 0) then
                	DE1 = DE1 + (new(j,i+1)%x - 2*polymer(j,i)%x + polymer(j,i-1)%x)**2   &
                    	+ (new(j,i+1)%y - 2*polymer(j,i)%y + polymer(j,i-1)%y)**2   &  
	                    + (new(j,i+1)%z - 2*polymer(j,i)%z + polymer(j,i-1)%z)**2   &
    	                - (polymer(j,i+1)%x - 2*polymer(j,i)%x + polymer(j,i-1)%x)**2   &
        	            - (polymer(j,i+1)%y - 2*polymer(j,i)%y + polymer(j,i-1)%y)**2   &  
            	        - (polymer(j,i+1)%z - 2*polymer(j,i)%z + polymer(j,i-1)%z)**2
           
			          end if   !endif i                                    
            else
                exit ! exit the rotate try when over 50 trial pivots
            end if

        end do      ! end of all chains trial move 
    
    end if
        
        if (change_1 == 0)then
            change = 0
        else
            change = 1
        	  DE2 = 0
            DE3 = 0
            do j = 1, n_chain
                do i= 0, Nm_pol 
                    r_radius = dsqrt( new(j,i)%x*new(j,i)%x + new(j,i)%y*new(j,i)%y )
                    ir_temp(j,i) = floor( r_radius / dr ) + 1
                    iz_temp(j,i) = floor( ( p_sphere_2 - new(j,i)%z ) / dz ) + 1  
                    DE2 = DE2 + w(ir_temp(j,i), iz_temp(j,i)) - w(ir(j,i), iz(j,i))
                end do
                DE3 = DE3 + eta(ir_temp(j,Nm_pol), iz_temp(j,Nm_pol)) - eta(ir(j,Nm_pol), iz(j,Nm_pol))
            end do

            r = ran2(seed)
            if ( r < dexp ( -epsilon*DE1 - deltaS*DE2 + DE3 ) )then
                do j = 1, n_chain
                    do i = 0,Nm_pol          
                        polymer(j,i)%x = new(j,i)%x
                        polymer(j,i)%y = new(j,i)%y
                        polymer(j,i)%z = new(j,i)%z
                        iz(j,i) = iz_temp(j,i)
                        ir(j,i) = ir_temp(j,i)
                    end do
                end do            

            else
                change = 0
            end if
            
		end if  ! check metroplis

    end subroutine rotate_sphere 


      FUNCTION ran2(idum)
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1, &
        IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,	 &
        NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
        end do
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=AM*iy
      if(ran2.gt.RNMX) ran2=RNMX
      return
      END function  


subroutine checkpolymer (flag_c)
use global_parameters
implicit none
integer :: i, j, flag_c
double precision :: r

if(i_azo_temp<nm)then
    do j=1,N_azo
        r = (azo(j,i_azo(j))%x - azo(j,i_azo(j)-1)%x)*(azo(j,i_azo(j)+1)%x - azo(j,i_azo(j))%x)    &
           +(azo(j,i_azo(j))%y - azo(j,i_azo(j)-1)%y)*(azo(j,i_azo(j)+1)%y - azo(j,i_azo(j))%y)    &
           +(azo(j,i_azo(j))%z - azo(j,i_azo(j)-1)%z)*(azo(j,i_azo(j)+1)%z - azo(j,i_azo(j))%z) 
        if( abs(r-dcos(2.0d0*pi/3.0d0))>1.0d-6 )then
            print*, j, i_azo(j), "i_azo is broken"    
            print*,azo(j,i_azo(j)-1)%x,azo(j,i_azo(j)-1)%y,azo(j,i_azo(j)-1)%z
            print*,azo(j,i_azo(j))%x,azo(j,i_azo(j))%y,azo(j,i_azo(j))%z
            print*,azo(j,i_azo(j)+1)%x,azo(j,i_azo(j)+1)%y,azo(j,i_azo(j)+1)%z
            stop        
        end if
    end do
end if

do j = 1, N_azo
	do i=1, Nm
		r = (azo(j,i)%x - azo(j,i-1)%x)**2    &
			+(azo(j,i)%y - azo(j,i-1)%y)**2    &
       		+(azo(j,i)%z - azo(j,i-1)%z)**2         
    if (abs(r-1)>1.d-5) then
       print*, "azobond", j,i, "length",abs(r-1)
       print*,azo(j,i)%x,azo(j,i-1)%x
       print*,azo(j,i)%y,azo(j,i-1)%y
       print*,azo(j,i)%z,azo(j,i-1)%z
       flag_c = 1
       stop
    end if
    enddo
end do

do j= 1, N_azo
	do i=1, Nm
    	r = dsqrt(azo(j,i)%z*azo(j,i)%z &
    			+ azo(j,i)%y*azo(j,i)%y &
    			+ azo(j,i)%x*azo(j,i)%x )
        if ( (r-r_sphere) <= -1.0d-5) then
           print*, "azomonoer",i,"on chain", j,"overlap sphere_1"
           print*,azo(j,i)%x
           print*,azo(j,i)%y
           print*,azo(j,i)%z
           flag_c = 1
           stop
        end if

        if ( azo(j,i)%z > p_sphere_2) then
           print*, "azomonoer",i,"on chain", j,"overlap substrate"
           print*,azo(j,i)%x
           print*,azo(j,i)%y
           print*,azo(j,i)%z
           flag_c = 1
           stop
        end if        
    end do
end do


do j = 1, N_chain
	do i=1, Nm_pol
		r = (polymer(j,i)%x - polymer(j,i-1)%x)**2    &
			+(polymer(j,i)%y - polymer(j,i-1)%y)**2    &
       		+(polymer(j,i)%z - polymer(j,i-1)%z)**2         
    if (abs(r-1)>1.d-5) then
       print*, "bond", j,i, "length",abs(r-1)
       print*,polymer(j,i)%x,polymer(j,i-1)%x
       print*,polymer(j,i)%y,polymer(j,i-1)%y
       print*,polymer(j,i)%z,polymer(j,i-1)%z
       flag_c = 1
       stop
    end if
    enddo
end do

do j= 1, N_chain
	do i=1, Nm_pol
    	r = dsqrt(polymer(j,i)%z*polymer(j,i)%z &
    			+ polymer(j,i)%y*polymer(j,i)%y &
    			+ polymer(j,i)%x*polymer(j,i)%x )
        if ( (r-r_sphere) <= -1.0d-5) then
           print*, "monoer",i,"on chain", j,"overlap sphere_1"
           print*,polymer(j,i)%x
           print*,polymer(j,i)%y
           print*,polymer(j,i)%z
           flag_c = 1
           stop
        end if

        if ( polymer(j,i)%z > p_sphere_2) then
           print*, "monoer",i,"on chain", j,"overlap substrate"
           print*,polymer(j,i)%x
           print*,polymer(j,i)%y
           print*,polymer(j,i)%z
           flag_c = 1
           stop
       end if        
    end do
end do
end subroutine checkpolymer
        

             
  SUBROUTINE simpson(g,h,rr)
  IMPLICIT NONE
  Double precision, INTENT(IN) :: h
  Double precision, INTENT(OUT) :: rr
  DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: g
  integer :: n,i
  rr=0.0d0
  n=size(g)

  do i=4,n-3
    rr=rr+g(i)
  end do
  rr=3.0d0*(g(1)+g(n))/8.0d0 &
     +7.0d0*(g(2)+g(n-1))/6.0d0 &
      +23.0d0*(g(3)+g(n-2))/24.0d0 &
       + rr
  rr=rr*h
  END SUBROUTINE simpson 
  

SUBROUTINE comformation_write()
use global_parameters
implicit none
integer :: i,j
open(unit=27,file='azo_ini.txt')
do j=1,n_azo
    do i=0,nm
        write(27,*) azo(j,i)%x, azo(j,i)%y,azo(j,i)%z
    end do
end do
close(27)
open(unit=27,file='polymer_ini.txt')
do j=1,n_azo
    do i=0,nm
        write(27,*) polymer(j,i)%x, polymer(j,i)%y,polymer(j,i)%z
    end do
end do
close(27)
END SUBROUTINE comformation_write

SUBROUTINE comformation_out()
use global_parameters
implicit none
integer :: i, flag_c,j
character*7 res
character res0,res1,res2

!call checkpolymer (flag_c)

do j=1,N_chain
        
    res0=achar(48+mod(j,10))
    res1=achar(48+mod(int(j/10),10))
    res2=achar(48+int(j/100))
    res=res2 // res1 // res0 // '.dat'
    open(4,file=res,status='new')
 
	do i=0,Nm_pol
    	!write(4,"(A4,27X,F7.3,1X,F7.3,1X,F7.3)") "ATOM",polymer(j,i)%x,polymer(j,i)%y,polymer(j,i)%z
    	write(4,*) polymer(j,i)%x,polymer(j,i)%y,polymer(j,i)%z
	end do
	close(4) 
end do
END SUBROUTINE comformation_out
SUBROUTINE comformation_azo_out()
use global_parameters
implicit none
integer :: i, flag_c,j
character*7 res
character res0,res1,res2

!call checkpolymer (flag_c)

do j=1,N_azo
        
    res0=achar(48+mod(j,10))
    res1=achar(48+mod(int(j/10),10))
    res2=achar(48+int(j/100))
    res=res2 // res1 // res0 // '.dat'
    open(4,file=res,status='new')
 
	do i=0,nm
    	!write(4,"(A4,27X,F7.3,1X,F7.3,1X,F7.3)") "ATOM",azo(j,i)%x,azo(j,i)%y,azo(j,i)%z
    	write(4,*) azo(j,i)%x,azo(j,i)%y,azo(j,i)%z
	end do
	close(4) 
end do
END SUBROUTINE comformation_azo_out

END MODULE utility_routines
