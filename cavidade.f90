!Pressure-Velocity coupling solution for the Navier-Stokes simulation in a bidimensional domain
!Undergraduate student: Rodrigo Daher
!Professor: Aristeu da Silveira Neto 

!MODULE TO DEFINE THE GLOBAL PARAMETERS NEEDED TO SIMULATION
module global
    implicit none 

    !!!========PhysicaL Parameters==========!!!

    real*16 :: A_x, D_x, A_y, D_y, P_x, P_y, rho, mi, gama, V_lid, teste
    double precision, dimension(20,20) :: P, P_linha
    double precision, dimension(20,19) :: u, u_int
    double precision, dimension(19,20) :: v, v_int 

    !!!=======Computational Parameters======!!!
    integer :: i, j, it, itf, m, pressure_step
    double precision :: dX, dY, dT, tolerance, CFL, nu, Lx, Ly
    double precision, dimension(20,20) :: Erro, Div 

end module global




!MODULE TO DEFINE THE PARAMETERS NEEDED TO WRITE THE DATA
module file 
    implicit none 
    
    !!!----------Files Parameters---------!!!
    character(len=70) :: fn, command
    integer, parameter :: outunit=400   

end module file




!MAIN PROGRAM
program cavidade
    use global
    use file    
    implicit none

    !Start simulation loop 

    call initial_conditions()
   ! call boundary_conditions()
    call write_file()

    print*, '!-------------------------------!'
    print*, '!------BEGIN OF SIMULATION------!'
    print*, '!-------------------------------!'
    print*, ''
    print*, 'Solving pressure-velocity coupling by Fractional Step Method'
    print*, 'The mesh used is an Unstaggered grid'
    print*, ''
    print*, 'Parameters of simulation:'
    print*, 'dt = ', dt
    print*, 'dx = ', dX 
    print*, 'dy = ', dY
    print*, 'Shape (U) = ', shape(U)
    print*, 'Shape (V) = ', shape(V)
    print*, 'Shape (P) = ', shape(P)
    print*, ''

    do while (it < itf)

        !call boundary_conditions()
        call velocity_prediction()
        call pressure_correction()
        call velocity_correction()
        call boundary_conditions()
        call continuity_equation()
        call pressure_update()
        call boundary_conditions()

        it = it + 1
        call write_file()
        print*, it
    end do 

    print*, '!-------------------------------!'
    print*, '!-------END OF SIMULATION-------!'
    print*, '!-------------------------------!'

end program cavidade




!SET INITIAL CONDITIONS OF THE SIMULATION
subroutine initial_conditions()
    use global
    implicit none
    
    !Initial Parameters!
    it = 0
    itf = 100                                                                                         !Number of iterations
    m  = 20                                                                                           !Dimension of the pressure matrix MxM
    Lx = 1.0d0
    Ly = 1.0d0
    dX = Lx/m                                                                                         !X Cell length 
    dY = Ly/m                                                                                         !Y Cell length
    CFL = 0.025d0                                                                                     !Relation between time and space 
    mi = 1.0533d-6                                                                                    !Fluid viscosity
    rho = 20.0d0                                                                                      !Fluid density  
    nu = mi/rho
    dT = 1.0d-4                                                                                       !Time step
    v  = 0.0d0                                                                                        !Inicial velocity in y direction 
    V_lid = 1.0d0                                                                                     !Lid velocity
    P  = 1.0d0                                                                                        !Inicial pressure
    u_int = 0.0d0                                                                                     !Intermediate velocity in x direction
    v_int = 0.0d0                                                                                     !Intermediate velocity in y direction
    P_linha = 0.0d0                                                                                   !Matrix of correction of pressure 
    tolerance = 1.d-10                                                                                !Increment for implicity Gaus-Seidel solutions
    Erro = 0.0d0
    do i = 1, m                                                                                       !Inicial velocity in x direction
        do j = 1, m-1                                
            if (i == m-1 .or. i == m) then
                u(i , j) = V_lid
            else 
                u(i , j) = 0.0d0
            end if
        end do 
    end do

end subroutine initial_conditions




!UPDATE GHOST CELLS
subroutine boundary_conditions()
    use global
    implicit none 

    do i = 1, m
        u(i , 1) = 0.0d0                                                                 !U velocity on the left wall
        u_int(i , 1) = 0.0d0                                                             !U* velocity on the left wall
        u(i , m-1) = 0.0d0                                                               !U velocity on the right wall
        u_int(i , m-1) = 0.0d0                                                           !U* velocity on the right wall

        v(1 , i) = 0.0d0                                                                 !V velocity on the up wall
        v_int(1 , i) = 0.0d0                                                             !v* velocity on the up wall
        v(m-1 , i) = 0.0d0                                                               !V velocity on the down wall
        v_int(m-1 , i) = 0.0d0                                                           !v* velocity on the down wall

        P_linha(1 , i) =  P_linha(2, i)                                                  !Pressure on the up wall
        P_linha(m , i) =  P_linha(m-1, i)                                                !Pressure on the down wall    
        P_linha(i , 1) =  P_linha(i, 2)                                                  !Pressure on the left wall
        P_linha(i , m) =  P_linha(i, m-1)                                                !Pressure on the right wall

        P(1 , i) =  P(2, i)                                                              !Pressure on the up wall
        P(m , i) =  P(m-1, i)                                                            !Pressure on the down wall    
        P(i , 1) =  P(i, 2)                                                              !Pressure on the left wall
        P(i , m) =  P(i, m-1)                                                            !Pressure on the right wall
    end do 

    do i = 1, m
        do j = 1, m-1 
            u(m , i) = 2 * V_lid - u(m-1 , j)                                            !U Velocity on the up wall
            u_int(m , i) = 2 * V_lid - u_int(m-1 , j)                                    !U* velocity on the up wall
        end do 
    end do 

    do i = 1, m-1
   
        u(1 , i) = - u(2 , i)                                                            !U velocity on the down wall
        u_int(1 , i) = - u_int(2 , i)                                                    !U* Velocity on the down wall

        v(i , 1) = - v(i , 2)                                                            !V velocity on the left wall
        v_int(i , 1) = - v(i , 2)                                                        !V* velocity on the left wall
        v(i , m) = - v(i , m-1)                                                          !V velocity on the right wall
        v_int(i , m) = - v_int(i , m-1)                                                  !V* Velocity on the right wall
    end do 

end subroutine boundary_conditions




!VELOCITY PREDICITION
subroutine velocity_prediction()
    use global    
    implicit none 

    !-----------------------------STEP 1: PREDICT U*, V*-----------------------------!

    do i = 2, m-1
        do j = 2, m-2
            A_x = - u(i, j) * (u(i+1, j) - u(i,j))/dX - (v(i, j) + v(i-1, j) + v(i, j+1) + v(i-1, j+1))/(4 * dX) &          !Advective term in the x direction
                * (u(i, j+1) - u(i, j)) 

            D_x = (u(i+1, j) + u(i-1, j) + u(i, j+1) + u(i, j-1) - 4 * u(i, j)) / (dX**2)                                   !Difusive term in the x direction

            P_x = (P(i, j) - P(i-1, j)) / dX                                                                                !Pressure Gradient in the x direction                                                                              

            u_int(i,j) = u(i, j) + dT * (A_x + nu * D_x  - (1/rho) * P_x)                                                   !U prediction  
        end do 
    end do

    call boundary_conditions()

    do i = 2, m-2
        do j = 2, m-1
            A_y = -(u(i, j) + u(i, j-1) + u(i+1, j) + u(i+1, j-1))/(4*dX) * (v(i+1, j) - v(i, j)) &                         !Advective term in the y direction
                - v(i, j) *(v(i, j+1) - v(i,j))/dY

            D_y = (v(i+1, j) + v(i-1, j) + v(i, j+1) + v(i, j-1) - 4 * v(i, j)) / (dY**2)                                   !Difusive term in the y direction

            P_y = (P(i, j) - P(i, j-1)) / dY                                                                                !Pressure Gradient in the y direction

            v_int(i, j) = v(i, j) + dT * (A_y + nu * D_y - (1/rho) * P_y)                                                   !V prediction   
        end do 
    end do

    call boundary_conditions() 

end subroutine velocity_prediction




!PRESSURE CORRECTION
subroutine pressure_correction()
    use global    
    implicit none 

    !----------------------------STEP 2: PRESSURE CORRECTION----------------------------!

    gama = (-rho * dX) / dT

    pressure_step = 0

    P_linha = 0.0d0
    Erro(10,10) = 1.0d0

    !-----------------------------Gaus-Seidel Method----------------------------!
    do while (MAXVAL(Erro) > tolerance)
        do i = 2, m-1
            do j = 2, m-1
                Erro(i , j) =  P_linha(i , j)
            end do 
        end do

       !P_linha(m/2 , m/2) = 0

        do i = 2, m-1
            do j = 2, m-1
                P_linha(i, j) = (gama * (u_int(i+1 , j) - u_int(i , j) + v_int(i , j+1) - v_int(i , j)) + &
                (P_linha(i-1 , j) + P_linha(i+1 , j) + P_linha(i , j-1) + P_linha(i , j+1))) /4
            end do 
        end do 

        !call boundary_conditions()

        do i = 2, m-1
            do j = 2, m-1
                Erro(i , j) = abs(Erro(i , j) - P_linha(i , j))
            end do 
        end do


        pressure_step = pressure_step + 1
        print*, maxval(p_linha), pressure_step
    end do

end subroutine pressure_correction




!VELOCITY CORRECTION
subroutine velocity_correction()
    use global
    implicit none

    !-----STEP 3: CORRECT U, V------!

    !====== X velocity correction ======!
    do i = 2, m-1
        do j = 2, m-2 
            u(i , j) = u_int(i , j) - (dT/rho) * (P_linha(i , j) - P_linha(i-1 , j))/dX                         !Correct the U velocity based on the pressure correction
        end do
    end do

    !======= Y Velocity correction ======!
    do i = 2, m-2
        do j = 2, m-1
            v(i , j) = v_int(i , j) - (dT/rho) * (P_linha(i , j) - P_linha(i , j-1))/dY                         !Correct the V velocity based on the pressure correction
        end do
    end do

end subroutine velocity_correction




!PRESSURE UPDATE
subroutine pressure_update()
    use global
    implicit none

    !===== Pressure update =====!
    do i = 1, m
        do j = 1, m
            P(i , j) = P(i , j) + P_linha(i , j)                                   !Update the pressure based on the pressure gradient correction on step 2
        end do 
    end do 

end subroutine pressure_update




!DIVERGENT CALCULATION
subroutine continuity_equation()
    use global    
    implicit none

    !===== Check if the continuity equation is OK =====!
    do i = 2, m-2
        do j = 2, m-2
            Div(i , j) = abs((u(i+1 , j) - u(i , j))/dX + (v(i , j) - v(i , j+1))/dY)                         !Calculates the divergent of the velocity 
        end do 
    end do 

    print*, maxval(Div)                                                                                 !Shows the divergent, is expected a number close to 0

end subroutine continuity_equation




!WRITE FILE WITH DATA OBTAINED IN THE simulation
subroutine write_file()
    use global
    use file 
    implicit none 

    !======Write the file.dat ====!
    !if (mod(it, 10) == 0) then
        write(fn,fmt='(i0,a)') it, '.dat'
        open(unit = outunit, file=fn, form='formatted')

        !Data write
        !write(outunit,*) "TITLE = " , '"Cavity_simu"'
        !write(outunit,*) 'Variables="X","Y","U"'
        !write(outunit,*) 'Zone I=', Lx ,', J=', Ly ,', F=POINT'

        do i = 2, m-1
            do j = 2, m-1
                write(outunit, *) j , i, P(i , j)
            end do 
        end do                 
        close(outunit)
    !end if
    !if (mod(it,10) == 0 ) then
    !    command='py plot.py '// fn
    !    call EXECUTE_COMMAND_LINE(command)
    !end if

end subroutine write_file