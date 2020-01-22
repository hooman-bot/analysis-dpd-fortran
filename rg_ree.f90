program main
    
    real, dimension (:), allocatable :: X, Y, Z, VX, VY, VZ, Rg, Ree, Rg2, Ree2
    integer, dimension (:), allocatable :: NID, NTYP
    real :: rc, dr, avRg, avRee

    integer :: I, K, P, NDUMP, NS, NB, NP
    
    NDUMP = 1500
    NB = 192000
    NP = 480 
    NS = NB - (NP*20)
    MA = 1
    MB = 1.51 
    
    allocate (Rg(NP), Ree(NP))
    
    open (unit = 10, file = '*insert file_name*')
    DO 15 K = 1, NDUMP
    print *, 'Evaluating time step: ', K*1000
    
    read(10, *)
    read(10, *)
    read(10, *)
    read(10, *) 
    read(10, *)
    read(10, *)
    read(10, *)
    read(10, *)
    read(10, *)
    
    DO I = 1, NS
    read(10, *)
    END DO
       
    IF (K.EQ.1) THEN
        allocate (X(NP), Y(NP), Z(NP), VX(NP), VY(NP), VZ(NP), NID(NP), NTYP(NP))
    END IF
    
        
    DO 13 I = 1, NP
    
    DO K = 1, 20
    read(10, *) NID(K), NTYP(K), X(K), Y(K), Z(K), VX(K), VY(K), VZ(K)
    END DO
    
    Ree(I) = sqrt((X(20)-X(1))**2+(Y(20)-Y(1))**2+(Z(20)-Z(1))**2)
           
    
    
    
end program main
