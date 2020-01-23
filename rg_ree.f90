program main
    
    real, dimension (:), allocatable :: X, Y, Z, VX, VY, VZ, Rg, Ree, Rg2, Ree2
    integer, dimension (:), allocatable :: NID, NTYP
    real :: rc, v, Rgv, avRg, avRee, comx, comy, comz, enavRg, enavRee

    integer :: I, K, P, NDUMP, NS, NB, NP
    
    !900 format (A)
    !901 format (1X, F9.3, 1X, F9.3, 1X, F9.3, 1X, F9.3)

    NDUMP = 500
    NB = 192000
    NP = 480 
    NS = NB - (NP*20)
    MA = 1
    MB = 1.51 
    
    allocate (Rg(NP), Ree(NP))
    
    open(unit = 39, file = 'rg.txt', action = 'write')
    open(unit = 45, file = 'ree.txt', action = 'write')
    write(39, *) '#Radius of gyration'
    write(39, *)
    write(45, *) '#End to End distance'
    write(45, *)
 
    open (unit = 10, file = 'dump_eq.dil')
    DO 15 P = 1, NDUMP
    print *, 'Evaluating time step: ', P*1000
    
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
       
    IF (P.EQ.1) THEN
        allocate (X(NP), Y(NP), Z(NP), VX(NP), VY(NP), VZ(NP), NID(NP), NTYP(NP))
    END IF
            
    DO 13 I = 1, NP
    comx = 0
    comy = 0
    comz = 0

    DO K = 1, 20
    px = 0
    py = 0
    pz = 0

    read(10, *) NID(K), NTYP(K), X(K), Y(K), Z(K), VX(K), VY(K), VZ(K)
    IF (K.GE.1 .and. K.LE.14) THEN
    px = MA*X(K)
    py = MA*Y(K)
    pz = MA*Z(K)
    ELSE IF (K.GE.15 .and. K.LE.20) THEN
    px = MB*X(K)
    py = MB*Y(K)
    pz = MB*Z(K)
    END IF
   
    comx = comx + px
    comy = comy + py
    comz = comz + pz

    END DO
    comx = comx/(MA*14+MB*6)
    comy = comy/(MA*14+MB*6)
    comz = comz/(MA*14+MB*6)
   
    Rg_i = 0
    DO K = 1, 20
    
    Rgv = 0
    Rgv = (X(K)-comx)**2 + (Y(K)-comy)**2 + (Z(K)-comz)**2
    Rg_i = Rg_i + Rgv
    END DO

    Rg(I) = Rg_i/20
    
    Ree(I) = (X(20)-X(1))**2+(Y(20)-Y(1))**2+(Z(20)-Z(1))**2
    
    avRg = avRg + Rg(I)
    avRee = avRee + Ree(I)

    13 continue

    avRg = avRg/NP
    avRee = avRee/NP
    
    write(39, *) (P*1000), avRg
    write(45, *) (P*1000), avRee


    15 continue
    
    close(10)
    
    close(39)
    close(45)

end program main
