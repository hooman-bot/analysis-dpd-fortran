program main
    real, dimension (:), allocatable :: X, Y, Z, VX, VY, VZ, Vol, unG, unGsol, unGpol, G, Gsol, Gpol
    integer, dimension (:), allocatable :: NID, NTYP
    real :: rc, dr

    integer :: I, K, P, NDUMP, NS, NB, NP, bn, kn


    900 format (A)
    901 format (1X, F9.3, 1X, F9.3, 1X, F9.3, 1X, F9.3)

    REAL, PARAMETER :: pi = 3.1415927
    REAL, PARAMETER :: Boltz = 1.38064852
    REAL, PARAMETER :: Avo = 6.0221409


    NDUMP = 1500
    NB = 192000
    !NB = 50
    L = 40
    rc = 1
    dr = 0.025*rc



    bn = int((3*rc/dr))

    allocate (Vol(bn), unG(bn), unGsol(bn), unGpol(bn), G(bn), Gsol(bn), Gpol(bn))

    DO 16 K = 1, 1

    open (unit = 10, file = 'dump_eq_restart.dil')

    read(10, *)
    read(10, *)
    read(10, *)
    read(10, *) !NB
    read(10, *)
    read(10, *)
    read(10, *)
    read(10, *)
    read(10, *)

    NP = 480
    NS = NB - NP
    IF (K.EQ.1) THEN
        allocate (X(NB), Y(NB), Z(NB), VX(NB), VY(NB), VZ(NB), NID(NB), NTYP(NB))
    END IF

    DO I = 1, NB
    
    read(10, *) NID(I), NTYP(I), X(I), Y(I), Z(I), VX(I), VY(I), VZ(I)

    END DO

    DO 19 I = 1, NB-1
    DO J=I+1, NB

    RIJ = sqrt((X(I)-X(J))**2+(Y(I)-Y(J))**2+(Z(I)-Z(J))**2)
    kn = int(RIJ/dr)

    IF (kn.LE.bn) THEN
        unG(kn) = unG(kn) + 1
        IF (I.GE.NS .and. J.GE.NS) THEN
            unGpol(kn) = unGpol(kn) + 1
        ELSE IF (I.LE.NS .and. J.LE.NS) THEN
            unGsol(kn) = unGsol(kn) + 1
        END IF
    END IF

    END DO

    19 continue

    DO 20 kn =1, bn

    Vol(kn) = 4*pi*dr*(dr*kn)**2+(pi*dr**3)/3
    G(kn) = abs(((L**3)*unG(kn))/(Vol(kn)*NB**2))
    Gsol(kn) =abs(((L**3)*unGsol(kn))/(Vol(kn)*NB**2))
    Gpol(kn) = abs(((L**3)*unGpol(kn))/(Vol(kn)*NB**2))

    !print *, Vol(kn), G(kn), Gsol(kn), Gpol(kn)
    20 continue

    open(unit = 15, file = 'gr_test_all.txt', action = 'write')

    open(unit = 23, file = 'gr_test_sol.txt', action = 'write')

    open(unit = 27, file = 'gr_test_pol.txt', action = 'write')

    write(15, 900) '#Test gr results: All'
    write(15, 900)

    write(23, 900) '#Test gr results: SOL'
    write(23, 900)

    write(27, 900) '#Test gr results: POL'
    write(27, 900)


    DO 25 I = 1, bn

    write(15, 901) (I*dr), G(I)
    write(23, 901) (I*dr), Gsol(I)
    write(27, 901) (I*dr), Gpol(I)

    25 continue

    close(15)

    close(10)

    16 continue 
end program main
