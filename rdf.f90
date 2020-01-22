program main
    ! Calculating the radial distribution function of the solvent
    real, dimension (:), allocatable :: X, Y, Z, VX, VY, VZ, Vol, unG, unGsol, unGpol, G, Gsol, Gpol, Gn, Gsoln, Gpoln
    integer, dimension (:), allocatable :: NID, NTYP
    real :: rc, dr

    integer :: I, K, P, NDUMP, NS, NB, NP, bn, kn

    900 format (A)
    901 format (1X, F9.3, 1X, F9.3, 1X, F9.3, 1X, F9.3)

    REAL, PARAMETER :: pi = 3.1415927
    
    NDUMP = 1500
    NB = 192000
    L = 40
    rc = 1
    dr = 0.025*rc
    
    bn = int((3*rc/dr))

    allocate (Vol(bn), unG(bn), unGsol(bn), unGpol(bn), G(bn), Gsol(bn), Gpol(bn), Gn(bn), Gsoln(bn), Gpoln(bn))
    
    Gn(:) = 0
    Gsoln(:) = 0
    Gpoln(:) = 0

    open (unit = 10, file = '*insert file_name*')
    
    DO 16 K = 1, 1
    
    print *, 'Evaluating time step: ', K*1000
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

    DO kn =1, bn

    Vol(kn) = 4*pi*dr*(dr*kn)**2+(pi*dr**3)/3
    G(kn) = abs(((L**3)*unG(kn))/(Vol(kn)*NB**2))
    Gsol(kn) =abs(((L**3)*unGsol(kn))/(Vol(kn)*NB**2))
    Gpol(kn) = abs(((L**3)*unGpol(kn))/(Vol(kn)*NB**2))

    END DO

    Gn = Gn + G
    Gsoln= Gsoln + Gsol
    Gpoln = Gpoln + Gpol

    close(10)

    16 continue
    
    Gn = Gn/NDUMP
    Gsoln = Gsoln/NDUMP
    Gpoln = Gpoln/NDUMP
    
    open(unit = 39, file = 'gr_ndump_all.txt', action = 'write')
    open(unit = 45, file = 'gr_ndump_sol.txt', action = 'write')
    open(unit = 67, file = 'gr_ndump_pol.txt', action = 'write')

    write(39, 900) '#Test gr results: All'
    write(39, 900)

    write(45, 900) '#Test gr results: sol'
    write(45, 900)

    write(67, 900) '#Test gr results: pol'
    write(67, 900)

    DO I = 1, bn

    write(39, 901) (I*dr), Gn(I)

    write(45, 901) (I*dr), Gsoln(I)

    write(67, 901) (I*dr), Gpoln(I)
    
    END DO

    close(39)
    close(45)
    close(67)

end program main
