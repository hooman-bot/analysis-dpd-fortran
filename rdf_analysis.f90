program main

    real, dimension (:), allocatable :: X, Y, Z, VX, VY, VZ, Vsphere, G, Gsol, Gpol
    real, dimension (:), allocatable :: G_k, Gsol_k, Gpol_k
    integer, dimension (:), allocatable :: NID, NTYP,  Nsam, Nsamsol, Nsampol
    real, dimension (:,:), allocatable :: G_k, Gsol_k, Gpol_k
    real, dimension (:), allocatable :: bin
    real :: rcut, dr, Nshell

    integer :: I, K, P, NDUMP, NS, NB, NP


    900 format (A)
    901 format (1X, F9.3, 1X, F9.3, 1X, F9.3, 1X, F9.3)


    NDUMP = 1500
    NB = 192000
    rcut = 1
    dr = 0.025*rcut

    Nmax = int((3*rcut/dr) + 0.5)

    allocate(bin(Nmax))

    DO I =1, Nmax
    bin(I) = dr*(I-1)
    END DO

    ! intialising the array
    allocate (G(Nmax), Gsol(Nmax), Vsphere(Nmax), Nsam(Nmax), Nsamsol(Nmax), Nsampol(Nmax))
    allocate (G_k(Nmax,NDUMP),Gsol_k(Nmax,NDUMP),Gpol_k(Nmax,NDUMP))


    DO 16 K = 1, 1

    open (unit = 10, file = 'dump_eq_restart.dil')

    read(10, *)
    read(10, *)
    read(10, *)
    read(10, *) NB
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
    ! make an array of the bins, for each value in the bin, loop through
    ! if rij is between r & r+dr, add 1 to the count.

    !DO 20 P = 1, Nmax 

    Nsam(:) = 0
    Nsamsol(:) = 0
    Nsampol(:) = 0

    DO 19 I = 1, NB

    DO J=(I+1), NB-1

    RIJ = sqrt((X(I)-X(J))**2+(Y(I)-Y(J))**2+(Z(I)-Z(J))**2)

    Nshell = int((RIJ/dr) + 0.5)

    IF (Nshell.LE.Nmax) THEN
        Nsam(Nshell) = Nsam(Nshell)+1
        IF (I.GE.NS .and. J.GE.NS) THEN
            Nsampol(Nshell) = Nsampol(Nshell) + 1
        ELSE IF (I.LE.NS .and. J.LE.NS) THEN
            Nsam(Nshell) = Nsamsol(Nshell) + 1
        END IF
    END IF

    print *, Nsam

    END DO

    19 continue
    
    DO 21 Nshell = 1, Nmax

    Vsphere(Nshell) = 4*pi*dr*(dr*Nshell)**2+(pi*dr**3)/3
    G_k(Nshell,K) = (Nsam(Nshell))/(Vsphere(Nshell))
    Gsol_k(Nshell,K) = (Nsamsol(Nshell))/(Vsphere(Nshell))
    Gpol_k(Nshell,K) = (Nsampol(Nshell))/(Vsphere(Nshell))

    21 continue
    16 continue

    G  = abs(2*BOXL**3*sum(sum(G_k, DIM = 2), DIM = 2)/(NDUMP*NPRT**2))
    Gsol  = abs(2*BOXL**3*sum(sum(Gsol_k, DIM = 2), DIM = 2)/(NDUMP*NSOL**2))
    Gpol  = abs(2*BOXL**3*sum(sum(Gpol_k, DIM = 2), DIM = 2)/(NDUMP*NPOL**2))

    open(unit = 200, file = 'GR.txt', STATUS='REPLACE', ACTION='WRITE' )

    write(200,900) '#Simulation results: Radial Distribuion'
    write(200,900)

    DO 155 I = 1, Nshell
    write(200,901) (I*dr/rcut), G(I), Gsol(I), Gpol(I)
    155   CONTINUE
    close(200, STATUS='KEEP' )

    close(10)

end program main
