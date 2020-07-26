SUBROUTINE INIT
!### subroutine to build the grid and the cells for U, V, P
!###  F.L. 5/2020 and initialize the flow field including the proper boundary conditions.
!###
USE DATA
IMPLICIT NONE
 
    NIM1 = NI - 1
    NJM1 = NJ - 1

    NWRITE = 0
    GREAT = 1.0E30   ! large constant
    R_ONLY = .TRUE.

    DO J = 1, NJ
      R(J) = Y(J)
      IF(INDCOS == 1) R(J) = 1.0 ! INDCOS == 1 -> depth is constant (no pipe)
    END DO

!---distance between nodal points in x-direction------------------------

    DXPW(1) = 0.0
    DXEP(NI) = 0.0
    DO I = 1, NIM1
      DXEP(I)   = X(I + 1) - X(I)
      DXPW(I+1) = DXEP(I)
    END DO

!---distance between nodal points in y-direction------------------------

    DYPS(1) = 0.0
    DYNP(NJ) = 0.0
    DO J = 1, NJM1
      DYNP(J) = Y(J + 1) - Y(J)
      DYPS(J + 1) = DYNP(J)
    END DO

!---create p cell-------------------------------------------------------

    ! width of p-cell (= width of v-cell)

    SEW(1) = 0.0
    SEW(NI) = 0.0
    DO I = 2, NIM1
      SEW(I) = 0.5 * (DXEP(I) + DXPW(I))
    END DO
    SEW(2) = 0.5 * (DXEP(2) + 2. * DXPW(2))
    SEW(NIM1) = 0.5 * (2. * DXEP(NIM1) + DXPW(NIM1))

    ! height of p-cell (= height of u-cell)

    SNS(1) = 0.0
    SNS(NJ) = 0.0
    DO J = 2, NJM1
      SNS(J) = 0.5 * (DYNP(J) + DYPS(J))
    END DO
    SNS(2) = 0.5 * (DYNP(2) + 2. * DYPS(2))
    SNS(NJM1) = 0.5 * (2. * DYNP(NJM1) + DYPS(NJM1))

!---create u cell-------------------------------------------------------

    ! new coordinate XU as location for u-velocity

     XU(1) = 0.0
     XU(2) = 0.0
     XU(NI) = X(NI)
    DO I = 3, NIM1
      XU(I) = 0.5*(X(I-1) +X(I))
    END DO
!   XU(1) = XU(2) -(XU(3) -XU(2))

    ! distance between loction of u velocity in x-direction

    DXPWU(1) = 0.0
    DXPWU(2) = 0.0
    DXEPU(1) = 0.0
    DXEPU(NI) = 0.0
    DO I = 2, NIM1
      DXEPU(I)   = XU(I+1) -XU(I)
      DXPWU(I+1) = DXEPU(I)
    END DO
!    DXEPU(NI) = DXPWU(NI)

    ! width of u-cell

    SEWU(1) = 0.0
    SEWU(2) = 0.0
    SEWU(NI) = 0.0
    DO I = 3, NIM1
      SEWU(I) = X(I) - X(I - 1)
    END DO

!---create v cell-------------------------------------------------------

    ! new coordinate YV as location for v-velocity

    YV(1) = 0.0
    YV(2) = 0.0
    DO J = 3, NJM1
      YV(J)=0.5 * (Y(J) + Y(J - 1))
    END DO
    YV(NJ) = Y(NJ)

    ! RV : radial coordinate in middle of v-cell

    RV(1) = 0.0
    RV(2) = 1.0
    DO J = 3, NJM1
      RV(J) = 0.5 * (R(J) + R(J-1))
    END DO
    RV(NJ) = R(NJ)

    ! RCV : radial coordinate at the bottem of v-cell

    RCV(NJ) = R(NJ)
    RCV(2) = 1.0
    DO J = 3, NJM1
      RCV(J) = R(J - 1)
    END DO
    RCV(3) = 1.0

    ! distance between locations for v-velocity (YV) in y-direction

    DYPSV(1) = 0.0
    DYPSV(2) = 0.0
    DYNPV(1) = 0.0
    DYNPV(NJ) = 0.0
    DO J = 2, NJM1
      DYNPV(J) = YV(J + 1) - YV(J)
      DYPSV(J + 1) = DYNPV(J)
    END DO

    ! height of v-cell

    SNSV(1) = 0.0
    SNSV(2) = 0.0
    SNSV(NJ) = 0.0
    DO J = 3, NJM1
      SNSV(J) = Y(J) - Y(J - 1)
    END DO

!---ratio of U or V cell distance to nodal distance---------------------

    ! x-direction

    DXP(2) = 0.0
    DXM(2) = 1.0
    DO I = 3, NIM1
      DXP(I) = SEWU(I) / (2. * DXEPU(I))
      DXM(I) = 1. - DXP(I)
    END DO
    DXP(NIM1) = 1.0
    DXM(NIM1) = 0.0

    SEWU(3) = SEWU(3) + DXPW(2)        ! bigger cell at boundary
    SEWU(NIM1) = SEWU(NIM1) + DXPW(NI)

    ! y-direction

    DYP(2) = 0.0
    DYM(2) = 1.0
    DO J = 3, NJM1
      DYP(J) = SNSV(J) / (2. * DYNPV(J))
      DYM(J) = 1.0 - DYP(J)
    END DO
    DYP(NJM1) = 1.0
    DYM(NJM1) = 0.0

    SNSV(3) = SNSV(3) + DYPS(2)        ! bigger cell at boundary
    SNSV(NJM1) = SNSV(NJM1) + DYPS(NJ)

!---set variables to zero-----------------------------------------------

    NIM2 = NIM1 - 1
    NJM2 = NJM1 - 1

        U    = 0.0
        V    = 0.0
        P    = 1.0
        BFX  = 0.
        BFY  = 0.
        PN   = 1.0
        PP   = 0.0
        PEE  = 0.0
        DEN  = DENSIT
        VIS  = VISCOS
        GAMH = VISCOS/PRANDT
        DU   = 0.0
        DV   = 0.0
        AN   = 0.0
        AS   = 0.0
        AE   = 0.0
        AW   = 0.0
        SU   = 0.0
        SP   = 0.0
        T    = 0.0
        UI   = 0.0
        VI   = 0.0
!###
!### initialize variable fields
!###
    DO I = 2, NIM1       ! temperature top and bottom boundary
      T(I, NJ) = TTOP
      T(I, 1) = TBOT
    END DO
    DO J = 2, NJM1       ! temperature left and right boundary
      T(1, J) = TLEFT
      T(NI, J) = TRITE
    END DO
!**    DO J = 2, JSTEP      ! temperature at wall behind step
!**      T(1, J) = TSTEP
!**    END DO
!###
!### specify inlet velocity profile and make all initial flow field the same
!###
  DO I=1,NI
  DO J=1,NJ
    U(I,J) = VELOC*(1. - (2.0*Y(J)/H)**2.)
    V(I,J) = 0.
  ENDDO
  ENDDO
  UN = U
  VN = V
  RSUN = 0.
  RSVN = 0.
  RSU  = 0.
  RSV  = 0.
  I = NI/2
  J = NJ/2
  DT  = CFL* SEW(I)*SNS(J)/SQRT((U(I,J)*SNS(J))**2 +(V(I,J)*SEW(I))**2)
! write(6,*) CFL,DT

! Set Quantities for Normalization of Residuals
    TIN = 1.0
    AMFLIN = DENSIT * (Y(NJ) -Y(1))*VELOC
    AMOMIN = AMFLIN * VELOC
    TFLOIN = AMFLIN* TIN
! Calculate Reynolds number
!   RYNDS = DENSIT * VELOC * 2.0 * RADIUS / VISCOS    ! Reynolds number
!
!    OPEN(UNIT = 11, FILE='mon.out', FORM='FORMATTED', ACCESS = 'SEQUENTIAL', STATUS = 'UNKNOWN')
!    write(11,1000)
!    DO I=1,NI
!      WRITE(11,1100) I, X(I),XU(I),DXPW(I),DXEP(I),SEW(I),DXPWU(I),DXEPU(I),SEWU(I)
!    ENDDO
!    write(11,2000)
!    DO J=1,NJ
!      WRITE(11,2100) J, Y(J),YV(J),DYPS(J),DYNP(J),SNS(J),DYPSV(J),DYNPV(J),SNSV(J)
!    ENDDO
!    CLOSE(11)
1000 FORMAT('I', 13x,  'X', 9x, 'XU', 8x, 'DXPW', 6x, 'DXEP', 6x, 'SEW', 6x, 'DXPWU',6x,'DXEPU',6x,'SEWU')
1100 FORMAT(I2,8X  8F10.4)
2000 FORMAT('J', 13x,  'Y', 9x, 'YV', 8x, 'DYPS', 6x, 'DYNP', 6x, 'SNS', 6x, 'DYPSV',6x,'DXNPV',6x,'SNSV')
2100 FORMAT(I2,8X  8F10.4)
END SUBROUTINE INIT
