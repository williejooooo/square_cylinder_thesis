SUBROUTINE GEOM
!### Set up geometry and control volumes
USE DATA
IMPLICIT NONE

    DX = W / FLOAT(NI - 2)  ! distance between grid points
    DX = XLEN / FLOAT(IB2-IB1)  ! distance between grid points
    DY = H / FLOAT(NJ - 2)
!  Original from SEE
!    X(1)=0.0
!      DO I=2,NI
!      X(I)=X(I-1)+DX
!    END DO
!
!    Y(1)=0.0
!    DO J=2,NJ
!      Y(J)=Y(J-1)+DY
!    END DO
!
!  FL Change: first and last cell will be half size of the interior.
!
    DO I=1,NI
      X(I) = (I-1.5)*DX
    END DO
    DO J=1,NJ
      Y(J) = (J-1.5)*DY
    END DO
    DO J=1,NJ    ! set y=0 at center of channel
      Y(J)  = Y(J) -H/2.
    ENDDO
    RADIUS = Y(NJ)   ! radius of the pipe
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

    RETURN
END SUBROUTINE GEOM
