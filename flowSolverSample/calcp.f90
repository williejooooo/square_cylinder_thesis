SUBROUTINE CALCP
USE DATA
IMPLICIT NONE
REAL :: DVA, DUA
    RESORM = 0.0

  !---assembly of coefficients-----------------------------------------

    DO I = 2, NIM1
      DO J = 2, NJM1

        IP1 = I + 1
        IM1 = I - 1
        JP1 = J + 1
        JM1 = J - 1

        ! compute areas and volume

        AREAN = RV(J + 1) * SEW(I)
        AREAS = RV(J) * SEW(I)
        AREAEW = R(J) * SNS(J)
        VOL = R(J) * SNS(J) * SEW(I)

        ! calculate coefficients

        DENN = 0.5*(DEN(I, J) + DEN(I, JP1))
        DENS = 0.5 * (DEN(I, J) + DEN(I, JM1))
        DENE = 0.5 * (DEN(I, J) + DEN(IP1, J))
        DENW = 0.5 * (DEN(I, J) + DEN(IM1, J))
        AN(I,J) = DENN * AREAN * DV(I, JP1)
        AS(I,J) = DENS * AREAS * DV(I, J)
        AE(I,J) = DENE * AREAEW * DU(IP1, J)
        AW(I,J) = DENW * AREAEW * DU(I, J)
!       DVA     = .5*(DV(I,JP1) +DV(I,J))
!       DUA     = .5*(DU(IP1,J) +DU(I,J))
!       AN(I,J) = DENN * AREAN * DVA
!       AS(I,J) = DENS * AREAS * DVA
!       AE(I,J) = DENE * AREAEW * DUA
!       AW(I,J) = DENW * AREAEW * DUA

        ! calculate source terms

        CN = DENN * V(I, JP1) * AREAN
        CS = DENS * V(I, J) * AREAS
        CE = DENE * U(IP1, J) * AREAEW
        CW = DENW * U(I, J) * AREAEW
        SMP = CN - CS + CE - CW
        SP(I, J) = 0.0
        SU(I, J) = -SMP

        ! compute sum of absolute mass sources

      END DO
    END DO

  !---problem modifications---------------------------------------------

!***************
    CALL MODP
!***************

  !---final coefficient assembly---------------------------------------

    DO I = 2, NIM1
      DO J = 2, NJM1
        RESORM = RESORM + ABS(SU(I, J))
        AP(I, J) = AN(I, J) + AS(I, J) + AE(I, J) + AW(I, J) - SP(I, J)
      END DO
    END DO

    RESORM = RESORM / Amflin

  !---solution of difference equations----------------------------------

    DO N = 1, NSWPP
!**********************************
      CALL SOLVEX(2,2,mx,my,PP)
      CALL LISOLV(2,2,mx,my,PP)
      CALL SOLVEX(2,2,mx,my,PP)
!**********************************
    END DO

  !---correct velocities and pressures----------------------------------

    ! velocities

    DO I = 2, NIM1
      DO J = 2, NJM1
        IM1 = I - 1
        JM1 = J - 1
!       IF(I /= 2) U(I, J) = U(I, J) + DU(I, J) * (PP(IM1, J) - PP(I, J))
!       IF(J /= 2) V(I, J) = V(I, J) + DV(I, J) * (PP(I, JM1) - PP(I, J))
        U(I, J) = U(I, J) + DU(I, J) * (PP(IM1, J) - PP(I, J))
        V(I, J) = V(I, J) + DV(I, J) * (PP(I, JM1) - PP(I, J))
      END DO
    END DO

    ! pressures (with provision for under-relaxation

    PPREF = PP(IPREF, JPREF)
!   DO J=2,NJM1
!    PP(IPREF,J) = PPREF
!   ENDDO
    DO I = 2, NIM1
      DO J = 2, NJM1
        P(I, J) = P(I, J) + URFP * (PP(I, J) - PPREF)
        PP(I, J) = 0.0
      END DO
    END DO
    DO J=2,NJM1
      U(NI,J) = U(NIM1,J)
!     U(2,J)  = U(3,J)
!     U(1,J)  = U(2,J)
    ENDDO
!   do not correct the velocities on the block, set them all to zero
!   blank out u on and inside block
    DO I = IB1,IB2
    DO J = JB1,JB2-1
       U(I, J)  = 0.0
    END DO
    END DO
!   blank out v on and inside block
    DO I = IB1,IB2-1
    DO J = JB1,JB2
       V(I, J)  = 0.
    END DO
    END DO
!### Set pp=0 and p=0 inside block just to mark the block
    DO I = IB1,IB2-1
    DO J = JB1,JB2-1
       PP(I, J)  = 0.0
       P(I, J)   = 0.0
    END DO
    END DO
END SUBROUTINE CALCP
