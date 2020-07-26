SUBROUTINE PROMOD
USE DATA
IMPLICIT NONE
REAL :: RHOA, SUM1, SUM2, DELU, AMAS
REAL :: AD


  !---properities-------------------------------------------------------

    !***********
    ENTRY MODPRO    ! no modifications for this problem
    !***********

    RETURN

  !---U momentum--------------------------------------------------------

    !*********
    ENTRY MODU
    !*********

    DO I = 2, NIM1
      U(I, 1)  = 0.0  !lower wall
      U(I, NJ) = 0.0  !upper wall
    END DO

!    DO I = 2, NIM1     ! symmtry plane at lower boundary
!      AS(I, 2) = 0.0
!    END DO
!   DO J = 2, NJM1
!     U(NI, J)   = U(NIM1, J)
!   END DO
!
!   blank out block
!
    DO I = IB1,IB2
    DO J = JB1,JB2-1
       SU(I, J) = 0.0
       AP(I, J) = -GREAT
       U(I, J)  = 0.0
       RSU(I,J) = 0.0
    END DO
    END DO
!
    DO J=2,NJM1
      U(NI,J) = U(NIM1,J)
    ENDDO
!
! Convective outflow boundary condition
!
!   AD = DT*VELOC/(X(NI) -X(NIM1))
!   DO J=2, NJM1
!     U(NI,J)  = (UN(NI,J) +AD*U(NIM1,J))/(1. +AD)
!   END DO

    RETURN

  !---V momentum--------------------------------------------------------

    !*********
    ENTRY MODV
    !*********

    DO I = 1, NI
      V(I, 2)  = 0.0  !lower wall
      V(I, NJ) = 0.0  !upper wall
    END DO
!    DO I = 2, NIM1    ! symmtry plane at lower boundary
!      AS(I, 3) = 0.0
!    END DO
!
!   blank out block
!
    DO I = IB1,IB2-1
    DO J = JB1,JB2
       SU(I, J)  = 0.0
       AP(I, J)  = -GREAT
       RSV(I,J)  = 0.0
       V(I, J)   = 0.
    END DO
    END DO
!
    DO J=2,NJM1
      V(NI,J) = V(NIM1,J)
    ENDDO
!
! Convective outflow boundary condition
!
!   AD = DT*VELOC/(X(NI) -X(NIM1))
!   DO J=2, NJM1
!     V(NI,J)  = (VN(NI,J) +AD*V(NIM1,J))/(1. +AD)
!   END DO
    RETURN


    !-------------------------------------------------------------------

    !***********
    ENTRY MODVEL
    !***********

    AMAS = 0.0

    DO J = 2, NJM1
      AMAS = AMAS + DEN(2, J) * U(2, J) * R(J) * SNS(J)
    END DO

    SUM1 = 0.0
    SUM2 = 0.0
    DELU = 0.0

    DO J = 2, NJM1
      RHOA = DEN(NIM1, J) * R(J) * SNS(J)
      SUM1 = SUM1 + RHOA
      SUM2 = SUM2 + RHOA * U(NIM1, J)
    END DO

    DELU = (AMAS - SUM2) / SUM1

    DO J = 2, NJM1
      U(NIM1, J) = U(NIM1, J) + DELU
      U(NI, J)   = U(NIM1, J)
    END DO

   ! RETURN

  !---pressure correction-----------------------------------------------

    !*********
    ENTRY MODP
    !*********

!   DO J = 2, NJM1
!      SU(NIM1, J) = 0.0
!      SP(NIM1, J) = -GREAT
!   END DO
    DO J = JB1, JB2-1
    DO I = IB1, IB2-1
       SU(I, J) = 0.0
       SP(I, J) = -GREAT
    END DO
    END DO
    DO J=1,NJ
!     P(NIM1,J) = P(NI-2,J)
      P(NI,J)   = P(NIM1,J)
    ENDDO
!---alternative B.C. for P'-------------------------------------------
!the following boundary condition simulates the condition dp/dx =0
!   DO J=2,NJM1
!     AE(NIM1, J) = 0.0
!     AW(2, J)    = 0.0
!   ENDDO
!   DO I=2,NIM1
!     AN(I,NJM1) = 0.0
!     AS(I,2)    = 0.0
!   ENDDO

    RETURN

  !---thermal energy----------------------------------------------------

    !*********
    ENTRY MODT          ! no modifications for this problem
    !*********

    RETURN


  END SUBROUTINE PROMOD
