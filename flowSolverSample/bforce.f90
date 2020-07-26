SUBROUTINE BFORCE
!### subroutine to calculate body force per unit volume
!###
USE DATA
IMPLICIT NONE
REAL FHZ,PI

!###
  FHZ = 0.5
  PI  = ATAN(1.0)*4.0
  DO I=IB1+5,IB2
  DO J=JB2,JB2+1
    BFX(I,J) = BFSCALE*sin(2.*PI*FHZ*time)
  ENDDO
  ENDDO
END SUBROUTINE BFORCE
