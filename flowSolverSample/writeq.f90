SUBROUTINE WRITEQ
USE DATA
IMPLICIT NONE


! Programed by FL 6/15/2020
!
! Interpolate the U, V solutions to the nodal points as UI, and VI and
! write the grid and solutions U,V,P at grid nodal points to files for Matlab and Tec360
! This version is for the block in a 2D channel.
!
!   Set boundary values for V
!
    DO I = 1, NI
      VI(I, 1)  = V(I,2)  
      VI(I, NJ) = V(I,NJ)
!### interpolate values for V first layer next to the J boundaries
       VI(I, 2)   = VI(I,1)  +(V(I, 3) -VI(I,1))*(Y(2)-Y(1))/(YV(3) -Y(1))
       VI(I,NJM1) = VI(I,NJ) +(V(I,NJM1) -VI(I,NJ))  &
                             *(Y(NJ)-Y(NJM1))/(Y(NJ)-YV(NJM1))
      DO J = 3, NJ-2
        VI(I, J) = 0.5*(V(I, J) +V(I, J+1))
      END DO
    END DO
!
!   interpolate interior values for V along j=const lines
!
!### Interpolate U-velocity
!### interpolate values for U at inlet and outlet
    DO J = 1, NJ
      UI(1,J)    = U(2,J)
      UI(NI,J)   = U(NIM1,J) 
      UI(2,J)    = UI(1,J)  +(U(3, J) -UI(1,J))*(X(2)-X(1))/(XU(3) -X(1))
      UI(NIM1,J) = UI(NI,J) +(U(NIM1,J) -UI(NI,J))  &
                           *(X(NI)-X(NIM1))/(X(NI)-XU(NIM1))
      DO I = 3, NI-2
        UI(I, J) = .5*(U(I, J) +U(I+1, J))  !interpolate interior values for U
      END DO
    END DO
!### write interpolated data on nodal points to tst.out
!   CALL PRINT (1, 1, NI, NJ, mx, my, X, Y, UI, 'U   ')
!   CALL PRINT (1, 1, NI, NJ, mx, my, X, Y, VI, 'V   ')
!   CALL PRINT (1, 1, NI, NJ, mx, my, X, Y, P, 'P   ')
!### write file for matlab use
!   OPEN(20, FILE='matlab.bin', FORM='UNFORMATTED')
!   WRITE(20) NI,NJ,(X(I),I=1,NI),(Y(J),J=1,NJ),((UI(I,J),I=1,NI),J=1,NJ), &
!             ((VI(I,J),I=1,NI),J=1,NJ), ((P(I,J),I=1,NI),J=1,NJ)
!   CLOSE(20)
!### WRITE TECPLOT FILE
!### WRITE TECPLOT GRID FILE
    IF(NWRITE==0) THEN
      write(21, *) 'TITLE = "FLOW PAST BLOCK"'
      write(21, 102)
      write(21, 103) NI, NJ
    ENDIF
      write(21, 204) NWRITE
      write(21, 105) ((X(I),I=1,NI),J=1,NJ)
      write(21, 105) ((Y(J),I=1,NI),J=1,NJ)
      write(21, 105) ((U(I,J),I=1,NI),J=1,NJ)
      write(21, 105) ((V(I,J),I=1,NI),J=1,NJ)
      write(21, 105) ((P(I,J),I=1,NI),J=1,NJ)
    
    flush(21)
    NWRITE = NWRITE +1
100 format('TITLE = "FLOW PAST BLOCK"')
102 format('VARIABLES = "X" "Y" "U" "V" "P"')
103 format('ZONE I=', I6, ', J=', I6, ', K=1', ', F= BLOCK')
105 format(10E14.5)
201 format('FILETYPE = SOLUTION')
202 format('VARIABLES = "U" "V" "P"')
204 format('T = "Time ', I4, '"')
    RETURN

!    IF(NWRITE==0) THEN
!      write(21, *) 'TITLE = "FLOW PAST BLOCK"'
!      write(21, *) 'FILETYPE = GRID'
!      write(21, 102)
!      write(21, 103) NI, NJ
!      write(21, 104)
!      write(21, 105) ((X(I),I=1,NI),J=1,NJ)
!      write(21, 105) ((Y(J),I=1,NI),J=1,NJ)
!    ELSE
!      write(21, 103) NI, NJ
!      write(21,*) 'VARSHARELIST = ([1-2]=1)'
!    ENDIF
!
!!### WRITE TECPLOT SOLUTION FILE
!    IF(NWRITE==0) THEN
!      write(22, 100)
!      write(22, 201)
!      write(22, 202)
!    ENDIF
!    write(22, 103) NI, NJ
!    write(22, 204)  NWRITE
!!     write(22, *) 'SOLUTIONTIME = ', float(NWRITE)
!    write(22, 104)
!    write(22, 105) ((U(I,J),I=1,NI),J=1,NJ)
!    write(22, 105) ((V(I,J),I=1,NI),J=1,NJ)
!    write(22, 105) ((P(I,J),I=1,NI),J=1,NJ)
!    flush(21)
!    flush(22)
!    NWRITE = NWRITE +1
!100 format('TITLE = "FLOW PAST BLOCK"')
!101 format('FILETYPE = GRID')
!102 format('VARIABLES = "X" "Y"')
!103 format('ZONE I=', I6, ', J=', I6, ', K=1')
!104 format('ZONETYPE = Ordered, DATAPACKING = BLOCK')
!105 format(10E14.5)
!201 format('FILETYPE = SOLUTION')
!202 format('VARIABLES = "U" "V" "P"')
!204 format('T = "Time ', I4, '"')
!    RETURN
END SUBROUTINE WRITEQ
