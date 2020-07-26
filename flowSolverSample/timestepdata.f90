SUBROUTINE TIMESTEPDATA(getLoopIterationNumber, getPrintEveryNumber)
 USE DATA
 IMPLICIT NONE
 INTEGER :: getLoopIterationNumber, getPrintEveryNumber
 character(len=70) :: filename
 integer :: outunit   !assign number to correspond to file



 outunit = getLoopIterationNumber + 5000

 IF (MOD(getLoopIterationNumber, getPrintEveryNumber) == 0) THEN
    WRITE(filename, '("timestep_",I7.7,".dat")') getLoopIterationNumber
    OPEN(unit = outunit, file = filename, form = 'formatted')

    WRITE(outunit, 654)
    WRITE(outunit, 655) NI, NJ

    ! Programmed by FL, implemented by WJ to this subroutine

    DO I = 1, NI
      VI(I, 1)  = 0.   ! zero v at lower wall
      VI(I, NJ) = 0.   ! zero v at pipe wall
    END DO
    DO J = 2, NJM1
      VI(1, J)  = 0.          ! zero v at inlet
      V(NI, J)  = V(NIM1,J)   ! v at exit extrapolated  from inside the pipe
    END DO
!
!   interpolate values for V first layer next to the J boundaries
!
    DO I = 2, NIM1
      VI(I, 2)   = VI(I,1)  +(V(I, 3) -VI(I,1))*(Y(2)-Y(1))/(YV(3) -Y(1))
      VI(I,NJM1) = VI(I,NJ) +(V(I,NJM1) -VI(I,NJ))  &
                            *(Y(NJ)-Y(NJM1)/(Y(NJ)-YV(NJM1)))
    END DO
!
!   interpolate interior values for V along j=const lines
!
    DO J = 3, NJ-2
      DO I = 2, NIM1
        VI(I, J) = (V(I, J) + V(I, J+1))/2.0
      END DO
    END DO
    DO J = 2, NJM1
      VI(NI, J) = VI(NIM1,J)   ! v at exit extrapolated  from inside the pipe
    END DO
!### Interpolate U-velocity
!### interpolate values for U at inlet and outlet
    DO J = 2, NJM1
      UI(1,J)    = U(1,J)
      UI(2,J)    = U(1,J)  +(U(3, J) -U(1,J))*(X(2)-X(1))/(XU(3) -X(1))
      UI(NIM1,J) = U(NIM1,J)  ! simple extrapolation at exit
      UI(NI,J)   = U(NIM1,J)  ! simple extrapolation at exit
      DO I = 3, NI-2
        UI(I, J) = .5*(U(I, J) +U(I+1, J))  !interpolate interior values for U
      END DO
    END DO
!###   Set boundary values for U at  j=1 and j= NJ
    DO I = 1, NI
      UI(I,1)  = 0.   ! zero u at lower wall
      UI(I,NJ) = 0.   ! zero u at upper wall
    END DO
!
! Extrapolate p
!
    DO J=1,NJ
!     P(1,J)  = P(2,J)
!     P(1,J)  = 1.0
      P(NI,J) = P(NIM1,J)
    ENDDO
    DO I=1,NI
      P(I,1)  = P(I,2)
      P(I,NJ) = P(I,NJM1)
    ENDDO

   
    write(outunit, 105) ((X(I),I=1,NI),J=1,NJ)
    write(outunit, 105) ((Y(J),I=1,NI),J=1,NJ)
    write(outunit, 105) ((U(I,J),I=1,NI),J=1,NJ)
    write(outunit, 105) ((V(I,J),I=1,NI),J=1,NJ)
    write(outunit, 105) ((P(I,J),I=1,NI),J=1,NJ)
    close(outunit)
    !  DO J = 2, NJ
    !   DO I = 1, NI
    !     WRITE (outunit, 1100)X(I)
    !   END DO
    ! END DO

    ! DO J = 2, NJ
    !   DO I=1, NI
    !     WRITE (outunit, 1100)Y(J)
    !   END DO
    ! END DO

    ! DO J = 2, NJ
    !   DO I = 1, NI
    !     WRITE (outunit, 1100) UI(I, J)
    !   END DO
    ! END DO
 
    ! DO J = 2, NJ
    !   DO I = 1, NI
    !     WRITE (outunit, 1100)VI(I, J)
    !   END DO
    ! END DO
 
    ! DO J = 2, NJ
    !   DO I = 1, NI
    !     WRITE (outunit, 1100) P(I, J)
    !   END DO
    ! END DO
    ! 
    ! CLOSE(outunit)
 
END IF

100 format('TITLE = "FLOW PAST BLOCK"')
102 format('VARIABLES = "X" "Y" "U" "V" "P"')
103 format('ZONE I=', I6, ', J=', I6, ', K=1', ', F= BLOCK')
105 format(10E14.5)
201 format('FILETYPE = SOLUTION')
202 format('VARIABLES = "U" "V" "P"')
204 format('T = "Time ', I4, '"')
    RETURN















 654 FORMAT('VARIABLES', 2X, '=', 2X,' "X", ', 2X, ' "Y", ', 2X, &
 ' "U", ', 2X, ' "V" ', 2X, ' "P" ')
 1100 FORMAT(1P13E10.3)
 655 FORMAT('ZONE I=', I6, ', J=', I6, ', K=1', ', F= BLOCK')



END SUBROUTINE TIMESTEPDATA

