SUBROUTINE MOD_TIMESTEPDATA(getLoopIterationNumber, getPrintEveryNumber)
 USE DATA
 IMPLICIT NONE
 INTEGER :: getLoopIterationNumber, getPrintEveryNumber
 character(len=70) :: filename
 integer :: outunit   !assign number to correspond to file
 REAL, DIMENSION (1:my) :: U_def, V_def !saves data from first column to separate array
 REAL, DIMENSION (1:mx, 1:my) :: U_NEW, V_NEW !alter U and V to these new array

 outunit = getLoopIterationNumber + 10000

 IF (MOD(getLoopIterationNumber, getPrintEveryNumber) == 0) THEN
    
    DO J = 1, NJ
        U_def(J) = U(1,J)
        V_def(J) = V(1,J)
    END DO     
         
    
    DO J = 1, NJ
        DO I = 1, NI
           U_NEW(I,J) = U(I,J) - U_def(J)
           V_NEW(I,J) = V(I,J) - V_def(J)
        END DO
    END DO

    WRITE(filename, '("mod_timestep_",I7.7,".dat")') getLoopIterationNumber
    OPEN(unit = outunit, file = filename, form = 'formatted')

    WRITE(outunit, 654)
    WRITE(outunit, 655) NI, NJ
    write(outunit, 105) ((X(I),I=1,NI),J=1,NJ)
    write(outunit, 105) ((Y(J),I=1,NI),J=1,NJ)
    write(outunit, 105) ((U_NEW(I,J),I=1,NI),J=1,NJ)
    write(outunit, 105) ((V_NEW(I,J),I=1,NI),J=1,NJ)
    write(outunit, 105) ((P(I,J),I=1,NI),J=1,NJ)
    close(outunit)
 
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



END SUBROUTINE MOD_TIMESTEPDATA
