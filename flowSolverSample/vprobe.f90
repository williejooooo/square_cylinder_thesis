SUBROUTINE vprobe(y_position)

     USE DATA
     IMPLICIT NONE

     INTEGER :: Y_Position
     CHARACTER(len = 70) :: filename
     INTEGER :: outunit

     outunit = 69

     WRITE(filename, '("probe_at_horiz_",I7.7,".csv")') Y_Position
     OPEN(unit = outunit, file = filename, form = 'formatted')

     WRITE(outunit, 6969)

     DO I = 2, NI
       WRITE (outunit, 1100) X(I), Y(Y_Position), UI(I,Y_Position), VI(I,Y_Position)
     END DO

     CLOSE(outunit)


     6969 FORMAT(' #X ', 2X, ' Y ', 2X,' U ', 2X,  'V')
     1100 FORMAT(F12.10, 2x, F12.10, 2x, F12.10, 2x, F12.10, 2x)





END SUBROUTINE
