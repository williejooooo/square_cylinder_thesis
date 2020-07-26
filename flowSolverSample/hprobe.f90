SUBROUTINE hprobe(y_pos, print_interval)

     USE DATA
     IMPLICIT NONE

     INTEGER :: y_pos, print_interval
     CHARACTER(len = 70) :: filename
     INTEGER :: outunit

     outunit = 69

     IF (MOD(NT,print_interval) == 0) THEN
     
        WRITE(filename, '("hprobe_",I7.7,"_timestep_",I7.7,".csv")') y_pos,NT
        OPEN(unit = outunit, file = filename, form = 'formatted')

        WRITE(outunit, 6969)

        DO I = 1, NI
          WRITE (outunit, 1100) X(I), Y(y_pos), UI(I,y_pos), VI(I,y_pos), P(I, y_pos)
        END DO

        CLOSE(outunit)
     
     END IF
     
     6969 FORMAT(' #X ', 2X, ' Y ', 2X,' U ', 2X,  'V', 2X, 'P')
     1100 FORMAT(F8.4, 2x, F12.10, 2x, F12.10, 2x, F12.10, 2x, F12.10, 2x)

END SUBROUTINE hprobe
