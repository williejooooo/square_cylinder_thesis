SUBROUTINE FORCEDATA
  USE DATA
  IMPLICIT NONE

  OPEN(23, FILE = 'Force_Data.csv', FORM = 'formatted')
  
  IF (NT == 1) THEN
    WRITE(23,300)
    WRITE(23,301)
  ENDIF
  
  WRITE(23,302) TIME, CL, CD, CLP, CDP, CLV, CDV 
  
  IF (NT == NTIME) THEN
    CLOSE(23)
  ENDIF


300 format('#Force Data')
301 format('#TIMESTEP, Cl, Cd, Clp, Cdp, Clv, Cdv')
302 format(10E14.5)

END SUBROUTINE FORCEDATA
