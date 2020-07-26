PROGRAM Laminar_Flow 
! Elliptic PDE Solver
!-----------------------------------------------------------------------
! This program calculates laminar elliptic flows in two dimensional geometries.
! The geometry and initial values are given in the subroutine "DATA".
! All subroutines and all variables are contained as modules
! in the file "subroutines.f90".
! Used subroutines: - DATA 
!                   - INIT 
!                   - PROPS             
!                   - CALCU             
!                   - CALCV            
!                   - CALCP             
!                   - CALCT             
!                   - LISOLV          
!                   - PROMOD                           
! By default the equations for x-velocity (U), y-velocity (V) and
! pressure correction (p') are solved. 
! The equation for temperature (T) (or Enthalpy) is disabled;
! it can be enabled in subroutine "DATA".
! To compile the program type "f90 subroutines.f90 main.f90".
! Output files: - tst.out
!-----------------------------------------------------------------------

  USE DATA     ! all variables defined in "Subroutines"
  IMPLICIT NONE
  REAL :: SORCE0, SORCE, ALPHA0

!***********************************************************************
!---MAIN PROGRAM--------------------------------------------------------
!-----------------------------------------------------------------------
!***********************************************************************

! assign output files
    OPEN(UNIT = 7, FILE = 'tst.out', ACCESS = 'SEQUENTIAL', STATUS = 'UNKNOWN')
    OPEN(21,FILE='grid.dat', FORM='formatted')
    OPEN(22,FILE='solution.dat', FORM='formatted')
!### read data,  set up geometry and initialize flow field
    CALL INPUT
    CALL GEOM
    CALL INIT 
!****************     
                                                                        
! calculate properties 

!*****************                                                        
    CALL PROPS 
!*****************                                                  

    WRITE(6, 700)
    WRITE(6, 710) NTIME, MAXIT, Reynolds, CFL, DT 
    WRITE(7, 700)
    WRITE(7, 710) NTIME, MAXIT, Reynolds, CFL, DT 
700 FORMAT(3x,'NTIME', 4x, 'MAXIT', 4x, 'Reynolds', 6x, 'CFL',10x,'DT')
710 FORMAT(I6,3x,I6,2x,3G14.5) 

!### Loop time steps
    TIME = 0.0
!   CALL WRITEQ
    DO NT=1,NTIME
      IF (NT==1) THEN
        MAXIT = MAXIT*2.
        DT = DT*GREAT
        ALPHA0  = ALPHA
        ALPHA  = 1.0
      ELSE IF (NT==2) THEN
        MAXIT = MAXIT/2.
        DT = DT/GREAT
        ALPHA = ALPHA0
      ENDIF
      ALPHAC = 1. -ALPHA
      PDAMP = PDAMP0*ABS(ALPHAC)
      IF (NT>1) TIME = TIME +DT
      R_ONLY = .FALSE. 
      CALL BFORCE
      DO NITER=1,MAXIT ! inner iteration loop                                 
        CALL CALCU  
        CALL CALCV          
        CALL CALCP  
        IF (INCALT) CALL CALCT
!       CALL MODVEL
        IF(INPRO) CALL PROPS
        SORCE = AMAX1(RESORU, RESORV, RESORM, RESORT)
        IF(NT==1.AND.MOD(NITER,100)==1) WRITE(6, 311) NITER, RESORU, RESORV, RESORM,&
           RESORT, U(IMON, JMON), V(IMON, JMON), P(IMON, JMON),T(IMON, JMON)
        IF(NITER==1) SORCE0=SORCE
        IF(NITER == 60 .AND. SORCE > 1.0E6 * SORMAX) exit
        IF(SORCE < SORMAX) exit
      ENDDO    ! end inner iteration
      CALL FORCE
      CALL FORCEDATA
      WRITE(6,601) NT, NITER, RESORU, RESORV, RESORM, SORCE0, CL, CD, TIME
      WRITE(7,601) NT, NITER, RESORU, RESORV, RESORM, SORCE0, CL, CD, TIME
601   format(2I6,4E10.2,2x,2E13.5,2x,E15.6)
      
      ! Use block below for normal use
      IF(NT>(NTIME-400)) THEN 
       CALL WRITEQ
       CALL TIMESTEPDATA(NT,1)
       CALL MOD_TIMESTEPDATA(NT,1)
       CALL HPROBE(41,10)
       CALL HPROBE(42,10)
      ENDIF
      
      ! Use block below for debugging 
      !CALL WRITEQ
      !CALL TIMESTEPDATA(NT,1)
      !CALL MOD_TIMESTEPDATA(NT,1)
      !CALL HPROBE(43,10)
      !CALL HPROBE(42,10)
           
      flush(21)
      flush(22)
      UN = U
      VN = V
      PN = P
      R_ONLY = .TRUE. 
      CALL CALCU
      CALL CALCV
      RSUN = RSU
      RSVN = RSV
    ENDDO    ! end time step
!  CALL WRITEQ
   CLOSE(7)
   CLOSE(21)
   CLOSE(22)
   210 FORMAT(35X, 'Laminar Pipe Flow.   RE=',F7.2//)
   310 FORMAT(/1X, 'ITER ', 'I----------NORMALIZED RESIDUAL SOURCE SUMS---   & 
       ----------I I-------FIELD VALUES AT MONITORING LOCATION', '(', I2,    & 
       ', ', I2, ')', '--------I' / 2X, 'NO.' ,5X, 'UMOM', 6X, 'VMOM', 6X,   &
       'MASS', 6X,'ENER',32X,'U',9X,'V',9X,'P',9X,'T') 
   311 FORMAT(1H ,I5, 3X, 1P4E10.3, 5X, 1P4E10.3)       
END PROGRAM Laminar_Flow
