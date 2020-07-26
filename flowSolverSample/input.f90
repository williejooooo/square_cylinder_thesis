SUBROUTINE INPUT
USE DATA
IMPLICIT NONE

!### Set defulat values

    VELOC    = 1.0    ! velocitiy at inlet
    DENSIT   = 1.0
    Reynolds = 100.
    PRANDT   = 1.0

    NI       = 302   ! number of grid lines in i and j-direction
    NJ       = 82
    IB1      = 102   !block locations: block wall in between nodal points designatied as IB
    IB2      = 112
    JB1      = 37
    JB2      = 47

    NTIME    = 400   ! Number of time steps to be computed
    MAXIT    = 800    ! maximum number of sub-iterations within each time step
    CFL      = 1.0
    SORMAX   = 1.0E-4    ! convergence criteria, max. residual source
    ALPHA    = 0.50
    NSWPU = 2   ! number of grid sweeps for u-velocity
    NSWPV = 1   ! number of grid sweeps for v-velocity
    NSWPP = 2   ! number of grid sweeps for pressure
    NSWPT = 1   ! number of grid sweeps for temperature

! ### Relaxation parameters
    URFU = 0.5     ! u-velocity
    URFV = 0.5     ! v-velocity
    URFP = 0.3     ! pressure
    URFT = 1.0     ! temperature

!   Original for steady pipe
!   URFU = 0.7     ! u-velocity
!   URFV = 0.7     ! v-velocity
!   URFP = 0.3     ! pressure
!   URFT = 1.0     ! temperature

    READ(5,*) 
    READ(5,*) NI,NJ,IB1,IB2,JB1,JB2
    READ(5,*) 
    READ(5,*) VELOC, DENSIT, XLEN, Reynolds, PRANDT, BFSCALE
    READ(5,*) 
    READ(5,*) NTIME, MAXIT, CFL, SORMAX, ALPHA, PDAMP0
    READ(5,*) 
    READ(5,*) NSWPU, NSWPV, NSWPP, NSWPT
    READ(5,*) 
    READ(5,*) URFU, URFV, URFP, URFT

    ALPHAC   = 1.0 -ALPHA
    PDAMP    = PDAMP0*ABS(ALPHAC) 

    INDCOS = 1 ! select coordinate system (1 = cartesian, 2 = cylindrical polar)
    W =  50.0  ! length of the channel
    H =   8.0  ! height from centerline to top (half height of the channel)

! equations to be solved ( .TRUE. -> solve equaation )
    INCALU = .TRUE.   ! u-velocity
    INCALV = .TRUE.   ! v-velocity
    INCALP = .TRUE.   ! pressure
    INCALT = .FALSE.  ! temperature
    INPRO  = .FALSE.   ! properties

    IF(Reynolds > 0) THEN
       VISCOS = DENSIT*VELOC*XLEN/Reynolds   ! properties of the fluid
       write(*,*) 'Read in Reynolds number Re=', Reynolds
    ELSE
       VISCOS = -Reynolds
       Reynolds = DENSIT*VELOC*XLEN/VISCOS
       write(*,*) 'Read in dynamic viscosity', VISCOS, ' Re=', Reynolds
    ENDIF

    TTOP  = 100.0    ! temperature at the boundaries
    TBOT  = TTOP
    TLEFT = 0.0
    TRITE = TLEFT
    TSTEP = 0.0

!### Reference Pressure point
    IPREF = 2       
    JPREF = 2

!### Flow Monitor point 
    IMON = 6      ! control variables for output
    JMON = 6

    RETURN
END SUBROUTINE INPUT
