MODULE DATA
!-----------------------------------------------------------------------
!"PROGRAM Laminar Flow Elliptic PDE Solver".
!All variables are also defined in this module
!-----------------------------------------------------------------------

!---sepecification part-------------------------------------------------

  IMPLICIT NONE
  INTEGER, PARAMETER :: mx = 502, my = 82      ! maximum size of arrays

  INTEGER :: I, IP1, IM1, IM2, J, JP1, JM2, JM1, II, JJ, &
             N, NIM1, NJM1, NIM2, NJM2, &
             ISTA                         ! counter variables
  INTEGER :: IPREF, JPREF                 ! NI-1 , NJ-1
  INTEGER :: NI, NJ                       ! number of grid lines
  INTEGER :: IB1, JB1, IB2, JB2           ! block location
  INTEGER :: NITER                        ! iteration counter
  INTEGER :: MAXIT                        ! maximum number of iterations
  INTEGER :: NSWPU, NSWPV, NSWPT, NSWPP   ! number of grid sweeps
  INTEGER :: ITEST, JTEST                 ! printoutput variables

  INTEGER :: IMON, JMON                  ! printoutput variables
  INTEGER :: INDCOS                       ! index of coordinate system

  REAL :: W, H, DX, DY, VOL, RADIUS               ! geometry
  REAL :: VISCOS, DENSIT, PRANDT,  EXCHAT  ! properties
  REAL :: XLEN, REYNOLDS                   ! Length scale and Reynolds number
  REAL :: URFU, URFV, URFP, URFT           ! under relaxation factors
  REAL :: RESORU, RESORV, RESORT, RESORM, RESOR   ! residual sources
  REAL :: SMP, CP, CPO                     ! source term coefficients
  REAL :: AMFLIN, AMOMIN, TFLOIN, TIN      ! normalized residuals
  REAL :: SORMAX                           ! convergence criterium
  REAL :: SORVOL                           ! volume source
  REAL :: TTOP, TBOT, TLEFT, TRITE, TSTEP  ! boundary values temperature
  REAL :: VELOC                            ! initial values
  REAL :: PPREF                            ! reference pressure
  REAL :: GREAT                            ! large constant
  REAL :: AREAN, AREAS, AREAEW             ! areas of control volume
  REAL :: DENS, DENN, DENE, DENW           ! density
  REAL :: VISN, VISS, VISE, VISW           ! viscosity
  REAL :: GAMN, GAMS, GAME, GAMW           ! diffusitivity
  REAL :: APN, APS, APE, APW               ! finite difference coeff.
  REAL :: GN, GS, GE, GW, GNW, GSE, GSW    ! convection coefficients
  REAL :: CN, CS, CE, CW                   ! convection coefficients
  REAL :: DN, DS, DE, DW                   ! diffusion coefficients
  REAL :: PEN, PES, PE, PEW                ! Peclet numbers
  CHARACTER*4 :: HEDAN, HEDAS, HEDAE, HEDAW, HEDSP, HEDSU, &
          HEDTAU, HEDVOL, HEDU, HEDV, HEDT, HEDP      ! for headlines

  REAL, DIMENSION(1:mx) :: X, DXEP, DXPW, SEW, XU,  &
                           DXEPU, DXPWU, SEWU, DXP, &
                           DXM          ! grid variables x-diretion
  REAL, DIMENSION(1:mx) :: Y, DYNP, DYPS, SNS, YV, R, RV, &
                           DYNPV, DYPSV, DYP, DYM, SNSV,  &
                           RCV          ! grid variables y-direction
  REAL, DIMENSION(1:mx, 1:my) :: U, V, P, T, PP   ! field variables at t=t_(n+1)
  REAL, DIMENSION(1:mx, 1:my) :: UI, VI     ! interpolated velocities
  REAL, DIMENSION(1:mx, 1:my) :: AP, AN, AS, AW, AE ! finite diff coeff
  REAL, DIMENSION(1:mx, 1:my) :: SU, SP   ! linearized source terms
  REAL, DIMENSION(1:mx, 1:my) :: DU, DV   ! main coefficients
  REAL, DIMENSION(1:mx, 1:my) :: PEE      ! Peclet number
  REAL, DIMENSION(1:mx, 1:my) :: DEN, VIS, GAMH  ! properties
  REAL, DIMENSION(1:mx, 1:my) :: UN, VN, PN, RSUN, RSVN, RSU, RSV  ! field variables at t=tn
  REAL, DIMENSION(1:mx, 1:my) :: BFX, BFY   ! body force per unit volume
  REAL :: DT, RDT, TIME, TTIME, CFL, ALPHA, ALPHAC,PDAMP,PDAMP0
  REAL :: BFSCALE, CL, CD, CLP, CDP, CLV, CDV
  INTEGER :: NTIME, NT

  LOGICAL :: INCALU, INCALV, INCALP, INCALT, INPRO ! equations to be solved
  LOGICAL :: R_ONLY
  INTEGER :: NWRITE !N-th call of WRITEQ
end module data
