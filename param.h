! Dimensions of Physical Grid and array dimensions:
!----------- -- -------- ---- --- ----- -----------
! LLm   Number of the internal points of the PHYSICAL grid in XI-
! MMm   and ETA-directions, excluding physical side boundary points,
!       peroodic ghost points, and MPI-margins (if any).
!
! Lm    Number of the internal points [see above] of array covering
! Mm    a single MPI-subdomain.  These are identical to LLm, MMm if
!       there is no MPI-partitioning.

      integer, parameter ::
#if defined BASIN
     &               LLm=180, MMm=140, N=12
#elif defined CANYON_A
     &               LLm=65,  MMm=48,  N=10
#elif defined CANYON_B
     &               LLm=66,  MMm=48,  N=12
#elif defined DAMEE_B
     &               LLm=128, MMm=128, N=20
c**  &               LLm=256, MMm=256, N=20
#elif defined EKMAN_SPIRAL
     &               LLm=2, MMm=2, N=30

#elif defined GRAV_ADJ
     &               LLm=160, MMm=2,   N=40

#elif defined COLD_FILAMENT
c     &               LLm=400, MMm=2,   N=40
     &                LLm=800, MMm=2,   N=80
#elif defined NJ_BIGHT
c*   &               LLm=98,  MMm=206, N=12
     &               LLm=24,  MMm=34,  N=10
#elif defined PACIFIC

c    &               LLm=384, MMm=224, N=30
c    &               LLm=392, MMm=288, N=30

c    &               LLm=432, MMm=320, N=32
c     &               LLm=488, MMm=360, N=40  ! PAC44
     &                 LLm=976, MMm=720, N=40 ! PAC22

c     &               LLm=460, MMm=240, N=50  ! PACSMA grid

#elif defined PACIFIC_25KM
     &     LLm=920, MMm=480, N=32

#elif defined PACIFIC_2D
     &               LLm=768, MMm=512, N=1
c    &               LLm=1520, MMm=1088, N=1

#elif defined OVERFLOW
     &               LLm=4,   MMm=128, N=20
#elif defined SEAMOUNT
     &               LLm=41,  MMm=41,  N=30
c     &               LLm=80,  MMm=80,  N=30
c    &               LLm=192, MMm=96,  N=22
#elif defined ISWAKE
c     &               LLm=192, MMm=84, N=20
c     &               LLm=384, MMm=168, N=20
c     &               LLm=720, MMm=160, N=2
c     &              LLm=768, MMm=192, N=2
     &              LLm=1920, MMm=384, N=2
c     &               LLm=1440, MMm=320, N=2
c     &               LLm=768, MMm=336, N=20
c     &              LLm=1536, MMm=672, N=20
#elif defined SHELFRONT
     &               LLm=4,   MMm=40,  N=12
#elif defined SOLITON
c     &               LLm=96,  MMm=32,  N=1
c     &                 LLm=97, MMm=33, N=1
c     &               LLm=192, MMm=64,  N=1
     &               LLm=384, MMm=128, N=1
c     &                LLm=576, MMm=192, N=1
c     &               LLm=768, MMm=256, N=1
#elif defined BRAZIL
     &               LLm=168, MMm=272, N=40

#elif defined TASMAN_SEA
     &               LLm=128, MMm=128, N=4
#elif defined RIVER
     &               LLm=40,  MMm=160,  N=16
c     &                LLm=60,  MMm=240,  N=24
#elif defined UPWELLING
c     &                 LLm=20,  MMm=80,  N=16
     &                 LLm=20,  MMm=128,  N=40


c     &               LLm=20,  MMm=80,  N=32
#elif defined CANBAS2
c     &                LLm=224, MMm=288, N=32      ! CanBas
     &                LLm=225, MMm=328, N=32      ! NEA_EXT
c     &                LLm=384, MMm=480, N=32      ! GranCan

#elif defined USW_2D
     &               LLm=240, MMm=3, N=32
#elif defined USWC10KM_2D
     &               LLm=120, MMm=3, N=32
#elif defined USWC12
     &               LLm=320, MMm=448, N=42
#elif defined BALTIC
     &                LLm=440, MMm=384, N=32
#elif defined USWC4
     &               LLm=375, MMm=625, N=42
#elif defined USWC42
     &               LLm=435, MMm=660, N=60
#elif defined USSWC1
     &               LLm=200, MMm=200, N=60
#elif defined FLAT_BOTTOM
     &               LLm=10, MMm=10, N=10
#elif defined WindFarm
     &               LLm=250, MMm=350, N=42
#elif defined CHABU
     &                LLm=1200, MMm=1400, N=50
#elif defined USWEST
# ifdef GRID_LEVEL
#  if GRID_LEVEL == 1
     &               LLm=83,  MMm=168, N=20,         ! Monteray Bay,
     &               imin_child=40,  imax_child=71,  ! Level 1, 15 km
     &               jmin_child=54,  jmax_child=117
#  elif GRID_LEVEL == 2
     &               LLm=93,  MMm=189, N=20          ! Level 2, 5 km
#  endif
# else


c>>  &               LLm=72,   MMm=240, N=32    ! PEC2 of Xavier
c**  &               LLm=62,   MMm=126, N=40    ! SCB L0 grid
c**  &               LLm=83,   MMm=168, N=20    ! MB_L1
c**  &               LLm=126,  MMm=254, N=20    ! USWEST grid 16
!     &               LLm=312,  MMm=512, N=32    ! USW51 - lev0
     &               LLm=62,   MMm=126, N=20    ! USWEST grid 15 (20 km)


# endif
#elif defined WAVE_RAD
     &              LLm=384,  MMm=384, N=1
#else
     &                LLm=??, MMm=??, N=??
#endif


!
! MPI related variables
! === ======= =========
!
      integer Lmmpi,Mmmpi,iminmpi,imaxmpi,jminmpi,jmaxmpi
      common /comm_setup_mpi/ Lmmpi,Mmmpi,
     &                    iminmpi,imaxmpi,jminmpi,jmaxmpi


! Domain subdivision parameters:
!------- ----------- -----------
! NNODES             total number of MPI processes (nodes);
! NP_XI,  NP_ETA     number of MPI subdomains in XI-, ETA-directions;
! NSUB_X, NSUB_E     number of shared memory subdomains (tiles) in
!                                             XI- and ETA-directions;
      integer, parameter ::
#ifdef MPI
#  if defined FLAT_BOTTOM
     &      NP_XI=2, NP_ETA=2, NSUB_X=1, NSUB_E=1
#  else 
     &      NP_XI=3, NP_ETA=3, NSUB_X=1, NSUB_E=1
#  endif
#else
     &      NSUB_X=1, NSUB_E=1
#endif

! Array dimensions and bounds of the used portions of sub-arrays

#ifdef MPI
      integer, parameter :: NNODES=NP_XI*NP_ETA,
     &    Lm=(LLm+NP_XI-1)/NP_XI, Mm=(MMm+NP_ETA-1)/NP_ETA

      integer ocean_grid_comm, mynode,  iSW_corn, jSW_corn,
     &                         iwest, ieast, jsouth, jnorth
# ifndef EW_PERIODIC
      logical west_exchng,  east_exchng
# endif
# ifndef NS_PERIODIC
      logical south_exchng, north_exchng
# endif
      common /mpi_comm_vars/  ocean_grid_comm, mynode,
     &     iSW_corn, jSW_corn, iwest, ieast, jsouth, jnorth
# ifndef EW_PERIODIC
     &                , west_exchng,  east_exchng
# endif
# ifndef NS_PERIODIC
     &                , south_exchng, north_exchng
# endif
#else
      integer, parameter :: Lm=LLm, Mm=MMm
#endif

! Derived dimension parameters, number of tracers and tracer
! identification indices:
      integer, parameter :: padd_X=(Lm+2)/2-(Lm+1)/2,
     &                      padd_E=(Mm+2)/2-(Mm+1)/2
#ifdef SOLVE3D
     &       , itemp=1
# ifdef SALINITY
     &       , isalt=2
     &       , ntrc_salt=1
#  ifdef BIOLOGY
     &       , iNO3_=3, iNH4_=4, iDet_=5, iPhyt=6, iZoo_=7
#  endif
# else
             , ntrc_salt=0
#  ifdef BIOLOGY
     &       , iNO3_=2, iNH4_=3, iDet_=4, iPhyt=5, iZoo_=6
#  endif
# endif
# ifdef PASSIVE_TRACER
#  ifdef TRACER_BE7
      integer, parameter :: ntrc_pas=1
#  else /* TRACER_BE7 */
      integer, parameter :: ntrc_pas=6
#  endif /* TRACER_BE7 */
      integer itpas(ntrc_pas)
      common /pass_tracer/ itpas
# else
      integer, parameter :: ntrc_pas=0
# endif

# ifdef BIOLOGY_BEC
      integer, parameter :: itrc_bio=itemp+ntrc_salt+ntrc_pas+1
      integer, parameter :: iPO4=itrc_bio,iNO3=iPO4+1, iSIO3=iPO4+2, 
     &     iNH4=iPO4+3,
     &     iFE=iPO4+4, iO2=iPO4+5, iDIC=iPO4+6,
     &     iALK=iPO4+7, iDOC=iPO4+8, iSPC=iPO4+9,
     &     iSPCHL=iPO4+10, iSPCACO3=iPO4+11, iDIATC=iPO4+12,
     &     iDIATCHL=iPO4+13, iZOOC=iPO4+14, iSPFE=iPO4+15,
     &     iDIATSI=iPO4+16, iDIATFE=iPO4+17, iDIAZC=iPO4+18,
     &     iDIAZCHL=iPO4+19, iDIAZFE=iPO4+20, iDON=iPO4+21,
     &     iDOFE=iPO4+22, iDOP=iPO4+23
#  ifdef DMS_CYCLE
      integer, parameter :: iDMSP=iDOP+1,iDMS=iDMSP+1
#  endif /* DMS_CYCLE */
#  ifdef NCYCLE_ANOXIC
      integer, parameter :: iNO2 = iDOP + 1
#   ifdef DMS_CYCLE
     &     + 2
#   endif /* DMS_CYCLE */
      integer, parameter :: iN2O = iNO2 + 1, iN2 = iN2O + 1
#  endif /* NCYCLE_ANOXIC */
#  ifdef N15_CYCLE
      integer, parameter ::iNO3_15N = iDOP + 1
#   ifdef DMS_CYCLE
     &     + 2
#   endif /* DMS_CYCLE */
#   ifdef NCYCLE_ANOXIC
     &     + 3
#   endif /* NCYCLE_ANOXIC */
      integer, parameter :: iNH4_15N = iNO3_15N + 1
      integer, parameter :: iSp_15N = iNH4_15N + 1
      integer, parameter :: iDiat_15N = iSp_15N + 1
      integer, parameter :: iZoo_15N = iDiat_15N + 1
      integer, parameter :: iDiaz_15N = iZoo_15N + 1
      integer, parameter :: iDON_15N = iDiaz_15N + 1
#   ifdef NCYCLE_ANOXIC
      integer, parameter :: iNO2_15N = iDON_15N + 1
      integer, parameter :: iN2O_15N = iNO2_15N + 1
      integer, parameter :: iN2_15N = iN2O_15N + 1
#   endif /* NCYCLE_ANOXIC */
#  endif /* N15_CYCLE */

      integer, parameter :: ntrc_bio = 24
#  ifdef DMS_CYCLE
     &     + 2
#  endif /* DMS_CYCLE */
#  ifdef NCYCLE_ANOXIC
     &     + 3
#  endif /* NCYCLE_ANOXIC */
#  ifdef N15_CYCLE
     &     + 7
#   ifdef NCYCLE_ANOXIC
     &     + 3
#   endif /* NCYCLE_ANOXIC */
#  endif /* N15_CYCLE */
#  ifdef SEDIMENT_BIOLOGY
      integer, parameter :: NT_sed = 3      
#   ifdef N15_CYCLE
     &     + 1
#   endif 
#  endif /* SEDIMENT_BIOLOGY */
# else
      integer, parameter :: ntrc_bio = 0
# endif /* BIOLOGY_BEC */
      integer, parameter :: NT=itemp+ntrc_salt+ntrc_pas+ntrc_bio
#endif /*SOLVE3D */

#ifdef ANA_PSOURCE
#  if defined FLAT_BOTTOM
      integer, parameter :: Msrc=4  ! Number of point sources
#  endif
      integer, parameter :: Stepsrc=233   ! Number of time steps of the point sources
#endif
!#ifdef PSOURCE
!      integer, parameter :: Msrc=1   ! Number of point sources
!      integer, parameter :: Stepsrc=233   ! Number of time steps of the point sources
!#endif
#ifdef STATIONS
      integer, parameter :: NS=5      ! Number of output stations
#endif
#ifdef FLOATS
      integer, parameter :: Mfloats=32000 ! Maximum number of floats
#endif
