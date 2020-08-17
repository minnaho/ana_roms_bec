!
!
#include "cppdefs.h"
#ifdef BIOLOGY_BEC

       subroutine init_scalars_bec(ierr)
!
! Set initial values for  globally accessible (stored in common
! blocks) scalar variables of the BEC model.
!
      implicit none
      integer ierr, i,j
      integer omp_get_num_threads
#include "param.h"
#include "scalars.h"
#include "ncvars.h"
#include "comm_vars.h"
#include "ncroms_bgcflux.h"
#if (defined DAILYPAR_PHOTOINHIBITION || defined DAILYPAR_BEC) && !defined DIURNAL_SRFLUX
#include "forces.h"
#endif
!
! Find out how many threads are created and check whether the number
! of threads exceeds maximum allowed, as well as check if the number
! of tiles is divisible by the number of threads, so the job can be
! evenly partitioned, complain and signal to terminate, if something
! is wrong. 
! 
C$OMP PARALLEL
C$OMP CRITICAL (isca_cr_rgn)
      numthreads=1
C$      numthreads=omp_get_num_threads()
C$OMP END CRITICAL (isca_cr_rgn)
C$OMP END PARALLEL
      mpi_master_only write(*,'(1x,A,3(1x,A,I3),A)') 'NUMBER',
     &    'OF THREADS:',numthreads,'BLOCKING:',NSUB_X,'x',NSUB_E,'.'

#ifdef MPI
# define NUMPROC NNODES
#else
# define NUMPROC 1
#endif

      if (numthreads.gt.NUMPROC) then
        mpi_master_only write(*,'(/1x,A,I3/)')
     &    'ERROR: Requested number of threads exceeds setting: ',NUMPROC
        ierr=ierr+1                                      !--> ERROR
      elseif (mod(NSUB_X*NSUB_E,numthreads).ne.0) then
        mpi_master_only write(*,
     &                '(/1x,A,1x,A,I3,4x,A,I3,4x,A,I4,A)') 'ERROR:',
     &                'wrong choice of numthreads =', numthreads,
     &                'NSUB_X =', NSUB_X, 'NSUB_E =', NSUB_E, '.'
        ierr=ierr+1                                      !--> ERROR
      endif
!
      do i=0,2*NT+1
        global_sum(i)=QuadZero
      enddo

!
! Define names of variables in NetCDF output files.
! Inner dimension is for variable type and outer is for variable
! attributes (name, long-name, units, field).
!
#ifdef SOLVE3D
      vname(1,indxPO4)='PO4                                       '
      vname(2,indxpO4)='Phosphate                                 '
      vname(3,indxpO4)='mMol P m-3                                '

CR      write(*,'(1x,A,I3)') 'indxNO3 =',indxNO3
      vname(1,indxNO3)='NO3                                       '
      vname(2,indxNO3)='Nitrate                                   '
      vname(3,indxNO3)='mMol N m-3                                '

      vname(1,indxSIO3)='SiO3                                      '
      vname(2,indxSIO3)='Silicate                                  '
      vname(3,indxSIO3)='mMol Si m-3                               '

CR      write(*,'(1x,A,I3)') 'indxNH4 =',indxNH4
      vname(1,indxNH4)='NH4                                       '
      vname(2,indxNH4)='Ammonium                                  '
      vname(3,indxNH4)='mMol N m-3                                '


      vname(1,indxFE)='Fe                                        '
      vname(2,indxFE)='Iron                                      '
      vname(3,indxFE)='mMol Fe m-3                               '

      vname(1,indxO2)='O2                                        '
      vname(2,indxO2)='Oxygen                                    '
      vname(3,indxO2)='mMol O2 m-3                               '

      vname(1,indxDIC)='DIC                                       '
      vname(2,indxDIC)='Dissolved inorganic carbon                '
      vname(3,indxDIC)='mMol C m-3                                '

      vname(1,indxALK)='Alk                                       '
      vname(2,indxALK)='Alkalinity                                '
      vname(3,indxALK)='mMol m-3                                '

      vname(1,indxDOC)='DOC                                       '
      vname(2,indxDOC)='Dissolved organic carbon                  '
      vname(3,indxDOC)='mMol C m-3                                '

      vname(1,indxSPC)='SPC                                       '
      vname(2,indxSPC)='Small phytoplankton carbon                '
      vname(3,indxSPC)='mMol C m-3                                '

      vname(1,indxSPCHL)='SPCHL                                     '
      vname(2,indxSPCHL)='Small phytoplankton chlorophyll           '
      vname(3,indxSPCHL)='mg Chl-a m-3                              '

      vname(1,indxSPCACO3)='SPCACO3                                   '
      vname(2,indxSPCACO3)='Small phytoplankton CaCO3                 '
      vname(3,indxSPCACO3)='mMol CaCO3 m-3                            '

      vname(1,indxDIATC)='DIATC                                     '
      vname(2,indxDIATC)='Diatom carbon                             '
      vname(3,indxDIATC)='mMol C m-3                                '

      vname(1,indxDIATCHL)='DIATCHL                                   '
      vname(2,indxDIATCHL)='Diatom chlorophyll                        '
      vname(3,indxDIATCHL)='mg Chl-a m-3                              '

      vname(1,indxZOOC)='ZOOC                                     '
      vname(2,indxZOOC)='Zooplankton                              '
      vname(3,indxZOOC)='mMol C m-3                               '

      vname(1,indxSPFE)='SPFE                                     '
      vname(2,indxSPFE)='Small phytoplankton iron                 '
      vname(3,indxSPFE)='mMol Fe m-3                              '

      vname(1,indxDIATSI)='DIATSI                                    '
      vname(2,indxDIATSI)='Diatom silicon                            '
      vname(3,indxDIATSI)='mMol Si m-3                               '

      vname(1,indxDIATFE)='DIATFE                                    '
      vname(2,indxDIATFE)='Diatom Iron                               '
      vname(3,indxDIATFE)='mMol Fe m-3                               '

#ifdef PHAEOCYSTIS
      vname(1,indxDiazc)='PHAEC                                    '
      vname(2,indxDiazc)='Phaeocystis carbon                       '
      vname(3,indxDiazc)='mMol C m-3                               '

      vname(1,indxDiazchl)='PHAECHL                                  '
      vname(2,indxDiazchl)='Phaeocystis chlorophyll                  '
      vname(3,indxDiazchl)='mg Chl-a m-3                             '

      vname(1,indxDiazfe)='PHAEFE                                   '
      vname(2,indxDiazfe)='Phaeocystis iron                         '
      vname(3,indxDiazfe)='mMol Fe m-3                              '
#else /* PHAEOCYSTIS */
      vname(1,indxDiazc)='DIAZC                                    '
      vname(2,indxDiazc)='Diazotroph carbon                        '
      vname(3,indxDiazc)='mMol C m-3                               '

      vname(1,indxDiazchl)='DIAZCHL                                  '
      vname(2,indxDiazchl)='Diazotroph chlorophyll                   '
      vname(3,indxDiazchl)='mg Chl-a m-3                             '

      vname(1,indxDiazfe)='DIAZFE                                   '
      vname(2,indxDiazfe)='Diazotroph iron                          '
      vname(3,indxDiazfe)='mMol Fe m-3                              '
#endif /* PHAEOCYSTIS */

      vname(1,indxDon)='DON                                      '
      vname(2,indxDon)='Dissolved organic nitrogen               '
      vname(3,indxDon)='mMol N m-3                               '

      vname(1,indxDofe)='DOFE                                     '
      vname(2,indxDofe)='Dissolved organic iron                   '
      vname(3,indxDofe)='mMol Fe m-3                              '

      vname(1,indxDop)='DOP                                      '
      vname(2,indxDop)='Dissolved organic phosphorus             '
      vname(3,indxDop)='mMol P m-3                               '

# ifdef DMS_CYCLE
      vname(1,indxDMSP)='DMSP                                     '
      vname(2,indxDMSP)='Dimethylsulfide particulate              '
      vname(3,indxDMSP)='mMol S m-3                               '

      vname(1,indxDMS)='DMS                                      '
      vname(2,indxDMS)='Dimethylsulfide                          '
      vname(3,indxDMS)='mMol S m-3                               '
# endif /* DMS_CYCLE */

# ifdef NCYCLE_ANOXIC
      vname(1,indxNO2)='NO2                                      '
      vname(2,indxNO2)='Nitrite                                  '
      vname(3,indxNO2)='mMol N m-3                               '

      vname(1,indxN2O)='N2O                                      '
      vname(2,indxN2O)='Nitrous oxide                            '
      vname(3,indxN2O)='mMol N2 m-3                              '

      vname(1,indxN2)='N2                                      '
      vname(2,indxN2)='N2 excess                                  '
      vname(3,indxN2)='mMol N2 m-3                              '
# endif /* NCYCLE_ANOXIC */

# ifdef N15_CYCLE
      vname(1,indxNO3_15N)='NO3_15N                                 '
      vname(2,indxNO3_15N)='Nitrate-15N                           '
      vname(3,indxNO3_15N)='mMol 15N m-3                          '

      vname(1,indxNH4_15N)='NH4_15N                                '
      vname(2,indxNH4_15N)='Ammonium-15N                            '
      vname(3,indxNH4_15N)='mMol 15N m-3                            '

      vname(1,indxSp_15n)='Sp_15N                                   '
      vname(2,indxSp_15n)='Small phytoplankton-15N                  '
      vname(3,indxSp_15n)='mMol 15N m-3                             '

      vname(1,indxDiat_15n)='Diat_15N                                 '
      vname(2,indxDiat_15n)='Diatom-15N                               '
      vname(3,indxDiat_15n)=' mMol 15N m-3                            '

      vname(1,indxZoo_15n)='Zoo_15N                                 '
      vname(2,indxZoo_15n)='Zooplankton-15N                         '
      vname(3,indxZoo_15n)='mMol 15N m-3                            '

#   ifdef PHAEOCYSTIS
      vname(1,indxDiaz_15n)='Phae_15N                                 '
      vname(2,indxDiaz_15n)='Phaeocystis-15N                           '
#   else
      vname(1,indxDiaz_15n)='Diaz_15N                                 '
      vname(2,indxDiaz_15n)='Diazotroph-15N                           '
#   endif
      vname(3,indxDiaz_15n)='mMol 15N m-3                             '

      vname(1,indxDON_15N)='DON_15N                                   '
      vname(2,indxDON_15N)='Dissolved organic nitrogen-15N           '
      vname(3,indxDON_15N)='mMol 15N m-3                             '

#  ifdef NCYCLE_ANOXIC
      vname(1,indxNO2_15N)='NO2_15N                                  '
      vname(2,indxNO2_15N)='Nitrite-15N                              '
      vname(3,indxNO2_15N)='mMol 15N m-3                             '

      vname(1,indxN2O_15N)='N2O_15N                                  '
      vname(2,indxN2O_15N)='Nitrous oxide-15N                        '
      vname(3,indxN2O_15N)='mMol 15N m-3                             '

      vname(1,indxN2_15N)='N2_15N                                  '
      vname(2,indxN2_15N)='N2-15N excess                           '
      vname(3,indxN2_15N)='mMol 15N m-3                             '
#  endif /* NCYCLE_ANOXIC */      
# endif /* N15_CYCLE */

#ifdef SEDIMENT_BIOLOGY
      vname(1,indxSedOrgC)='Sed_POC'
      vname(2,indxSedOrgC)='Sedimented Particulate Organic Carbon'
      vname(3,indxSedOrgC)='mMol C m-2'

      vname(1,indxSedCaCO3)='Sed_CaCO3'
      vname(2,indxSedCaCO3)='Sedimented Particulate Calcium Carbonate'
      vname(3,indxSedCaCO3)='mMol C m-2'

      vname(1,indxSedSi)='Sed_Si'
      vname(2,indxSedSi)='Sedimented Particulate Silicate'
      vname(3,indxSedSi)='mMol Si m-2'

# ifdef N15_CYCLE
      vname(1,indxSedOrgN15)='Sed_PON15'
      vname(2,indxSedOrgN15)='Sedimented Particulate Nitrogen-15'
      vname(3,indxSedOrgN15)='mMol 15N m-2'
# endif
#endif /* SEDIMENT_BIOLOGY */

CDL: pH is treated here for rst files;
CDL: and in init_scalars_becflux for flux files:
      vname(1,indxPH_rst)='PH'
      vname(2,indxPH_rst)='pH value'
      vname(3,indxPH_rst)=' '
!DL: some more variables which are to be output also to
!DL: averages/history/restart files:
      vname(1,indxPCO2_rst)='pCO2                                 '
      vname(2,indxPCO2_rst)='Surface water pCO2                   '
      vname(3,indxPCO2_rst)='ppm                                  '

      vname(1,indxPCO2air_rst)='pCO2air                           '
      vname(2,indxPCO2air_rst)='Atmospheric pCO2                  '
      vname(3,indxPCO2air_rst)='ppm                               '

      vname(1,indxPARinc_rst)='PARinc                                '
      vname(2,indxPARinc_rst)='Incoming photosynth. active radiation '
      vname(3,indxPARinc_rst)='W m-2                                 '

      vname(1,indxPAR_rst)='PAR                                   '
      vname(2,indxPAR_rst)='Photosynthetically active radiation   '
      vname(3,indxPAR_rst)='W m-2                                 '
# ifdef OUTPUT_SPECDOM
      vname(1,indxFreqDomSP_sfc)='FreqDomSp_sfc                        '
      vname(2,indxFreqDomSP_sfc)='Frequency of species dom. (sp, sfc)  '
      vname(3,indxFreqDomSP_sfc)='-                                    '
      vname(1,indxFreqDomDIAT_sfc)='FreqDomDiat_sfc                    '
      vname(2,indxFreqDomDIAT_sfc)=
     &     'Frequency of species dom. (diat, sfc)'
      vname(3,indxFreqDomDIAT_sfc)='-                                 '
#  ifdef PHAEOCYSTIS
      vname(1,indxFreqDomDIAZ_sfc)='FreqDomPhae_sfc                    '
      vname(2,indxFreqDomDIAZ_sfc)=
     &     'Frequency of species dom. (phae, sfc)'
#  else
      vname(1,indxFreqDomDIAZ_sfc)='FreqDomDiaz_sfc                    '
      vname(2,indxFreqDomDIAZ_sfc)=
     &     'Frequency of species dom. (diaz, sfc)'
#  endif
      vname(3,indxFreqDomDIAZ_sfc)='-                                  '
!
      vname(1,indxFreqDomSP_int)='FreqDomSp_int                        '
      vname(2,indxFreqDomSP_int)='Frequency of species dom. (sp, int)  '
      vname(3,indxFreqDomSP_int)='-                                    '
      vname(1,indxFreqDomDIAT_int)='FreqDomDiat_int                    '
      vname(2,indxFreqDomDIAT_int)=
     &     'Frequency of species dom. (diat, int)'
      vname(3,indxFreqDomDIAT_int)='-                                  '
#  ifdef PHAEOCYSTIS
      vname(1,indxFreqDomDIAZ_int)='FreqDomPhae_int                    '
      vname(2,indxFreqDomDIAZ_int)=
     &     'Frequency of species dom. (phae, int)'
#  else
      vname(1,indxFreqDomDIAZ_int)='FreqDomDiaz_int                    '
      vname(2,indxFreqDomDIAZ_int)=
     &     'Frequency of species dom. (diaz, int)'
#  endif
      vname(3,indxFreqDomDIAZ_int)='-                                  '
      if (indxFreqDomDIAZ_int .gt. 44+NT) stop 'vname too small'
# endif /* OUTPUT_SPECDOM */
#endif /* SOLVE3D */
      return
      end

      subroutine init_scalars_bec_part2(ierr)
      implicit none
      integer ierr
#if (defined DAILYPAR_PHOTOINHIBITION || defined DAILYPAR_BEC) && !defined DIURNAL_SRFLUX
#include "param.h"
#include "scalars.h"      
#include "forces.h"      
      iptr_srflx_day = 0
      iptr_srflx_day_set = 0
      PARinc_rst_read = 0
      num_srflx_day = n_srflx_day * 3600.0/dt  ! dt is in seconds
      if (num_srflx_day .gt. MAX_NUM_SRFLX_DAY) then
         print *, 'increase MAX_NUM_SRFLX_DAY to ', num_srflx_day
         ierr = ierr + 1
      end if
      srflx_frac = 1.0/num_srflx_day
      srflx_day = 0.0
#endif /* (DAILYPAR_PHOTOINHIBITION ||DAILYPAR_BEC) && ! DIURNAL_SRFLUX */
      end subroutine init_scalars_bec_part2
#else /* BIOLOGY_BEC */
      subroutine init_scalars_bec_empty()
      return
      end
#endif /* BIOLOGY_BEC */
