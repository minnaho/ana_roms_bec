! This is "ncvars.h":   indices for character array "vname" to keep
!----- -- -----------   names of netCDF variables and their attributes
! indxTime      time
! indxZ         free-surface
! indxUb,Vb     vertically integrated 2D U,V-momentum components
! indxU,V       3D u- and v-velocity components
! indxT,S,...,Zoo tracers (temperature, salinity, biotracers).
! indxO         "omega" vertical velocity (flux)
! indxW         "true" vertical velocity
! indxR         density anomaly
! indxAkv,Akt,Aks vertical viscosity/diffusivity coefficients
! indxHbls      depth of planetary boundary layer in KPP model
! indxHbbl      depth of bottom boundary layer in BKPP model

! indxAi        fraction of cell covered by ice
! indxUi,Vi     U,V-components of sea ice velocity
! indxHi,HS     depth of ice cover and depth of snow cover
! indxTIsrf     temperature of ice surface

! indxBSD,BSS bottom sediment grain Density and Size.
! indxWWA,WWD,WWP  wind induced wave Amplitude, Direction,and Period

      integer, parameter :: indxTime=1, indxZ=2, indxUb=3, indxVb=4
#ifdef SOLVE3D
     &                    , indxU=5, indxV=6, indxO=7, indxW=8
     &                    , indxR=9, indxT=10
# ifdef SALINITY
     &                    , indxS=indxT+1
# endif
# ifdef BIOLOGY
#  ifdef SALINITY
     &                    , indxNO3=indxS+1
#  else
     &                    , indxNO3=indxT+1
#  endif
     &                    , indxNH4 =indxNO3+1, indxChla=indxNO3+2
     &                    , indxPhyt=indxNO3+3, indxZoo =indxNO3+4
     &                    , indxSDet=indxNO3+5, indxLDet=indxNO3+6
# endif
# ifdef BIOLOGY_BEC
       integer, parameter :: indxPO4=indxT+ntrc_salt+ntrc_pas+1,
     &           indxNO3 =indxPO4+1, indxSio3=indxPO4+2,
     &           indxNH4 =indxPO4+3, indxFe=indxPO4+4,
     &           indxO2 =indxPO4+5, indxDic=indxPO4+6,
     &           indxAlk =indxPO4+7, indxDoc=indxPO4+8,
     &           indxSpc =indxPO4+9, indxSpchl=indxPO4+10,
     &           indxSpcaco3 =indxPO4+11, indxDiatc=indxPO4+12,
     &           indxDiatchl =indxPO4+13, indxZooc=indxPO4+14,
     &           indxSpfe =indxPO4+15, indxDiatsi=indxPO4+16,
     &           indxDiatfe =indxPO4+17, indxDiazc=indxPO4+18,
     &           indxDiazchl =indxPO4+19, indxDiazfe=indxPO4+20,
     &           indxDon =indxPO4+21, indxDofe=indxPO4+22,
     &           indxDop =indxPO4+23
#  ifdef DMS_CYCLE
       integer, parameter :: indxDMSP = indxDOP+1, indxDMS = indxDMSP+1
#  endif /* DMS_CYCLE */
#  ifdef NCYCLE_ANOXIC
       integer, parameter :: indxNO2 = indxDOP + 1
#   ifdef DMS_CYCLE
     &      + 2
#   endif /* DMS_CYCLE */
      integer, parameter :: indxN2O = indxNO2 + 1 
      integer, parameter :: indxN2 = indxN2O + 1 
#  endif /* NCYCLE_ANOXIC */
#  ifdef N15_CYCLE
       integer, parameter :: indxNO3_15N = indxDOP + 1
#   ifdef DMS_CYCLE
     &      + 2
#   endif /* DMS_CYCLE */
#   ifdef NCYCLE_ANOXIC
     &      + 3
#   endif /* NCYCLE_ANOXIC */
       integer, parameter :: indxNH4_15N = indxNO3_15N + 1
       integer, parameter :: indxSp_15n = indxNH4_15N + 1
       integer, parameter :: indxDiat_15n = indxSp_15n + 1
       integer, parameter :: indxZoo_15n = indxDiat_15n + 1 
       integer, parameter :: indxDiaz_15n = indxZoo_15n + 1
       integer, parameter :: indxDON_15N = indxDiaz_15n + 1
#   ifdef NCYCLE_ANOXIC
       integer, parameter :: indxNO2_15N = indxDON_15N + 1
       integer, parameter :: indxN2O_15N = indxNO2_15N + 1
       integer, parameter :: indxN2_15N = indxN2O_15N + 1
#   endif /* NCYCLE_ANOXIC */
#  endif /* N15_CYCLE */

       integer, parameter :: indxPH_rst = indxDOP+1
#  ifdef DMS_CYCLE
     &      + 2
#  endif /* DMS_CYCLE */
#  ifdef NCYCLE_ANOXIC
     &      + 3
#  endif /* NCYCLE_ANOXIC */
#  ifdef N15_CYCLE
     &      + 7
#   ifdef NCYCLE_ANOXIC
     &      + 3
#   endif /* NCYCLE_ANOXIC */
#  endif /* N15_CYCLE */

       integer, parameter :: indxPCO2_rst = indxPH_rst+1,
     &      indxPCO2air_rst = indxPCO2_rst+1,
     &      indxPARinc_rst = indxPCO2air_rst+1,
     &      indxPAR_rst = indxPARinc_rst+1
#  ifdef SEDIMENT_BIOLOGY
      integer, parameter :: indxSedOrgC = indxPAR_rst + 1
      integer, parameter :: indxSedCaCO3 = indxSedOrgC + 1
      integer, parameter :: indxSedSi = indxSedCaCO3 + 1
#   ifdef N15_CYCLE
      integer, parameter :: indxSedOrgN15 = indxSedSi + 1
#   endif
#  endif /* SEDIMENT_BIOLOGY */
#  ifdef OUTPUT_SPECDOM
      integer, parameter :: indxFreqDomSP_sfc = indxPAR_rst + 1
#   ifdef SEDIMENT_BIOLOGY
     &     + NT_sed
#   endif
      integer, parameter :: indxFreqDomDIAT_sfc = indxFreqDomSP_sfc + 1
      integer, parameter :: indxFreqDomDIAZ_sfc = indxFreqDomDIAT_sfc+1
      integer, parameter :: indxFreqDomSP_int = indxFreqDomDIAZ_sfc + 1
      integer, parameter :: indxFreqDomDIAT_int = indxFreqDomSP_int + 1
      integer, parameter :: indxFreqDomDIAZ_int = indxFreqDomDIAT_int+1
#  endif /* OUTPUT_SPECDOM */
# endif /* BIOLOGY_BEC */

# ifdef SEDIMENT_BIOLOGY
#  ifdef BIOLOGY_NPZDOC
       integer, parameter :: indxSedFirst = indxSedOrgN
#  elif defined BIOLOGY_BEC
       integer, parameter :: indxSedFirst = indxSedOrgC
#  endif
# endif /* SEDIMENT_BIOLOGY */

       integer, parameter :: indxAkv=indxT+NT
# ifdef BIOLOGY_BEC
     &      + 5 ! for indxPH, indxPCO2_rst, indxPCO2air_rst, indxPARinc_rst, indxPAR_rst
#  ifdef SEDIMENT_BIOLOGY
     &     + NT_sed
#  endif 
#  ifdef OUTPUT_SPECDOM
     &      + 6
#  endif
# endif /* BIOLOGY_BEC */
       integer, parameter :: indxAkt=indxAkv+1
# ifdef SALINITY
     &                    , indxAks=indxAkt+1
# endif
# ifdef LMD_KPP
#  ifdef SALINITY
     &                    , indxHbls=indxAks+1
#  else
     &                    , indxHbls=indxAkt+1
#  endif
# endif
# ifdef LMD_BKPP
     &                    , indxHbbl=indxHbls+1
# endif

# ifdef DIAGNOSTICS_EDDY
      integer indxeddyuu,indxeddyvv,indxeddyuv,indxeddyub,
     &        indxeddyvb,indxeddywb,indxeddynotused1,indxeddynotused2,
     &        indxeddynotused3
#  if defined LMD_BKPP
      parameter (indxeddyuu=indxHbbl+1,
#  else
      parameter (indxeddyuu=indxHbls+1
#  endif
     &           indxeddyvv=indxeddyuu+1,
     &           indxeddyuv=indxeddyvv+1,
     &           indxeddyub=indxeddyuv+1,
     &           indxeddyvb=indxeddyub+1,
     &           indxeddywb=indxeddyvb+1,
     &           indxeddynotused1=indxeddywb+1,
     &           indxeddynotused2=indxeddynotused1+1,
     &           indxeddynotused3=indxeddynotused2+1)
# endif



# ifdef SG_BBL96
#  ifndef ANA_WWAVE
     &                    , indxWWA=???          indxWWD=indxWWA+1
     &                    , indxWWP=indxWWA+2
#  endif
# endif
#endif
#ifdef ICE
      integer, parameter :: indxAi=????,     indxUi=indxAi+1,
     &                    , indxVi=indxAi+2, indxHi=indxAi+3,
     &                      indxHS=indxAi+4, indxTIsrf=indxAi+5
#endif

! Length of netCDF variable "time_step"

      integer, parameter :: iaux=6


! Naming conventions for indices, variable IDs, etc...
!------- ----------- --- -------- -------- ---- ------
! prefix ncid_  means netCDF ID for netCDF file
!        nrec_  record number in netCDF file since initialization
!        nrpf_  maximum number of records per file  (output netCDF
!                                                       files only)
! prefix/ending rst_/_rst refers to restart  netCDF file
!               his_/_his           history
!               avg_/_avg           averages
!               stn_/_stn           stations
!                    _frc           forcing
!                    _clm           climatology
!
! endings refer to:  ___Time  time [in seconds]
!                    ___Tstep time step numbers and record numbers
! all objects with   ___Z     free-surface
! these endings are  ___Ub    vertically integrated 2D U-momentum
! either:            ___Vb    vertically integrated 2D V-momentum
!
!  netCDF IDs, if    ___U     3D U-momentum
!  occur with prefix ___V     3D V-momentum
!  rst/his/avg/stn   ___T(NT) tracers
!                    ___R     density anomaly
! or                 ___O     omega vertical velocity
!                    ___W     true vertical velocity
!  parameter indices 
!  if combined with  ___Akv   vertical viscosity coefficient
!  prefix "indx"     ___Akt   vertical T-diffusion coefficient
!  (see above).      ___Aks   vertical S-diffusion coefficient
!                    ___Hbl   depth of KPP surface boundary layer
!                    ___Hbbl  thickness of KPP bottom boundary layer
!
! vname    character array for variable names and attributes;

#ifdef BULK_FLUX

      integer, parameter :: max_blk_file=8

      integer max_blk, ncidbulk(max_blk_file), nrst, ncidrst, nrecrst,
     &      nrrec, nrpfrst, ncidclm,nwrt, ncidhis, nrechis,
     &      nrpfhis
# ifdef BIOLOGY_BEC
     &     , ntdust, ntiron , ntnitrate, ntammonium, ntalkatm , ntanthpco2

# endif
      common /ncvars/ max_blk, ncidbulk, nrst, ncidrst, nrecrst,
     &      nrrec, nrpfrst, ncidclm,nwrt, ncidhis, nrechis,
     &      nrpfhis
# ifdef BIOLOGY_BEC
     &     , ntdust, ntiron , ntnitrate, ntammonium, ntalkatm , ntanthpco2
# endif
# ifdef BULK_DIAGS
      integer hisSustr_blk,   hisSvstr_blk
     &      , hisShflx_net,   hisShflx_rad
     &      , hisSwflx_emp
     &      , hisShflx_lat,   hisShflx_sen
     &      , hisShflx_wwk
     &      , hisSurf_u, hisSurf_v
      common /ncvars/ hisSustr_blk,   hisSvstr_blk
     &      , hisShflx_net,   hisShflx_rad
     &      , hisSwflx_emp
     &      , hisShflx_lat,    hisShflx_sen
     &      , hisShflx_wwk
     &      , hisSurf_u, hisSurf_v
# endif

#else /* BULK_FLUX */
      integer, parameter :: max_frc_file=9
      integer max_frc, ncidfrc(max_frc_file), nrst, ncidrst, nrecrst,
     &      nrrec, nrpfrst, ncidclm, nwrt, ncidhis, nrechis, nrpfhis
# ifdef BIOLOGY_BEC
     &     , ntdust, ntiron , ntnitrate, ntammonium, ntalkatm , ntanthpco2
# endif
      common /ncvars/       max_frc, ncidfrc, nrst, ncidrst, nrecrst,
     &      nrrec, nrpfrst, ncidclm, nwrt, ncidhis, nrechis, nrpfhis
# ifdef BIOLOGY_BEC
     &     , ntdust, ntiron , ntnitrate, ntammonium, ntalkatm , ntanthpco2
# endif
#endif /* BULK_FLUX */

#ifdef AVERAGES
      integer ntsavg,  navg
# ifdef OUTPUT_AVG_PHYS_BIO
     &      , nts_bio_avg, n_bio_avg
# endif
# ifdef OUTPUT_SST_AVG
     &      , nts_sst_avg, n_sst_avg
# endif
      common /ncvars/ ntsavg,  navg
# ifdef OUTPUT_AVG_PHYS_BIO
     &      , nts_bio_avg, n_bio_avg
# endif
# ifdef OUTPUT_SST_AVG
     &      , nts_sst_avg, n_sst_avg
# endif
#endif /* AVERAGES */
#ifdef STATIONS
      integer nsta
      common /ncvars/ nsta
#endif
#ifdef FLOATS
      integer nflt
      common /ncvars/ nflt
#endif
      

# ifdef DIAGNOSTICS_EDDY
      integer nwrtdiags_eddy
      common /ncvars/ nwrtdiags_eddy
#ifdef AVERAGES
      integer ntsdiags_eddy_avg, nwrtdiags_eddy_avg
      common /ncvars/ ntsdiags_eddy_avg, nwrtdiags_eddy_avg
#endif
#endif
! NetCFD ids for model variables 

      integer rstTime, rstTstep,      rstZ,   rstUb,  rstVb,
     &        hisTime, hisTstep,      hisZ,   hisUb,  hisVb
      common /ncvars/
     &        rstTime, rstTstep,      rstZ,   rstUb,  rstVb,
     &        hisTime, hisTstep,      hisZ,   hisUb,  hisVb
#ifdef SOLVE3D
# ifdef EXACT_RESTART
      integer rst_DU_avg2, rst_DV_avg2
      common /ncvars/ rst_DU_avg2, rst_DV_avg2
# endif
      integer rstU, rstV, rstT(NT+1), hisO,   hisW,   hisR,
     &        hisU, hisV, hisT(NT+1), hisAkv, hisAkt, hisAks
# if defined BIOLOGY_NPZDOC || defined BIOLOGY_BEC
     &      , rstPH, rstPCO2, rstPCO2air, rstPARinc, rstPAR
     &      , hisPH, hisPCO2, hisPCO2air, hisPARinc, hisPAR
     &      , avgPH, avgPCO2, avgPCO2air, avgPARinc, avgPAR
#  if defined AVERAGES && defined OUTPUT_AVG_PHYS_BIO
     &     , rstTstepBA
#  endif
# endif /* BIOLOGY_NPZDOC || BIOLOGY_BEC */
# if defined AVERAGES && defined OUTPUT_SST_AVG
     &     , rstTstepSA
# endif
# if defined BGC_FLUX_ANALYSIS || defined PHYS_FLUX_ANALYSIS
     &      , rstTstepFA
# endif
# ifdef SEDIMENT_BIOLOGY
     &      , rstTsed(NT_sed), hisTsed(NT_sed)
# endif /* SEDIMENT_BIOLOGY */
      common /ncvars/
     &        rstU, rstV, rstT,       hisO,   hisW,   hisR,
     &        hisU, hisV, hisT,       hisAkv, hisAkt, hisAks
# if defined BIOLOGY_NPZDOC || defined BIOLOGY_BEC
     &      , rstPH, rstPCO2, rstPCO2air, rstPARinc, rstPAR
     &      , hisPH, hisPCO2, hisPCO2air, hisPARinc, hisPAR
     &      , avgPH, avgPCO2, avgPCO2air, avgPARinc, avgPAR
#  if defined AVERAGES && defined OUTPUT_AVG_PHYS_BIO
     &     , rstTstepBA
#  endif
# endif /* BIOLOGY_NPZDOC || BIOLOGY_BEC */
# if defined AVERAGES && defined OUTPUT_SST_AVG
     &     , rstTstepSA
# endif
# if defined BGC_FLUX_ANALYSIS || defined PHYS_FLUX_ANALYSIS
     &      , rstTstepFA
# endif
# ifdef SEDIMENT_BIOLOGY
     &      , rstTsed, hisTsed
# endif /* SEDIMENT_BIOLOGY */

# ifdef LMD_KPP
      integer rstHbls, hisHbls
      common /ncvars/ rstHbls, hisHbls
# endif
# ifdef LMD_BKPP
      integer rstHbbl, hisHbbl
      common /ncvars/ rstHbbl, hisHbbl
# endif
# ifdef BULK_DIAGS
      integer avgSustr_blk,   avgSvstr_blk
     &      , avgShflx_net,   avgShflx_rad
     &      , avgSwflx_emp
     &      , avgShflx_lat,   avgShflx_sen
     &      , avgShflx_wwk
     &      , avgSurf_u, avgSurf_v
      common /ncvars/ avgSustr_blk,   avgSvstr_blk
     &      , avgShflx_net,   avgShflx_rad
     &      , avgSwflx_emp
     &      , avgShflx_lat,   avgShflx_sen
     &      , avgShflx_wwk
     &      , avgSurf_u, avgSurf_v
# endif

# ifdef DIAGNOSTICS_EDDY
      real timediags_eddy
      common /ncvars/
     &        timediags_eddy


      integer nciddiags_eddy, nrecdiags_eddy, nrpfdiags_eddy
     &      , diags_eddyTime, diags_eddyTstep
     &      , diags_eddyuu(2), diags_eddyvv(2), diags_eddyuv(2)
     &      , diags_eddyub(2), diags_eddyvb(2), diags_eddywb(2)
     &      , diags_eddynotused1(2), diags_eddynotused2(2)
     &      , diags_eddynotused3(2)
      common /ncvars/   nciddiags_eddy, nrecdiags_eddy, nrpfdiags_eddy
     &      , diags_eddyTime, diags_eddyTstep
     &      , diags_eddyuu, diags_eddyvv, diags_eddyuv, diags_eddyub
     &      , diags_eddyvb, diags_eddywb, diags_eddynotused1, diags_eddynotused2
     &      , diags_eddynotused3
# endif

#endif

#ifdef AVERAGES
      integer ncidavg, nrecavg, nrpfavg, avgTime, avgTstep, avgZ, avgUb, 
     &     avgVb
      common /ncvars/  ncidavg, nrecavg, nrpfavg,
     &                                   avgTime, avgTstep, avgZ, avgUb, 
     &     avgVb
# ifdef SOLVE3D


#  ifdef DIAGNOSTICS_EDDY
      real timediags_eddy_avg
      common /ncvars/
     &        timediags_eddy_avg

       integer nciddiags_eddy_avg, nrecdiags_eddy_avg
     &      , nrpfdiags_eddy_avg
     &      , diags_eddy_avgTime, diags_eddy_avgTstep
     &      , diags_eddyuu_avg(2), diags_eddyvv_avg(2), diags_eddyuv_avg(2)
     &      , diags_eddyub_avg(2), diags_eddyvb_avg(2), diags_eddywb_avg(2)
     &      , diags_eddynotused1_avg(2), diags_eddynotused2_avg(2)
     &      , diags_eddynotused3_avg(2)
      common /ncvars/   nciddiags_eddy_avg, nrecdiags_eddy_avg
     &      , nrpfdiags_eddy_avg
     &      , diags_eddy_avgTime, diags_eddy_avgTstep
     &      , diags_eddyuu_avg, diags_eddyvv_avg, diags_eddyuv_avg
     &      , diags_eddyub_avg, diags_eddyvb_avg, diags_eddywb_avg
     &      , diags_eddynotused1_avg, diags_eddynotused2_avg
     &      , diags_eddynotused3_avg
#  endif


#  ifdef OUTPUT_SST_AVG
      integer ncid_sst_avg, nrec_sst_avg, nrpf_sst_avg
     &     , avg_sst_time, avg_sst_Tstep, avgSST
      common /ncvars_sst/ ncid_sst_avg, nrec_sst_avg, nrpf_sst_avg
     &     , avg_sst_time, avg_sst_Tstep, avgSST
#  endif
#  ifdef OUTPUT_AVG_PHYS_BIO
      integer ncid_bio_avg, nrec_bio_avg, nrpf_bio_avg
     &     , avg_bio_time, avg_bio_Tstep, avg_bio_Z
      common /ncvars_bio/ ncid_bio_avg, nrec_bio_avg, nrpf_bio_avg
     &     , avg_bio_time, avg_bio_Tstep, avg_bio_Z
#  endif
      integer avgU,  avgV,  avgT(NT+1),  avgR,    avgO,    avgW,
     &                                   avgAkv,  avgAkt,  avgAks
#  ifdef SEDIMENT_BIOLOGY
     &      , avgTsed(NT_sed)
#  endif /* SEDIMENT_BIOLOGY */
      common /ncvars/ avgU, avgV, avgT,  avgR,    avgO,    avgW,
     &                                   avgAkv,  avgAkt,  avgAks
#  ifdef SEDIMENT_BIOLOGY
     &      , avgTsed
#  endif /* SEDIMENT_BIOLOGY */
#  ifdef OUTPUT_SPECDOM
      integer avgFreqDomSP_sfc, avgFreqDomDIAT_sfc, 
     &     avgFreqDomDIAZ_sfc, avgFreqDomSP_int,
     &     avgFreqDomDIAT_int, avgFreqDomDIAZ_int
      common /ncvars/ avgFreqDomSP_sfc, avgFreqDomDIAT_sfc, 
     &     avgFreqDomDIAZ_sfc, avgFreqDomSP_int,
     &     avgFreqDomDIAT_int, avgFreqDomDIAZ_int
#  endif /* OUTPUT_SPECDOM */
#  ifdef LMD_KPP
      integer avgHbls
      common /ncvars/ avgHbls
#  endif
#  ifdef LMD_BKPP
      integer avgHbbl
      common /ncvars/ avgHbbl
#  endif
# endif
#endif

#ifdef STATIONS
      integer nstation,  ispos(NS), jspos(NS),
     &        ncidstn, nrecstn, nrpfstn, stnTime, stnZ, stnUb, stnVb
      common /ncvars/ nstation, ispos,   jspos,
     &        ncidstn, nrecstn, nrpfstn, stnTime, stnZ, stnUb, stnVb
# ifdef SOLVE3D
      integer stnU,  stnV,  stnT(NT+1),  stnR,    stnO,    stnW,
     &                                   stnAkv,  stnAkt,  stnAks
      common /ncvars/ stnU, stnV, stnT,  stnR,    stnO,    stnW,
     &                                   stnAkv,  stnAkt,  stnAks
#  ifdef LMD_KPP
      integer stnHbls 
      common /ncvars/ stnHbls 
#  endif
#  ifdef LMD_BKPP
      integer stnHbbl 
      common /ncvars/ stnHbbl 
#  endif
# endif
#endif

#ifdef SOLVE3D
# define NWRTHIS 100+NT
#else
# define NWRTHIS 14      
#endif
      logical ldefhis, wrthis(NWRTHIS)
      common /ncvars/ ldefhis, wrthis
#ifdef DIAGNOSTICS_EDDY
      logical ldefdiags_eddy, wrtdiags_eddy(3)
      common /ncvars/ ldefdiags_eddy, wrtdiags_eddy
#endif
#ifdef AVERAGES
      logical wrtavg(NWRTHIS)
      common /ncvars/ wrtavg
# ifdef DIAGNOSTICS_EDDY
      logical ldefdiags_eddy_avg, wrtdiags_eddy_avg(3)
      common /ncvars/ ldefdiags_eddy_avg, wrtdiags_eddy_avg
# endif
#endif
# if defined AVERAGES && defined PHYS_FLUX_ANALYSIS && defined PFA_BY_TRACER
      logical wrt_pfa(NT)
      common /ncvars/ wrt_pfa
# endif /* AVERAGES && PHYS_FLUX_ANALYSIS && PFA_BY_TRACER */
#ifdef FLOATS
      logical ldefflt
      common /ncvars/ ldefflt
#endif
#ifdef STATIONS
      logical wrtsta(NWRTHIS)
      common /ncvars/ wrtsta
#endif


! Grid Type Codes:  r2dvar....w3hvar are codes for netCDF array types
!----- ---- ------  which are set according to the rule:
!
!              type = horiz_grid_type + 4*vert_grid_type
!
! where horiz_grid_type=0,1,2,3 for RHO-,U-,V-,PSI-points respectively,
! and vert_grid_type=0 for 2D fields; 1 for vertical RHO-points, 2 for
! vertical W-points.

      integer, parameter :: r2dvar=0, u2dvar=1, v2dvar=2, p2dvar=3,
     &            r3dvar=4, u3dvar=5, v3dvar=6, p3dvar=7, w3dvar=8

!            Horizontal array dimensions in netCDF files. In the case 
! xi_rho     of MPI code with PARALLEL_FILES activated these dimensions
! xi_u       depend on the corresponding sizes of individual subdomains
! eta_rho    rather than the whole physical grid, therefore they become 
! eta_v      live variables placed in common block and assigned values
!            in "mpi_setup" rather tnan be parameters defined here.

#if defined MPI && defined PARALLEL_FILES
      integer xi_rho, xi_u,   eta_rho, eta_v
      common /ncvars/ xi_rho, xi_u,   eta_rho, eta_v
#else
      integer, parameter :: xi_rho=LLm+2, eta_rho=MMm+2,
     &                     xi_u=xi_rho-1, eta_v=eta_rho-1
#endif

      integer, parameter :: max_name_size=64
      character date_str*44, title*80
      character*(max_name_size) ininame, grdname,
     &                 hisname, rstname
#if defined BULK_FLUX     
     &    , blkfile(max_blk_file)
#else
     &    , frcfile(max_frc_file)
# endif

      common /cncvars/ date_str, title,  ininame,
     &        grdname, hisname, rstname 
#if defined BULK_FLUX     
     & , blkfile
#else
     & , frcfile
#endif
#ifdef AVERAGES
      character*(max_name_size) avgname
      common /cncvars/ avgname
#ifdef OUTPUT_AVG_PHYS_BIO
      character*(max_name_size) bio_avgname
      common /cncvars_avg_bio/ bio_avgname
#endif
#ifdef OUTPUT_SST_AVG
      character*(max_name_size) sst_avgname
      common /cncvars_sst_avg/ sst_avgname
#endif
#endif /* AVERAGES */
#ifdef DIAGNOSTICS_EDDY
      character*(max_name_size) diags_eddyname
      common /cncvars/ diags_eddyname
# ifdef AVERAGES
      character*(max_name_size) diags_eddy_avgname
      common /cncvars/ diags_eddy_avgname
# endif
#endif
#if (defined TCLIMATOLOGY && !defined ANA_TCLIMA) || !defined ANA_SSH
      character*(max_name_size) clm_file
      common /cncvars/ clm_file
#endif
#if defined T_FRC_BRY  || defined M2_FRC_BRY || \
    defined M3_FRC_BRY || defined Z_FRC_BRY
      character*(max_name_size) bry_file 
      common /cncvars/ bry_file
#endif

#if (defined PSOURCE || defined Q_PSOURCE) && !defined ANA_PSOURCE
      character*(max_name_size) psrc_file
      common /source_file/ psrc_file
#endif

#ifdef STATIONS
      character*(max_name_size) staname
      common /cncvars/ staname
#endif
#ifdef ASSIMILATION
      character*(max_name_size) aparnam, assname
      common /cncvars/ aparnam, assname
#endif
      character*42  vname(3,45+NT
#ifdef SEDIMENT_BIOLOGY
     &     + NT_sed
#endif
     &     )
      common /cncvars/ vname
