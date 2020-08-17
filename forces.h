! Surface momentum flux (wind stress):
!======== ======== ==== ===== ========
! sustr   XI- and ETA-components of kinematic surface momentum
! svstr   flux (wind stresses) at defined horizontal U- and
!         V-points,  dimensioned as [m^2/s^2].
! tsms    time of surface momentum stresses.
! sustrg  two-time level gridded data for XI- anf ETA-componets
! svstrg  of wind stess (normally assumed to be in [Newton/m^2].

      real sustr(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE sustr(BLOCK_PATTERN) BLOCK_CLAUSE
      real svstr(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE svstr(BLOCK_PATTERN) BLOCK_CLAUSE
      common /forces_sustr/sustr /forces_svstr/svstr
#ifdef WIND_STRESS_AMP
      real amptau(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE amptau(BLOCK_PATTERN) BLOCK_CLAUSE
      common /forces_amptau/ amptau
#endif

#ifndef BULK_FLUX 


# ifndef ANA_SMFLUX
#  if defined SMFLUX_DATA || defined ALL_DATA
#   undef SMFLUX_DATA
      real sustrg(GLOBAL_2D_ARRAY,2)
CSDISTRIBUTE_RESHAPE sustrg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real svstrg(GLOBAL_2D_ARRAY,2)
CSDISTRIBUTE_RESHAPE svstrg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      common /smsdat_sustrg/sustrg /smsdat_svstrg/svstrg
#   ifdef WIND_STRESS_AMP
      real amptaug(GLOBAL_2D_ARRAY,2)
CSDISTRIBUTE_RESHAPE amptau(BLOCK_PATTERN) BLOCK_CLAUSE
      common /forces_amptaug/ amptaug
#   endif
 
      real sms_cycle, sms_time(2)
      integer sms_ncycle,   sms_rec,   itsms, ntsms,
     &        sms_file_id,  sms_tid,   susid, svsid
      common /smsdat/ sms_cycle,       sms_time,
     &        sms_ncycle,   sms_rec,   itsms, ntsms,
     &        sms_file_id,  sms_tid,   susid, svsid

#   ifdef WIND_STRESS_AMP
      real atau_cycle, atau_time(2)
      integer atau_ncycle,  atau_rec,  itatau, ntatau,
     &        atau_file_id, atau_tid,  atauid
      common /smsdat/ atau_cycle,      atau_time,
     &        atau_ncycle,  atau_rec,  itatau, ntatau,
     &        atau_file_id, atau_tid,  atauid
#   endif
#  endif /* SMFLUX_DATA */
# endif /* !ANA_SMFLUX */
#else
      integer ntrad,nttra,ntprec,ntwnd
      common /blk_nt/ntrad,nttra,ntprec,ntwnd

      integer tra_file_id, rad_file_id, prec_file_id,
     &        wnd_file_id
      common /blk_id/ tra_file_id,rad_file_id,
     &                prec_file_id,wnd_file_id
#endif


#ifdef SOLVE3D

! Surface tracer fluxes:
!======== ====== =======
!  stflx   kinematic surface fluxes of tracer type variables at
!          horizontal RHO-points. Physical dimensions [degC m/s]
!          - temperature; [PSU m/s] - salinity.
!  stflxg  two-time level surface tracer flux grided data.
!  tstflx  time of surface tracer flux.

      real stflx(GLOBAL_2D_ARRAY,NT)
CSDISTRIBUTE_RESHAPE stflx(BLOCK_PATTERN,*) BLOCK_CLAUSE
      common /forces_stflx/stflx
# ifndef BULK_FLUX

#  if !defined ANA_STFLUX || !defined ANA_SSFLUX
#   if defined STFLUX_DATA || defined ALL_DATA
#    undef STFLUX_DATA

      real stflxg(GLOBAL_2D_ARRAY,2,NT)
CSDISTRIBUTE_RESHAPE stflxg(BLOCK_PATTERN,*,*) BLOCK_CLAUSE
      common /stfdat_stflxg/stflxg
 
      real stf_cycle(NT),      stf_time(2,NT)
      integer stf_ncycle(NT),  stf_rec(NT), itstf(NT),  ntstf(NT),
     &        stf_file_id(NT), stf_id(NT),  stf_tid(NT)

      common /stfdat/          stf_cycle,   stf_time,
     &        stf_ncycle,      stf_rec,     itstf,      ntstf,
     &        stf_file_id,     stf_id,      stf_tid
 
#   endif /*  STFLUX_DATA */
#  endif /* !ANA_STFLUX || !ANA_SSFLUX */
 
#  if defined QCORRECTION && !defined ANA_SST
#   if defined SST_DATA || defined ALL_DATA
#    undef SST_DATA

! Heat flux correction:
!===== ==== ===========
! sstg    two-time-level  grided data for sea surface temperature
! dqdtg   SST[deg C] and net surface heat flux sensitivity to SST
!         dQdSST [Watts/m^2/Celsius]
! tsst    time of sea surface temperature data.

      real sstg(GLOBAL_2D_ARRAY,2)
CSDISTRIBUTE_RESHAPE  sstg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real dqdtg(GLOBAL_2D_ARRAY,2)
CSDISTRIBUTE_RESHAPE dqdtg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real sssg(GLOBAL_2D_ARRAY,2)
CSDISTRIBUTE_RESHAPE  sssg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      common /sst_dat/sstg /dqdtg_dat/dqdtg /sss_dat/sssg
 
      real sst_cycle, sst_time(2), sss_cycle,  sss_time(2),
     &                             dqdt_cycle, dqdt_time(2)

      integer sst_ncycle,   sst_rec,  itsst,    ntsst,
     &        sst_file_id,  sst_id,   sst_tid,
     &        dqdt_ncycle,  dqdt_rec, itdqdt,   ntdqdt,
     &        dqdt_file_id, dqdt_id,  dqdt_tid,
     &        sss_ncycle,   sss_rec,  itsss,    ntsss,
     &        sss_file_id,  sss_id,   sss_tid

      common /qcorr_data/ sst_cycle,  sst_time, sss_cycle,
     &        sss_time,   dqdt_cycle, dqdt_time,
     &        sst_ncycle,   sst_rec,  itsst,    ntsst,
     &        sst_file_id,  sst_id,   sst_tid,
     &        dqdt_ncycle,  dqdt_rec, itdqdt,   ntdqdt,
     &        dqdt_file_id, dqdt_id,  dqdt_tid,
     &        sss_ncycle,   sss_rec,  itsss,    ntsss, 
     &        sss_file_id,  sss_id,   sss_tid

#  endif /* SST_DATA */
# endif /* QCORRECTION && !ANA_SST */
#else

      real ust(GLOBAL_2D_ARRAY,2)
CSDISTRIBUTE_RESHAPE ust(BLOCK_PATTERN,*) BLOCK_CLAUSE
      common /bulk_ust/ust
      real tst(GLOBAL_2D_ARRAY,2)
CSDISTRIBUTE_RESHAPE tst_sl(BLOCK_PATTERN,*) BLOCK_CLAUSE
      common /bulk_tst/tst
      real qst(GLOBAL_2D_ARRAY,2)
CSDISTRIBUTE_RESHAPE qst(BLOCK_PATTERN,*) BLOCK_CLAUSE
      common /bulk_qst/qst


      real prec_scale,tra_scale,wnd_scale,srf_scale
      common /blk_scale/prec_scale,tra_scale,wnd_scale,srf_scale

!
!  BULK FORMULATION
!--------------------------------------------------------------------

      real tair(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE tair(BLOCK_PATTERN) BLOCK_CLAUSE
      real qair(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE qair(BLOCK_PATTERN) BLOCK_CLAUSE
      real rain(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE rain(BLOCK_PATTERN) BLOCK_CLAUSE
      real radlw(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE radlw(BLOCK_PATTERN) BLOCK_CLAUSE
      real radsw(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE radsw(BLOCK_PATTERN) BLOCK_CLAUSE
      real uwnd(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE uwnd(BLOCK_PATTERN) BLOCK_CLAUSE
      real vwnd(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE vwnd(BLOCK_PATTERN) BLOCK_CLAUSE

      common /bulk_tair/tair
     &       /bulk_qair/qair
     &       /bulk_rain/rain
     &       /bulk_radlw/radlw
     &       /bulk_radsw/radsw
     &       /bulk_uwnd/uwnd
     &       /bulk_vwnd/vwnd

      real tairg(GLOBAL_2D_ARRAY,2)
CSDISTRIBUTE_RESHAPE tairg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real qairg(GLOBAL_2D_ARRAY,2)
CSDISTRIBUTE_RESHAPE qairg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real raing(GLOBAL_2D_ARRAY,2)
CSDISTRIBUTE_RESHAPE raing(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real radlwg(GLOBAL_2D_ARRAY,2)
CSDISTRIBUTE_RESHAPE radlwg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real radswg(GLOBAL_2D_ARRAY,2)
CSDISTRIBUTE_RESHAPE radswg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real uwndg(GLOBAL_2D_ARRAY,2)
CSDISTRIBUTE_RESHAPE uwndg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real vwndg(GLOBAL_2D_ARRAY,2)
CSDISTRIBUTE_RESHAPE vwndg(BLOCK_PATTERN,*) BLOCK_CLAUSE

      common /bulkdat_tairg/tairg
     &       /bulkdat_qairg/qairg
     &       /bulkdat_raing/raing
     &       /bulkdat_radlwg/radlwg
     &       /bulkdat_radswg/radswg
     &       /bulk_uwndg/uwndg
     &       /bulk_vwndg/vwndg

      real    rad_time(2), rad_cycle
      real    tra_time(2), tra_cycle
      real    wnd_time(2), wnd_cycle
      real    prec_time(2),prec_cycle
      real    ztref, zref

      common /bulkdat2/
     &        rad_time, rad_cycle,
     &        tra_time, tra_cycle,
     &        wnd_time, wnd_cycle,
     &        prec_time, prec_cycle
     &        ,ztref, zref


      integer tair_id,uwnd_id,radsw_id,prec_id
      integer qair_id,vwnd_id,radlw_id

      integer itrad,rad_ncycle,rad_rec,rad_tid
      integer ittra,tra_ncycle,tra_rec,tra_tid
      integer itwnd,wnd_ncycle,wnd_rec,wnd_tid
      integer itprec,prec_ncycle,prec_rec,prec_tid

      common /bulkdat1/
     &        tair_id,uwnd_id,radsw_id,prec_id,
     &        qair_id,vwnd_id,radlw_id,
     &        itrad,rad_ncycle,rad_rec,rad_tid,
     &        ittra,tra_ncycle,tra_rec,tra_tid,
     &        itwnd,wnd_ncycle,wnd_rec,wnd_tid,
     &        itprec,prec_ncycle,prec_rec,prec_tid
#  if defined SALINITY && defined SALINITY_RESTORING
      real sssg(GLOBAL_2D_ARRAY,2)
CSDISTRIBUTE_RESHAPE  sssg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      common /sss_dat/sssg
      real sss_cycle, sss_time(2)
      integer sss_ncycle,   sss_rec,  itsss,    ntsss,
     &        sss_file_id,  sss_id,   sss_tid
      common /sssrest_data/ sss_cycle,sss_time,sss_ncycle,
     &  sss_rec,itsss,ntsss,sss_file_id,sss_id, sss_tid
#  endif

# endif /* BULK_FLUX */



#ifdef BULK_DIAGS
      real sustr_blk(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE sustr_blk(BLOCK_PATTERN) BLOCK_CLAUSE
      common /frc_sustr_blk/sustr_blk

      real svstr_blk(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE svstr_blk(BLOCK_PATTERN) BLOCK_CLAUSE
      common /frc_svstr_blk/svstr_blk

      real shflx_net(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE shflx_net(BLOCK_PATTERN) BLOCK_CLAUSE
      common /frc_shflx_net/shflx_net
/* Add Latent and Sensible heat fluxes output */
      real shflx_lat(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE shflx_lat(BLOCK_PATTERN) BLOCK_CLAUSE
      common /frc_shflx_lat/shflx_lat

      real shflx_sen(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE shflx_sen(BLOCK_PATTERN) BLOCK_CLAUSE
      common /frc_shflx_sen/shflx_sen
/* */

      real shflx_rad(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE shflx_rad(BLOCK_PATTERN) BLOCK_CLAUSE
      common /frc_shflx_rad/shflx_rad

      real swflx_emp(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE swflx_emp(BLOCK_PATTERN) BLOCK_CLAUSE
      common /frc_swflx_emp/swflx_emp

      real shflx_wwk(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE shflx_wwk(BLOCK_PATTERN) BLOCK_CLAUSE
      common /frc_shflx_wwk/shflx_wwk

      real surf_u(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE surf_u(BLOCK_PATTERN) BLOCK_CLAUSE
      common /frc_surf_u/surf_u

      real surf_v(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE surf_v(BLOCK_PATTERN) BLOCK_CLAUSE
      common /frc_surf_v/surf_v

#ifdef BULK_DIAGS
      real time_bulk_diags_his
      common /t_bulk_diags_his/time_bulk_diags_his

      real zeta_bulk_diags_his(GLOBAL_2D_ARRAY)
      common /zeta_bulk_diags_his/zeta_bulk_diags_his


      logical new_bulk_diags_his
      integer n_bulk_diags_his
      common /nc_bulk_diags_his/n_bulk_diags_his,new_bulk_diags_his

      integer nrpf_bulk_diags_his,ncid_bulk_diags_his
     &       , nrec_bulk_diags_his,bulk_diags_hisTstep, bulk_diags_hisZ
     &       , bulk_diags_hisTime

      common /nc_bulk_diags_his/ nrpf_bulk_diags_his
     &     , ncid_bulk_diags_his, nrec_bulk_diags_his
     &     , bulk_diags_hisTstep, bulk_diags_hisZ
     &     , bulk_diags_hisTime

      character*80 bulk_diags_his_name
      common /nc_bulk_diags_his/bulk_diags_his_name


# ifdef AVERAGES
      real time_bulk_diags_avg
      common /t_bulk_diags_avg/time_bulk_diags_avg

      real zeta_bulk_diags_avg(GLOBAL_2D_ARRAY)
      common /zeta_bulk_diags_avg/zeta_bulk_diags_avg


      logical new_bulk_diags_avg
      integer nts_bulk_diags_avg, n_bulk_diags_avg
      common /nc_bulk_diags_avg/nts_bulk_diags_avg
     $        ,n_bulk_diags_avg,new_bulk_diags_avg

      integer nrpf_bulk_diags_avg,ncid_bulk_diags_avg
     &       , nrec_bulk_diags_avg,bulk_diags_avgTstep, bulk_diags_avgZ
     &       , bulk_diags_avgTime

      common /nc_bulk_diags_avg/ nrpf_bulk_diags_avg
     &     , ncid_bulk_diags_avg, nrec_bulk_diags_avg
     &     , bulk_diags_avgTstep, bulk_diags_avgZ
     &     , bulk_diags_avgTime

      character*80 bulk_diags_avg_name
      common /nc_bulk_diags_avg/bulk_diags_avg_name
# endif

#endif

#endif




#ifdef BIOLOGY_BEC
! dust flux

         real dust(GLOBAL_2D_ARRAY) 
CSDISTRIBUTE_RESHAPE  dust(BLOCK_PATTERN,*) BLOCK_CLAUSE
           common /forces_dust/dust
       real dustg(GLOBAL_2D_ARRAY,2)
CSDISTRIBUTE_RESHAPE  dustg(BLOCK_PATTERN,*) BLOCK_CLAUSE
       common /dustdat_dustg/dustg

        real dustp(2), dust_time(2),dust_cycle, scldqdt
        integer itdust,dust_id,ldustgrd ,dust_ncycle,
     &  dust_rec,dust_tid,dust_file_id,iron_file_id,
     &  nitrate_file_id,ammonium_file_id,alkatm_file_id,
     &  anthpco2_file_id
       common/dustdat/itdust,dust_id,ldustgrd,
     &  dust_ncycle,dust_rec,dust_tid,dust_file_id,iron_file_id,
     &  nitrate_file_id,ammonium_file_id,alkatm_file_id,anthpco2_file_id
       common/dustdat1/dustp,dust_time,dust_cycle,scldqdt


! iron flux
       real iron(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE  iron(BLOCK_PATTERN,*) BLOCK_CLAUSE
       common /forces_iron/iron
       real irong(GLOBAL_2D_ARRAY,2)
CSDISTRIBUTE_RESHAPE irong(BLOCK_PATTERN,*) BLOCK_CLAUSE
       common /irondat_irong/irong
     
       real ironp(2),iron_time(2),iron_cycle
       integer itiron,iron_id,lirongrd,iron_ncycle,
     &  iron_rec,iron_tid

       common/irondat/ironp,iron_time,iron_cycle
       common/irondat1/itiron,iron_id,lirongrd,
     &  iron_ncycle,iron_rec,iron_tid

! nitrate flux
       real nitrate(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE  nitrate(BLOCK_PATTERN,*) BLOCK_CLAUSE
       common /forces_nitrate/nitrate
       real nitrateg(GLOBAL_2D_ARRAY,2)
CSDISTRIBUTE_RESHAPE nitrateg(BLOCK_PATTERN,*) BLOCK_CLAUSE
       common /nitratedat_nitrateg/nitrateg

       real nitratep(2),nitrate_time(2),nitrate_cycle
       integer itnitrate,nitrate_id,lnitrategrd,nitrate_ncycle,
     &  nitrate_rec,nitrate_tid

       common/nitratedat/nitratep,nitrate_time,nitrate_cycle
       common/nitratedat1/itnitrate,nitrate_id,lnitrategrd,
     &  nitrate_ncycle,nitrate_rec,nitrate_tid


! ammonium flux
       real ammonium(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE  ammonium(BLOCK_PATTERN,*) BLOCK_CLAUSE
       common /forces_ammonium/ammonium
       real ammoniumg(GLOBAL_2D_ARRAY,2)
CSDISTRIBUTE_RESHAPE ammoniumg(BLOCK_PATTERN,*) BLOCK_CLAUSE
       common /ammoniumdat_ammoniumg/ammoniumg

       real ammoniump(2),ammonium_time(2),ammonium_cycle
       integer itammonium,ammonium_id,lammoniumgrd,ammonium_ncycle,
     &  ammonium_rec,ammonium_tid

       common/ammoniumdat/ammoniump,ammonium_time,ammonium_cycle
       common/ammoniumdat1/itammonium,ammonium_id,lammoniumgrd,
     &  ammonium_ncycle,ammonium_rec,ammonium_tid

! alkatm flux
       real alkatm(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE  alkatm(BLOCK_PATTERN,*) BLOCK_CLAUSE
       common /forces_alkatm/alkatm
       real alkatmg(GLOBAL_2D_ARRAY,2)
CSDISTRIBUTE_RESHAPE alkatmg(BLOCK_PATTERN,*) BLOCK_CLAUSE
       common /alkatmdat_alkatmg/alkatmg

       real alkatmp(2),alkatm_time(2),alkatm_cycle
       integer italkatm,alkatm_id,lalkatmgrd,alkatm_ncycle,
     &  alkatm_rec,alkatm_tid

       common/alkatmdat/alkatmp,alkatm_time,alkatm_cycle
       common/alkatmdat1/italkatm,alkatm_id,lalkatmgrd,
     &  alkatm_ncycle,alkatm_rec,alkatm_tid


! anthpco2 flux
       real anthpco2(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE  anthpco2(BLOCK_PATTERN,*) BLOCK_CLAUSE
       common /forces_anthpco2/anthpco2
       real anthpco2g(GLOBAL_2D_ARRAY,2)
CSDISTRIBUTE_RESHAPE anthpco2g(BLOCK_PATTERN,*) BLOCK_CLAUSE
       common /anthpco2dat_anthpco2g/anthpco2g

       real anthpco2p(2),anthpco2_time(2),anthpco2_cycle
       integer itanthpco2,anthpco2_id,lanthpco2grd,anthpco2_ncycle,
     &  anthpco2_rec,anthpco2_tid

       common/anthpco2dat/anthpco2p,anthpco2_time,anthpco2_cycle
       common/anthpco2dat1/itanthpco2,anthpco2_id,lanthpco2grd,
     &  anthpco2_ncycle,anthpco2_rec,anthpco2_tid


#endif


!
! Solar short-wave radiation flux:
!====== ===== ==== ========= =====
! srflx   kinematic surface shortwave solar radiation flux
!              in [degC m/s] at horizontal RHO-points
! srflxg  two-time-level grided data for surface [Watts/m^2]
! tsrflx  time of solar shortwave radiation flux.

      real srflx(GLOBAL_2D_ARRAY)
#if defined DAILYPAR_PHOTOINHIBITION || defined DAILYPAR_BEC
      real srflx_dailyavg(GLOBAL_2D_ARRAY)
# ifndef DIURNAL_SRFLUX
! additional variables to re-create daily average from diurnal variations

! value should be 24 / "number of hours between new radsw inputs"
      integer, parameter :: n_srflx_day = 24
      integer :: num_srflx_day ! will be set to n_srflx_day * (3600/dt (in sec))
!      integer, parameter :: MAX_NUM_SRFLX_DAY = 576 ! WARNING: dt ! >= 600 sec
      integer, parameter :: MAX_NUM_SRFLX_DAY = 30000 ! WARNING: dt ! >= 600 sec
      integer iptr_srflx_day(GLOBAL_2D_ARRAY)
      integer iptr_srflx_day_set(GLOBAL_2D_ARRAY)
      integer PARinc_rst_read
      real srflx_day(GLOBAL_2D_ARRAY,MAX_NUM_SRFLX_DAY) 
      real srflx_frac   ! will be set to: 1/num_srflx_day
# endif
#endif /* DAILYPAR_PHOTOINHIBITION || DAILYPAR_BEC */
CSDISTRIBUTE_RESHAPE srflx(BLOCK_PATTERN) BLOCK_CLAUSE
      common /forces_srflx/srflx
#if defined DAILYPAR_PHOTOINHIBITION || defined DAILYPAR_BEC
     &       , srflx_dailyavg
# ifndef DIURNAL_SRFLUX
     &     , srflx_day, srflx_frac 
      common /i_forces_srflx/ iptr_srflx_day, PARinc_rst_read,
     &     iptr_srflx_day_set,  num_srflx_day
# endif
#endif /* DAILYPAR_PHOTOINHIBITION || DAILYPAR_BEC */
# ifndef ANA_SRFLUX
#  if defined SRFLUX_DATA || defined ALL_DATA
#   undef SRFLUX_DATA

      real srflxg(GLOBAL_2D_ARRAY,2)
CSDISTRIBUTE_RESHAPE srflxg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      common /srfdat_srflxg/srflxg
 
      real srf_cycle, srf_time(2)
      integer srf_ncycle,  srf_rec, itsrf, ntsrf,
     &        srf_file_id, srf_tid, srf_id

      common /srfdat/ srf_cycle, srf_time,
     &        srf_ncycle,  srf_rec, itsrf, ntsrf,
     &        srf_file_id, srf_tid, srf_id
 
#  endif /* SRFLUX_DATA */
# endif /* !ANA_SRFLUX */
# if defined SG_BBL96 && !defined ANA_WWAVE
#  if defined WWAVE_DATA || defined ALL_DATA

!  WIND INDUCED WAVES:
!----------------------------------------------------------
!  wwag  |  Two-time-level       | wave amplitude [m]
!  wwdg  |  gridded data         | wave direction [radians]
!  wwpg  |  for wind induced     ! wave period [s]
!
!  wwap  |  Two-time-level       | wave amplitude [m]
!  wwdp  |  point data           | wave direction [radians]
!  wwpp  |  for wind induced     ! wave period [s]
!
!  tww      Time of wind induced waves.

      real wwag(GLOBAL_2D_ARRAY,2)
CSDISTRIBUTE_RESHAPE wwag(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real wwdg(GLOBAL_2D_ARRAY,2)
CSDISTRIBUTE_RESHAPE wwdg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real wwpg(GLOBAL_2D_ARRAY,2)
CSDISTRIBUTE_RESHAPE wwpg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      common /wwf_wwag/wwag /wwf_wwdg/wwdg /wwf_wwpg/wwpg
 
      real ww_tintrp(2), wwap(2), wwdp(2),  wwpp(2), tww(2), tsww,
     &        ww_tstart, ww_tend, sclwwa,   sclwwd,  sclwwp, wwclen
      integer itww,      twwindx, wwaid,    wwdid,   wwpid,  wwtid
      logical lwwgrd,    wwcycle, ww_onerec
      common /wwfdat/
     &        ww_tintrp, wwap,    wwdp,     wwpp,    tww,    tsww,
     &        ww_tstart, ww_tend, sclwwa,   sclwwd,  sclwwp, wwclen,
     &        itww,      twwindx, wwaid,    wwdid,   wwpid,  wwtid,
     &        lwwgrd,    wwcycle, ww_onerec
 
#   undef WWAVE_DATA
#  endif /* WWAVE_DATA */
# endif /* SG_BBL96 && !ANA_WWAVE */

# ifdef READ_ATM_PCO2
      real pCO2_atm
! value is read in
      common /pco2_atm/ pCO2_atm
# elif defined READ_VARIABLE_ATM_PCO2
      character(len=100) pco2_atm_file
      common /pco2_atm_file/ pco2_atm_file
# elif ! defined VARIABLE_ANN_ATM_PCO2 
! value is hard-coded (in ppm)
      real, parameter :: pCO2_atm = 360.0
# endif /* READ_ATM_PCO2 */
#endif /* SOLVE3D */
