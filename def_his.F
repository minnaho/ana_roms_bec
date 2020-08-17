#ifndef AVRH
# include "cppdefs.h"
 
      subroutine def_his (total_rec, ierr)
#else
      subroutine def_avg (total_rec, ierr)
#endif

! Create/open averages/history netCDF file. In the case when a new
! netCDF file is created, define all variables, their dimensions and
! attributes. In the case when a previously existing netCDF file is
! to be opened for addition of new data, verify that all dimensions
! of the file are consistent with the present model configuration
! and all necessary variables exist. Save netCDF IDs for all needed
! variables. Also determine size of the unlimited dimension.

! The difference between def_his and def_avg is as follows: they
! have different netCDF file name (hisname/avgname); netCDF file ID
! (passed as argument); time record index (hisindx/avgindx); array
! of switches which variables to write (wrthis/wrtavg); and different
! sets of netCDF variable IDs (hisTime...hisHbl/avgTime...avgHbl);
! and the first attribute of each variable, long_name, has prefix
! 'averaged'.  Because most of the code is identical for both
! routines, the second one is generated from the first entirely
! by CPP.

      implicit none
      logical create_new_file
      integer total_rec, ierr, rec, lncn,lvar,lenstr, timedim
     &      , r2dgrd(3), u2dgrd(3), v2dgrd(3), auxil(2), checkdims
#ifdef SOLVE3D
     &      , r3dgrd(4), u3dgrd(4), v3dgrd(4), w3dgrd(4), itrc
#endif
      integer my_nf_def_dim  ! <-- instrumented netCDF function.
      character*60 text

#include "param.h"
#include "scalars.h"
#include "ncvars.h"
#include "netcdf.inc"


#ifdef MASK_LAND_DATA
# ifdef HIS_DOUBLE
      real*8, parameter :: spv_set=1.D+33
# else
      real*4, parameter :: spv_set=1.E+33
# endif
#endif


#ifndef AVRH
# define ncid ncidhis
# define ncname hisname
# define rec_per_file nrpfhis
# define wrt wrthis
# define vidTime hisTime
# define vidZ hisZ
# define vidUb hisUb
# define vidVb hisVb
# define vidU hisU
# define vidV hisV
# define vidT hisT
# define vidR hisR
# define vidO hisO
# define vidW hisW
# define vidAkv hisAkv
# define vidAkt hisAkt
# define vidAks hisAks
# define vidHbls hisHbls
# define vidHbbl hisHbbl
# if defined BULK_FLUX
#  define vidSustr_blk hisSustr_blk
#  define vidSvstr_blk hisSvstr_blk
#  define vidShflx_net hisShflx_net
#  define vidShflx_lat hisShflx_lat
#  define vidShflx_sen hisShflx_sen
#  define vidShflx_rad hisShflx_rad
#  define vidSwflx_emp hisSwflx_emp
#  define vidShflx_lat hisShflx_lat
#  define vidShflx_wwk hisShflx_wwk
# endif
# if defined BIOLOGY_NPZDOC || defined BIOLOGY_BEC
#  define vidPH hisPH
#  define vidPCO2 hisPCO2
#  define vidPCO2air hisPCO2air
#  define vidPARinc hisPARinc
#  define vidPAR hisPAR
#  ifdef SEDIMENT_BIOLOGY
#   define vidTsed hisTsed
#  endif
# endif /* BIOLOGY_NPZDOC || BIOLOGY_BEC*/
# ifdef WKB_WWAVE
#  define vidWKB hisWKB
# endif
# ifdef MRL_WCI
#  define vidSUP hisSUP
#  define vidUST2D hisUST2D
#  define vidVST2D hisVST2D
#  define vidUST hisUST
#  define vidVST hisVST
#  define vidWST hisWST
#  define vidAkb hisAkb
#  define vidAkw hisAkw
#  define vidKVF hisKVF
#  define vidCALP hisCALP
#  define vidKAPS hisKAPS
# endif
#else
# define ncid ncidavg
# define ncname avgname
# define rec_per_file nrpfavg
# define wrt wrtavg
# define vidTime avgTime
# define vidZ avgZ
# define vidUb avgUb
# define vidVb avgVb
# define vidU avgU
# define vidV avgV
# define vidT avgT
# define vidR avgR
# define vidO avgO
# define vidW avgW
# define vidAkv avgAkv
# define vidAkt avgAkt
# define vidAks avgAks
# define vidHbls avgHbls
# define vidHbbl avgHbbl
# if defined BULK_FLUX
#  define vidSustr_blk avgSustr_blk
#  define vidSvstr_blk avgSvstr_blk
#  define vidShflx_net avgShflx_net
#  define vidShflx_lat avgShflx_lat
#  define vidShflx_sen avgShflx_sen
#  define vidShflx_rad avgShflx_rad
#  define vidSwflx_emp avgSwflx_emp
#  define vidShflx_wwk avgShflx_wwk
# endif
# if defined BIOLOGY_NPZDOC || defined BIOLOGY_BEC
#  define vidPH avgPH
#  define vidPCO2 avgPCO2
#  define vidPCO2air avgPCO2air
#  define vidPARinc avgPARinc
#  define vidPAR avgPAR
#  ifdef SEDIMENT_BIOLOGY
#    define vidTsed avgTsed
#  endif
#  ifdef OUTPUT_SPECDOM
#   define vidFreqDomSP_sfc avgFreqDomSP_sfc
#   define vidFreqDomDIAT_sfc avgFreqDomDIAT_sfc
#   define vidFreqDomDIAZ_sfc avgFreqDomDIAZ_sfc
#   define vidFreqDomSP_int avgFreqDomSP_int
#   define vidFreqDomDIAT_int avgFreqDomDIAT_int
#   define vidFreqDomDIAZ_int avgFreqDomDIAZ_int
#  endif /* OUTPUT_SPECDOM */
# endif /* BIOLOGY_NPZDOC || BIOLOGY_BEC */
# ifdef WKB_WWAVE
#  define vidWKB avgWKB
# endif
# ifdef MRL_WCI
#  define vidSUP avgSUP
#  define vidUST2D avgUST2D
#  define vidVST2D avgVST2D
#  define vidUST avgUST
#  define vidWST avgWST
#  define vidVST avgVST
#  define vidAkb avgAkb
#  define vidAkw avgAkw
#  define vidKVF avgKVF
#  define vidCALP avgCALP
#  define vidKAPS avgKAPS
# endif
#endif

! Put time record index into file name. In  the case when model
! output is to be arranged into sequence of named files, the naming
! convention is as follows: 'his_root.INDEX.[MPI_node.]nc', where
! INDEX is an integer number such that (i) it is divisible by the
! specified number of records per file; and (ii)
!
!      INDEX + record_within_the_file = total_record
!
! where, 1 =< record_within_the_file =< records_per_file, so that
! total_record changes continuously throughout the sequence of files.

      ierr=0
      lncn=lenstr(ncname)
      if (rec_per_file.gt.0) then
        lvar=total_rec-(1+mod(total_rec-1, rec_per_file))
        call insert_time_index (ncname, lncn, lvar, ierr)
        if (ierr .ne. 0) goto 99
      endif

! Decide whether to create a new file, or open an existing one. 
! Overall the code below is globally organized as a 3-way switch,
!
! 10  if (create_new_file) then
!        ... create new file, save netCDF ids for all variables;
!     elseif (ncid.eq.-1) then
!        ... try to open existing file and check its dimensions
!       if (cannot be opened or rejected) then
!         create_new_file=.true.
!         goto 10
!       endif   and prepare
!        .... prepare the file for adding new data,
!        .... find and save netCDF ids for all variables
!     else
!        .... just open, no checking, all ids are assumed to be
!        .... already known (MPI single file output only).
!     endif
!
! which is designed to implement a flexible opening policy: if 
! ldefhis=.true., it forces creation of a new file [hence, if the 
! file already exists, it will be overwritten]; conversely, if
! ldefhis=.false. it is assumed that the file already exists and an 
! attempt to open it is made; if the attempt is successful, the file 
! is prepared for appending hew data; if it fails, a new file will 
! be created.

      create_new_file=ldefhis            !<-- initialize
      if (ncid .ne. -1)  create_new_file=.false.
#if defined MPI & !defined PARALLEL_FILES
      if (mynode .gt. 0) create_new_file=.false.
#endif

! Create new history/averages file:    Put global attributes and
! ====== === ======= ======== =====    define all the variables.

  10  if (create_new_file) then
        ierr=nf_create (ncname(1:lncn), nf_clobber+nf_64bit_offset,
     &                                                        ncid)
        if (ierr .ne. nf_noerr) then
#ifdef AVRH
          write(*,'(/1x,4A/12x,A/)')  '### ERROR: def_avg :: ',
#else
          write(*,'(/1x,4A/12x,A/)')  '### ERROR: def_his :: ',
#endif
     &          'Cannot create netCDF file ''', ncname(1:lncn),
     &                              '''.',   nf_strerror(ierr)
          goto 99
        endif
        if (rec_per_file .eq. 0) total_rec=0

! Put global attributes.
! --- ------ -----------

        call put_global_atts (ncid, ierr)

! Define dimensions of staggered fields.
! ------ ---------- -- --------- -------

        ierr=my_nf_def_dim (ncid, 'xi_rho',  xi_rho,  r2dgrd(1))
        ierr=my_nf_def_dim (ncid, 'xi_u',    xi_u,    u2dgrd(1))
        ierr=my_nf_def_dim (ncid, 'eta_rho', eta_rho, r2dgrd(2))
        ierr=my_nf_def_dim (ncid, 'eta_v',   eta_v,   v2dgrd(2))
#ifdef SOLVE3D
        ierr=my_nf_def_dim (ncid, 's_rho',   N,       r3dgrd(3))
        ierr=my_nf_def_dim (ncid, 's_w',     N+1,     w3dgrd(3))
#endif
        ierr=my_nf_def_dim (ncid, 'time', nf_unlimited, timedim)
        ierr=my_nf_def_dim (ncid, 'auxil',   iaux,     auxil(1))

        auxil(2)=timedim
 
        r2dgrd(3)=timedim           ! Free surface
 
        u2dgrd(2)=r2dgrd(2)         ! 2D UBAR-type
        u2dgrd(3)=timedim
 
        v2dgrd(1)=r2dgrd(1)         ! 2D VBAR-type
        v2dgrd(3)=timedim
 
#ifdef SOLVE3D
        r3dgrd(1)=r2dgrd(1)         !
        r3dgrd(2)=r2dgrd(2)         ! 3D RHO-type
        r3dgrd(4)=timedim           !
 
        u3dgrd(1)=u2dgrd(1)         !
        u3dgrd(2)=r2dgrd(2)         ! 3D U-type
        u3dgrd(3)=r3dgrd(3)         !
        u3dgrd(4)=timedim
 
        v3dgrd(1)=r2dgrd(1)         !
        v3dgrd(2)=v2dgrd(2)         ! 3D V-type
        v3dgrd(3)=r3dgrd(3)         !
        v3dgrd(4)=timedim
 
        w3dgrd(1)=r2dgrd(1)         !
        w3dgrd(2)=r2dgrd(2)         ! 3D W-type
        w3dgrd(4)=timedim           !
#endif
#if (defined PUT_GRID_INTO_HISTORY && !defined AVRH)\
 || (defined PUT_GRID_INTO_AVERAGES && defined AVRH)

! Define grid variables.
! ------ ---- ----------

        if (total_rec.le.1) call def_grid (ncid, r2dgrd)
#endif

! Define evolving model variables.
! ------ -------- ----- ----------

! Time-step number and time-record indices: (history file only, this
! may be needed in the event when a history record is used to restart
! the current model run);
#ifndef AVRH
        ierr=nf_def_var (ncid, 'time_step', nf_int, 2, auxil, hisTstep)
        ierr=nf_put_att_text (ncid, hisTstep, 'long_name', 48,
     &              'time step and record numbers from initialization')
#endif

! Time.
        lvar=lenstr(vname(1,indxTime))
        ierr=nf_def_var (ncid, vname(1,indxTime)(1:lvar), NF_FTYPE,
     &                                         1, timedim, vidTime)
#ifdef AVRH
        text='averaged '/ /vname(2,indxTime)
#else
        text=vname(2,indxTime)
#endif
        lvar=lenstr(text)
        ierr=nf_put_att_text (ncid, vidTime, 'long_name', lvar,
     &                                             text(1:lvar))
        lvar=lenstr(vname(3,indxTime))
        ierr=nf_put_att_text (ncid, vidTime, 'units',  lvar,
     &                                vname(3,indxTime)(1:lvar))

! Free-surface and 2D barotropic velocity components

        if (wrt(indxZ)) then
          lvar=lenstr(vname(1,indxZ))
          ierr=nf_def_var (ncid, vname(1,indxZ)(1:lvar), NF_FOUT,
     &                                           3, r2dgrd, vidZ)
#ifdef AVRH
          text='averaged '/ /vname(2,indxZ)
#else
          text=vname(2,indxZ)
#endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidZ, 'long_name', lvar,
     &                                            text(1:lvar))
          lvar=lenstr(vname(3,indxZ))
          ierr=nf_put_att_text (ncid, vidZ, 'units',     lvar,
     &                                  vname(3,indxZ)(1:lvar))
#ifdef MASK_LAND_DATA
          ierr=nf_put_att_FOUT (ncid, vidZ, '_FillValue', NF_FOUT,
     &                                                 1, spv_set)
          if (ierr.ne.nf_noerr) then
             write(*,*)  'nf_put_att_XXX:', nf_strerror(ierr)
          endif    
#endif

        endif

        if (wrt(indxUb)) then
          lvar=lenstr(vname(1,indxUb))
          ierr=nf_def_var (ncid, vname(1,indxUb)(1:lvar), NF_FOUT,
     &                                           3, u2dgrd, vidUb)
#ifdef AVRH
          text='averaged '/ /vname(2,indxUb)
#else
          text=vname(2,indxUb)
#endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidUb, 'long_name', lvar,
     &                                             text(1:lvar))
          lvar=lenstr(vname(3,indxUb))
          ierr=nf_put_att_text (ncid, vidUb, 'units',     lvar,
     &                                  vname(3,indxUb)(1:lvar))
#ifdef MASK_LAND_DATA
          ierr=nf_put_att_FOUT (ncid, vidUb, '_FillValue', NF_FOUT,
     &                                                  1, spv_set)
#endif
        endif
 
        if (wrt(indxVb)) then
          lvar=lenstr(vname(1,indxVb))
          ierr=nf_def_var (ncid, vname(1,indxVb)(1:lvar), NF_FOUT,
     &                                           3, v2dgrd, vidVb)
#ifdef AVRH
          text='averaged '/ /vname(2,indxVb)
#else
          text=vname(2,indxVb)
#endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidVb, 'long_name', lvar,
     &                                             text(1:lvar))
          lvar=lenstr(vname(3,indxVb))
          ierr=nf_put_att_text (ncid, vidVb, 'units',     lvar,
     &                                  vname(3,indxVb)(1:lvar))
#ifdef MASK_LAND_DATA
          ierr=nf_put_att_FOUT (ncid, vidVb, '_FillValue', NF_FOUT,
     &                                                  1, spv_set)
#endif
        endif
#ifdef SOLVE3D

! 3D momenta in XI- and ETA-directions.

        if (wrt(indxU)) then
          lvar=lenstr(vname(1,indxU))
          ierr=nf_def_var (ncid, vname(1,indxU)(1:lvar), NF_FOUT,
     &                                           4, u3dgrd, vidU)
# ifdef AVRH
          text='averaged '/ /vname(2,indxU)
# else
          text=vname(2,indxU)
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidU, 'long_name', lvar,
     &                                            text(1:lvar))
          lvar=lenstr(vname(3,indxU))
          ierr=nf_put_att_text (ncid, vidU, 'units',     lvar,
     &                                  vname(3,indxU)(1:lvar))
# ifdef MASK_LAND_DATA
          ierr=nf_put_att_FOUT (ncid, vidU, '_FillValue', NF_FOUT,
     &                                                 1, spv_set)
# endif

        endif
 
        if (wrt(indxV)) then
          lvar=lenstr(vname(1,indxV))
          ierr=nf_def_var (ncid, vname(1,indxV)(1:lvar), NF_FOUT,
     &                                           4, v3dgrd, vidV)
# ifdef AVRH
          text='averaged '/ /vname(2,indxV)
# else
          text=vname(2,indxV)
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidV, 'long_name', lvar,
     &                                            text(1:lvar))
          lvar=lenstr(vname(3,indxV))
          ierr=nf_put_att_text (ncid, vidV, 'units',     lvar,
     &                                  vname(3,indxV)(1:lvar))
# ifdef MASK_LAND_DATA
          ierr=nf_put_att_FOUT (ncid, vidV, '_FillValue', NF_FOUT,
     &                                                 1, spv_set)
# endif
        endif

! Tracer variables.

# if defined OUTPUT_AVG_PHYS_BIO && defined AVRH
        do itrc=1,1+ntrc_salt
# else
        do itrc=1,NT
# endif
          if (wrt(indxT+itrc-1)) then
            lvar=lenstr(vname(1,indxT+itrc-1))
            ierr=nf_def_var (ncid, vname(1,indxT+itrc-1)(1:lvar),
     &                             NF_FOUT, 4, r3dgrd, vidT(itrc))
# ifdef AVRH
            text='averaged '/ /vname(2,indxT+itrc-1)
# else
            text=vname(2,indxT+itrc-1)
# endif
            lvar=lenstr(text)
            ierr=nf_put_att_text (ncid, vidT(itrc), 'long_name',
     &                                          lvar, text(1:lvar))
            lvar=lenstr(vname(3,indxT+itrc-1))
            ierr=nf_put_att_text (ncid, vidT(itrc), 'units', lvar,
     &                               vname(3,indxT+itrc-1)(1:lvar))
# ifdef MASK_LAND_DATA
            ierr=nf_put_att_FOUT (ncid, vidT(itrc), '_FillValue',
     &                                      NF_FOUT, 1, spv_set)
# endif
          endif
        enddo

#if defined SEDIMENT_BIOLOGY && (!defined OUTPUT_AVG_PHYS_BIO || !defined AVRH)
!
! Sediment tracer variables.
!
        do itrc=1,NT_sed
          if (wrt(indxSedFirst+itrc-1)) then
            lvar=lenstr(vname(1,indxSedFirst+itrc-1))
            ierr=nf_def_var (ncid, 
     &           vname(1,indxSedFirst+itrc-1)(1:lvar),
     &           NF_FOUT, 3, r2dgrd, vidTsed(itrc))
# ifndef AVRH
            lvar=lenstr(vname(2,indxSedFirst+itrc-1))
            ierr=nf_put_att_text (ncid, vidTsed(itrc), 'long_name',
     &           lvar, vname(2,indxSedFirst+itrc-1)(1:lvar))
# else
            text='averaged '/ /vname(2,indxSedFirst+itrc-1)
            lvar=lenstr(text)
            ierr=nf_put_att_text (ncid, vidTsed(itrc), 'long_name',
     &           lvar, text(1:lvar))
# endif
            lvar=lenstr(vname(3,indxSedFirst+itrc-1))
            ierr=nf_put_att_text (ncid, vidTsed(itrc), 'units', lvar,
     &           vname(3,indxSedFirst+itrc-1)(1:lvar))
# ifdef MASK_LAND_DATA
            ierr=nf_put_att_FOUT (ncid, vidTsed(itrc), '_FillValue',
     &                                      NF_FOUT, 1, spv_set)
# endif
          endif
        enddo
# endif /* SEDIMENT_BIOLOGY...  */

! Density anomaly.

        if (wrt(indxR)) then
          lvar=lenstr(vname(1,indxR))
          ierr=nf_def_var (ncid, vname(1,indxR)(1:lvar), NF_FOUT,
     &                                           4, r3dgrd, vidR)
# ifdef AVRH
          text='averaged '/ /vname(2,indxR)
# else
          text=vname(2,indxR)
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidR, 'long_name', lvar,
     &                                            text(1:lvar))
          lvar=lenstr(vname(3,indxR))
          ierr=nf_put_att_text (ncid, vidR, 'units',     lvar,
     &                                  vname(3,indxR)(1:lvar))
# ifdef MASK_LAND_DATA
          ierr=nf_put_att_FOUT (ncid, vidR, '_FillValue', NF_FOUT,
     &                                                 1, spv_set)
# endif
        endif

! S-coordinate "omega" vertical velocity.

        if (wrt(indxO)) then
          lvar=lenstr(vname(1,indxO))
          ierr=nf_def_var (ncid, vname(1,indxO)(1:lvar), NF_FOUT,
     &                                           4, w3dgrd, vidO)
# ifdef AVRH
          text='averaged '/ /vname(2,indxO)
# else
          text=vname(2,indxO)
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidO, 'long_name', lvar,
     &                                            text(1:lvar))
          lvar=lenstr(vname(3,indxO))
          ierr=nf_put_att_text (ncid, vidO, 'units',     lvar,
     &                                  vname(3,indxO)(1:lvar))
# ifdef MASK_LAND_DATA 
          ierr=nf_put_att_FOUT (ncid, vidO, '_FillValue', NF_FOUT,
     &                                                 1, spv_set)
# endif
        endif

! True W-vertical velocity.

        if (wrt(indxW)) then
          lvar=lenstr(vname(1,indxW))
          ierr=nf_def_var (ncid, vname(1,indxW)(1:lvar), NF_FOUT,
     &                                           4, r3dgrd, vidW)
# ifdef AVRH
          text='averaged '/ /vname(2,indxW)
# else
          text=vname(2,indxW)
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidW, 'long_name', lvar,
     &                                            text(1:lvar))
          lvar=lenstr(vname(3,indxW))
          ierr=nf_put_att_text (ncid, vidW, 'units',     lvar,
     &                                  vname(3,indxW)(1:lvar))
# ifdef MASK_LAND_DATA 
          ierr=nf_put_att_FOUT (ncid, vidW, '_FillValue', NF_FOUT,
     &                                                 1, spv_set)
# endif
        endif

! Vertical viscosity coefficient.

        if (wrt(indxAkv)) then
          lvar=lenstr(vname(1,indxAkv))
          ierr=nf_def_var (ncid, vname(1,indxAkv)(1:lvar), NF_FOUT,
     &                                           4, w3dgrd, vidAkv)
# ifdef AVRH
          text='averaged '/ /vname(2,indxAkv)
# else
          text=vname(2,indxAkv)
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidAkv, 'long_name', lvar,
     &                                              text(1:lvar))
          lvar=lenstr(vname(3,indxAkv))
          ierr=nf_put_att_text (ncid, vidAkv, 'units',     lvar,
     &                                  vname(3,indxAkv)(1:lvar))
# ifdef MASK_LAND_DATA 
          ierr=nf_put_att_FOUT (ncid, vidAkv, '_FillValue', NF_FOUT,
     &                                                   1, spv_set)
# endif
        endif

! Vertical thermal conductivity coefficient.

        if (wrt(indxAkt)) then
          lvar=lenstr(vname(1,indxAkt))
          ierr=nf_def_var (ncid, vname(1,indxAkt)(1:lvar), NF_FOUT,
     &                                           4, w3dgrd, vidAkt)
# ifdef AVRH
          text='averaged '/ /vname(2,indxAkt)
# else
          text=vname(2,indxAkt)
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidAkt, 'long_name', lvar,
     &                                              text(1:lvar))
          lvar=lenstr(vname(3,indxAkt))
          ierr=nf_put_att_text (ncid, vidAkt, 'units',     lvar,
     &                                  vname(3,indxAkt)(1:lvar))
# ifdef MASK_LAND_DATA 
          ierr=nf_put_att_FOUT (ncid, vidAkt, '_FillValue', NF_FOUT,
     &                                                   1, spv_set)
# endif
        endif
# ifdef SALINITY
! Vertical diffusion coefficient for salinity.

        if (wrt(indxAks)) then
          lvar=lenstr(vname(1,indxAks))
          ierr=nf_def_var (ncid, vname(1,indxAks)(1:lvar), NF_FOUT,
     &                                           4, w3dgrd, vidAks)
#  ifdef AVRH
          text='averaged '/ /vname(2,indxAks)
#  else
          text=vname(2,indxAks)
#  endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidAks, 'long_name', lvar,
     &                                              text(1:lvar))
          lvar=lenstr(vname(3,indxAks))
          ierr=nf_put_att_text (ncid, vidAks, 'units',     lvar,
     &                                  vname(3,indxAks)(1:lvar))
#  ifdef MASK_LAND_DATA 
          ierr=nf_put_att_FOUT (ncid, vidAks, '_FillValue', NF_FOUT,
     &                                                   1, spv_set)
#  endif
        endif
# endif /* SALINITY */
# ifdef LMD_KPP
! Depth of surface planetary boundary layer.

        if (wrt(indxHbls)) then
          lvar=lenstr(vname(1,indxHbls))
          ierr=nf_def_var (ncid, vname(1,indxHbls)(1:lvar),
     &                           NF_FOUT, 3, r2dgrd, vidHbls)
#  ifdef AVRH
          text='averaged '/ /vname(2,indxHbls)
#  else
          text=vname(2,indxHbls)
#  endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidHbls, 'long_name', lvar,
     &                                              text(1:lvar))
          lvar=lenstr(vname(3,indxHbls))
          ierr=nf_put_att_text (ncid, vidHbls, 'units',     lvar,
     &                                 vname(3,indxHbls)(1:lvar))
#  ifdef MASK_LAND_DATA 
          ierr=nf_put_att_FOUT (ncid, vidHbls, '_FillValue',
     &                                       NF_FOUT, 1, spv_set)
#  endif
        endif
# endif
# ifdef LMD_BKPP
! Thickness of bottom boundary layer.

        if (wrt(indxHbbl)) then
          lvar=lenstr(vname(1,indxHbbl))
          ierr=nf_def_var (ncid, vname(1,indxHbbl)(1:lvar),
     &                            NF_FOUT, 3, r2dgrd, vidHbbl)
#  ifdef AVRH
          text='averaged '/ /vname(2,indxHbbl)
#  else
          text=vname(2,indxHbbl)
#  endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidHbbl, 'long_name', lvar,
     &                                              text(1:lvar))
          lvar=lenstr(vname(3,indxHbbl))
          ierr=nf_put_att_text (ncid, vidHbbl, 'units',     lvar,
     &                                 vname(3,indxHbbl)(1:lvar))
#  ifdef MASK_LAND_DATA
          ierr=nf_put_att_FOUT (ncid, vidHbbl, '_FillValue',
     &                                       NF_FOUT, 1, spv_set)
#  endif
        endif
# endif

# if (defined BIOLOGY_BEC || (defined BIOLOGY_NPZDOC && \
     defined OXYGEN && defined CARBON) ) && \
         (! defined OUTPUT_AVG_PHYS_BIO || ! defined AVRH)
!
! pH
!
        lvar=lenstr(vname(1,indxPH_rst))
        ierr=nf_def_var (ncid, vname(1,indxPH_rst)(1:lvar), NF_FOUT,
     &                                           3, r2dgrd, vidPH)
#  ifdef AVRH
        text='averaged '/ /vname(2,indxPH_rst)
#  else
        text=vname(2,indxPH_rst)
#  endif
        lvar=lenstr(text)
        ierr=nf_put_att_text (ncid, vidPH, 'long_name', lvar,
     &                                              text(1:lvar))
        lvar=lenstr(vname(3,indxPH_rst))
        ierr=nf_put_att_text (ncid, vidPH, 'units',     lvar,
     &                                  vname(3,indxPH_rst)(1:lvar))
#ifdef MASK_LAND_DATA 
        ierr=nf_put_att_FOUT (ncid, vidPH, '_FillValue',
     &                                  NF_FOUT, 1, spv_set)
#endif
!
! pCO2:
!
        lvar=lenstr(vname(1,indxPCO2_rst))
        ierr=nf_def_var (ncid, vname(1,indxPCO2_rst)(1:lvar), NF_FOUT,
     &                                           3, r2dgrd, vidPCO2)
#  ifdef AVRH
        text='averaged '/ /vname(2,indxPCO2_rst)
#  else
        text=vname(2,indxPCO2_rst)
#  endif
        lvar=lenstr(text)
        ierr=nf_put_att_text (ncid, vidPCO2, 'long_name', lvar,
     &                                              text(1:lvar))
        lvar=lenstr(vname(3,indxPCO2_rst))
        ierr=nf_put_att_text (ncid, vidPCO2, 'units',     lvar,
     &                                  vname(3,indxPCO2_rst)(1:lvar))
#ifdef MASK_LAND_DATA 
        ierr=nf_put_att_FOUT (ncid, vidPCO2, '_FillValue',
     &                                  NF_FOUT, 1, spv_set)
#endif
!
! pCO2air:
!
        lvar=lenstr(vname(1,indxPCO2air_rst))
        ierr=nf_def_var (ncid, vname(1,indxPCO2air_rst)(1:lvar), NF_FOUT,
     &                                           3, r2dgrd, vidPCO2air)
#  ifdef AVRH
        text='averaged '/ /vname(2,indxPCO2air_rst)
#  else
        text=vname(2,indxPCO2air_rst)
#  endif
        lvar=lenstr(text)
        ierr=nf_put_att_text (ncid, vidPCO2air, 'long_name', lvar,
     &                                              text(1:lvar))
        lvar=lenstr(vname(3,indxPCO2air_rst))
        ierr=nf_put_att_text (ncid, vidPCO2air, 'units',     lvar,
     &       vname(3,indxPCO2air_rst)(1:lvar))
#ifdef MASK_LAND_DATA 
        ierr=nf_put_att_FOUT (ncid, vidPCO2air, '_FillValue',
     &                                  NF_FOUT, 1, spv_set)
#endif

# endif /* BIOLOGY_BEC || (BIOLOGY_NPZDOC && OXYGEN && CARBON) ... */

# if (defined BIOLOGY_BEC || defined BIOLOGY_NPZDOC) && (! defined OUTPUT_AVG_PHYS_BIO || ! defined AVRH)
!
! PARinc:
!
        lvar=lenstr(vname(1,indxPARinc_rst))
        ierr=nf_def_var (ncid, vname(1,indxPARinc_rst)(1:lvar), NF_FOUT,
     &       3, r2dgrd, vidPARinc)
#  ifdef AVRH
        text='averaged '/ /vname(2,indxPARinc_rst)
#  else
        text=vname(2,indxPARinc_rst)
#  endif
        lvar=lenstr(text)
        ierr=nf_put_att_text (ncid, vidPARinc, 'long_name', lvar,
     &                                              text(1:lvar))
        lvar=lenstr(vname(3,indxPARinc_rst))
        ierr=nf_put_att_text (ncid, vidPARinc, 'units',     lvar,
     &                                  vname(3,indxPARinc_rst)(1:lvar))
#ifdef MASK_LAND_DATA 
        ierr=nf_put_att_FOUT (ncid, vidPARinc, '_FillValue',
     &                                  NF_FOUT, 1, spv_set)
#endif
!
! PAR:
!
        lvar=lenstr(vname(1,indxPAR_rst))
        ierr=nf_def_var (ncid, vname(1,indxPAR_rst)(1:lvar), NF_FOUT,
     &                                           4, r3dgrd, vidPAR)
#  ifdef AVRH
        text='averaged '/ /vname(2,indxPAR_rst)
#  else
        text=vname(2,indxPAR_rst)
#  endif
        lvar=lenstr(text)
        ierr=nf_put_att_text (ncid, vidPAR, 'long_name', lvar,
     &                                              text(1:lvar))
        lvar=lenstr(vname(3,indxPAR_rst))
        ierr=nf_put_att_text (ncid, vidPAR, 'units',     lvar,
     &                                  vname(3,indxPAR_rst)(1:lvar))
#ifdef MASK_LAND_DATA 
        ierr=nf_put_att_FOUT (ncid, vidPAR, '_FillValue',
     &                                  NF_FOUT, 1, spv_set)
#endif
# if defined OUTPUT_SPECDOM && defined AVRH && !defined OUTPUT_AVG_PHYS_BIO
! frequency of small phytoplankton dominance at the surface
        lvar=lenstr(vname(1,indxFreqDomSP_sfc))
        ierr=nf_def_var (ncid, vname(1,indxFreqDomSP_sfc)(1:lvar), 
     &       NF_FOUT, 3, r2dgrd, vidFreqDomSP_sfc)
#  ifdef AVRH
        text='averaged '/ /vname(2,indxFreqDomSP_sfc)
#  else
        text=vname(2,indxFreqDomSP_sfc)
#  endif
        lvar=lenstr(text)
        ierr=nf_put_att_text (ncid, vidFreqDomSP_sfc, 'long_name', lvar,
     &                                              text(1:lvar))
        lvar=lenstr(vname(3,indxFreqDomSP_sfc))
        ierr=nf_put_att_text (ncid, vidFreqDomSP_sfc, 'units',     lvar,
     &       vname(3,indxFreqDomSP_sfc)(1:lvar))
#ifdef MASK_LAND_DATA 
        ierr=nf_put_att_FOUT (ncid, vidFreqDomSP_sfc, '_FillValue',
     &                                  NF_FOUT, 1, spv_set)
#endif
! frequency of diatom dominance at the surface
        lvar=lenstr(vname(1,indxFreqDomDIAT_sfc))
        ierr=nf_def_var (ncid, vname(1,indxFreqDomDIAT_sfc)(1:lvar), 
     &       NF_FOUT, 3, r2dgrd, vidFreqDomDIAT_sfc)
#  ifdef AVRH
        text='averaged '/ /vname(2,indxFreqDomDIAT_sfc)
#  else
        text=vname(2,indxFreqDomDIAT_sfc)
#  endif
        lvar=lenstr(text)
        ierr=nf_put_att_text (ncid, vidFreqDomDIAT_sfc, 'long_name',
     &       lvar, text(1:lvar))
        lvar=lenstr(vname(3,indxFreqDomDIAT_sfc))
        ierr=nf_put_att_text (ncid, vidFreqDomDIAT_sfc, 'units', lvar,
     &       vname(3,indxFreqDomDIAT_sfc)(1:lvar))
#ifdef MASK_LAND_DATA 
        ierr=nf_put_att_FOUT (ncid, vidFreqDomDIAT_sfc, '_FillValue',
     &                                  NF_FOUT, 1, spv_set)
#endif
! frequency of diazotroph dominance at the surface
        lvar=lenstr(vname(1,indxFreqDomDIAZ_sfc))
        ierr=nf_def_var (ncid, vname(1,indxFreqDomDIAZ_sfc)(1:lvar), 
     &       NF_FOUT, 3, r2dgrd, vidFreqDomDIAZ_sfc)
#  ifdef AVRH
        text='averaged '/ /vname(2,indxFreqDomDIAZ_sfc)
#  else
        text=vname(2,indxFreqDomDIAZ_sfc)
#  endif
        lvar=lenstr(text)
        ierr=nf_put_att_text (ncid, vidFreqDomDIAZ_sfc, 'long_name',
     &       lvar, text(1:lvar))
        lvar=lenstr(vname(3,indxFreqDomDIAZ_sfc))
        ierr=nf_put_att_text (ncid, vidFreqDomDIAZ_sfc, 'units', lvar,
     &       vname(3,indxFreqDomDIAZ_sfc)(1:lvar))
#ifdef MASK_LAND_DATA 
        ierr=nf_put_att_FOUT (ncid, vidFreqDomDIAZ_sfc, '_FillValue',
     &                                  NF_FOUT, 1, spv_set)
#endif
! frequency of small phytoplankton dominance, integrated
        lvar=lenstr(vname(1,indxFreqDomSP_int))
        ierr=nf_def_var (ncid, vname(1,indxFreqDomSP_int)(1:lvar), 
     &       NF_FOUT, 3, r2dgrd, vidFreqDomSP_int)
#  ifdef AVRH
        text='averaged '/ /vname(2,indxFreqDomSP_int)
#  else
        text=vname(2,indxFreqDomSP_int)
#  endif
        lvar=lenstr(text)
        ierr=nf_put_att_text (ncid, vidFreqDomSP_int, 'long_name', lvar,
     &                                              text(1:lvar))
        lvar=lenstr(vname(3,indxFreqDomSP_int))
        ierr=nf_put_att_text (ncid, vidFreqDomSP_int, 'units',     lvar,
     &       vname(3,indxFreqDomSP_int)(1:lvar))
#ifdef MASK_LAND_DATA 
        ierr=nf_put_att_FOUT (ncid, vidFreqDomSP_int, '_FillValue',
     &                                  NF_FOUT, 1, spv_set)
#endif
! frequency of diatom dominance, integrated
        lvar=lenstr(vname(1,indxFreqDomDIAT_int))
        ierr=nf_def_var (ncid, vname(1,indxFreqDomDIAT_int)(1:lvar), 
     &       NF_FOUT, 3, r2dgrd, vidFreqDomDIAT_int)
#  ifdef AVRH
        text='averaged '/ /vname(2,indxFreqDomDIAT_int)
#  else
        text=vname(2,indxFreqDomDIAT_int)
#  endif
        lvar=lenstr(text)
        ierr=nf_put_att_text (ncid, vidFreqDomDIAT_int, 'long_name',
     &       lvar, text(1:lvar))
        lvar=lenstr(vname(3,indxFreqDomDIAT_int))
        ierr=nf_put_att_text (ncid, vidFreqDomDIAT_int, 'units', lvar,
     &       vname(3,indxFreqDomDIAT_int)(1:lvar))
#ifdef MASK_LAND_DATA 
        ierr=nf_put_att_FOUT (ncid, vidFreqDomDIAT_int, '_FillValue',
     &                                  NF_FOUT, 1, spv_set)
#endif
! frequency of diazotroph dominance, integrated
        lvar=lenstr(vname(1,indxFreqDomDIAZ_int))
        ierr=nf_def_var (ncid, vname(1,indxFreqDomDIAZ_int)(1:lvar), 
     &       NF_FOUT, 3, r2dgrd, vidFreqDomDIAZ_int)
#  ifdef AVRH
        text='averaged '/ /vname(2,indxFreqDomDIAZ_int)
#  else
        text=vname(2,indxFreqDomDIAZ_int)
#  endif
        lvar=lenstr(text)
        ierr=nf_put_att_text (ncid, vidFreqDomDIAZ_int, 'long_name',
     &       lvar, text(1:lvar))
        lvar=lenstr(vname(3,indxFreqDomDIAZ_int))
        ierr=nf_put_att_text (ncid, vidFreqDomDIAZ_int, 'units', lvar,
     &       vname(3,indxFreqDomDIAZ_int)(1:lvar))
#ifdef MASK_LAND_DATA 
        ierr=nf_put_att_FOUT (ncid, vidFreqDomDIAZ_int, '_FillValue',
     &                                  NF_FOUT, 1, spv_set)
#endif
# endif /* OUTPUT_SPECDOM */
#endif /* BIOLOGY_BEC || BIOLOGY_NPZDOC ... */
#endif /* SOLVE3D */

! Leave definition mode.
! ----- ---------- -----

        ierr=nf_enddef(ncid)
#ifdef MPI_SILENT_MODE
        if (mynode.eq.0) then
#endif
#ifdef AVRH
         write(*,'(7x,4A,I4)') 'def_avg :: created new netCDF file ''',
#else
         write(*,'(7x,4A,I4)') 'def_his :: created new netCDF file ''',
#endif
     &                          ncname(1:lncn), '''.' MYID
#ifdef MPI_SILENT_MODE
        endif
#endif

! Open existing file and prepare for appending data:   Inquire about
! ==== ======== ==== === ======= === ========= ====    dimensions and
! variables. Check for consistency with model dimensions. In the case
! when file is rejected (whether it cannot be opened, or something is
! wrong with its dimensions) create a new file.  After that verify
! that all necessary variables are already defined, find and save
! their netCDF IDs.

      elseif (ncid.eq.-1) then  !--> alt. to "create_new_file"
        ierr=nf_open (ncname(1:lncn), nf_write, ncid)
        if (ierr. eq. nf_noerr) then
#ifdef MPI_SILENT_MODE
          if (mynode.eq.0) then 
#endif
#ifdef AVRH
            write(*,'(7x,5A,I4)')    'def_avg :: opened existing ',
#else
            write(*,'(7x,5A,I4)')    'def_his :: opened existing ',
#endif
     &                     'file ''', ncname(1:lncn),  '''.' MYID
#ifdef MPI_SILENT_MODE
          endif
#endif
          ierr=checkdims (ncid, ncname, rec)
          if (ierr .eq. nf_noerr) then
            if (rec_per_file.eq.0) then
              ierr=rec+1 - total_rec
            else
              ierr=rec+1 - (1+mod(total_rec-1, rec_per_file))
            endif
            if (ierr.gt.0) then
              mpi_master_only write(*,
     &                '(/1x,A,I5,1x,3A/21x,2(A,I5),1x,A/21x,A/)')
#ifdef AVRH
     &            'WARNING: def_avg :: The actual number of records',
#else
     &            'WARNING: def_his :: The actual number of records',
#endif
     &             rec, 'present in file ''',  ncname(1:lncn),  '''',
     &            'exceeds record', rec+1-ierr, '/', total_rec,
     &            'specified by restart initial conditions.',
     &        'All records beyond this number will be overwritten.'
              rec=rec-ierr
            elseif (rec_per_file.eq.0) then
              total_rec=rec+1           ! <-- set to the next record
#if defined MPI & !defined PARALLEL_FILES
              if (mynode.gt.0) total_rec=total_rec-1
#endif
            endif
            ierr=nf_noerr
          endif
        endif
 
        if (ierr. ne. nf_noerr) then
#if defined MPI & !defined PARALLEL_FILES
          if (mynode.eq.0) then
            create_new_file=.true.
            goto 10
          else
# ifdef AVRH
            write(*,'(/1x,5A,I4/)')    '### ERROR: def_avg :: ',
# else
            write(*,'(/1x,5A,I4/)')    '### ERROR: def_his :: ',
# endif
     &        'Cannot open file ''', ncname(1:lncn), '''.' MYID
            goto 99                                       !--> ERROR
          endif
#else
          create_new_file=.true.
          goto 10
#endif
        endif

! Find netCDF IDs of evolving model variables:
! ---- ------ --- -- -------- ----- ----------

! Time-step indices, ocean time:
#ifndef AVRH
        ierr=nf_inq_varid (ncid, 'time_step', hisTstep)
        if (ierr .ne. nf_noerr) then
          write(*,1) 'time_step', ncname(1:lncn)
          goto 99                                         !--> ERROR
        endif
#endif

        lvar=lenstr(vname(1,indxTime))
        ierr=nf_inq_varid (ncid,vname(1,indxTime)(1:lvar),vidTime)
        if (ierr .ne. nf_noerr) then
          write(*,1) vname(1,indxTime)(1:lvar), ncname(1:lncn)
          goto 99                                         !--> ERROR
        endif

! Free-surface and 2D barotropic velocity components

        if (wrt(indxZ)) then
          lvar=lenstr(vname(1,indxZ))
          ierr=nf_inq_varid (ncid, vname(1,indxZ)(1:lvar), vidZ)
          if (ierr .ne. nf_noerr) then
            write(*,1) vname(1,indxZ)(1:lvar), ncname(1:lncn)
            goto 99                                       !--> ERROR
          endif
        endif

        if (wrt(indxUb)) then
          lvar=lenstr(vname(1,indxUb))
          ierr=nf_inq_varid (ncid, vname(1,indxUb)(1:lvar), vidUb)
          if (ierr .ne. nf_noerr) then
            write(*,1) vname(1,indxUb)(1:lvar), ncname(1:lncn)
            goto 99                                       !--> ERROR
          endif
        endif
 
        if (wrt(indxVb)) then
          lvar=lenstr(vname(1,indxVb))
          ierr=nf_inq_varid (ncid, vname(1,indxVb)(1:lvar), vidVb)
          if (ierr .ne. nf_noerr) then
            write(*,1) vname(1,indxVb)(1:lvar), ncname(1:lncn)
            goto 99                                       !--> ERROR
          endif
        endif
 
#ifdef SOLVE3D

! 3D momenta in XI- and ETA-directions.

        if (wrt(indxU)) then
          lvar=lenstr(vname(1,indxU))
          ierr=nf_inq_varid (ncid, vname(1,indxU)(1:lvar), vidU)
          if (ierr .ne. nf_noerr) then
            write(*,1) vname(1,indxU)(1:lvar), ncname(1:lncn)
            goto 99                                       !--> ERROR
          endif
        endif
        if (wrt(indxV)) then
          lvar=lenstr(vname(1,indxV))
          ierr=nf_inq_varid (ncid, vname(1,indxV)(1:lvar), vidV)
          if (ierr .ne. nf_noerr) then
            write(*,1) vname(1,indxV)(1:lvar), ncname(1:lncn)
            goto 99                                       !--> ERROR
          endif
        endif

! Tracer variables.

# if defined OUTPUT_AVG_PHYS_BIO && defined AVRH
        do itrc=1,1+ntrc_salt
# else
        do itrc=1,NT
# endif
          if (wrt(indxT+itrc-1)) then
            lvar=lenstr(vname(1,indxT+itrc-1))
            ierr=nf_inq_varid (ncid, vname(1,indxT+itrc-1)(1:lvar),
     &                                                  vidT(itrc))
            if (ierr .ne. nf_noerr) then
              write(*,1) vname(1,indxT+itrc-1)(1:lvar), ncname(1:lncn)
              goto 99                                     !--> ERROR
            endif
          endif
        enddo

# if defined SEDIMENT_BIOLOGY && (!defined OUTPUT_AVG_PHYS_BIO || !defined AVRH)
!
! Sediment tracer variables.
!
        do itrc=1,NT_sed
          if (wrt(indxSedFirst+itrc-1)) then
            lvar=lenstr(vname(1,indxSedFirst+itrc-1))
            ierr=nf_inq_varid (ncid, 
     &           vname(1,indxSedFirst+itrc-1)(1:lvar),vidTsed(itrc))
            if (ierr .ne. nf_noerr) then
              write(*,1) vname(1,indxSedFirst+itrc-1)(1:lvar),
     &                                       ncname(1:lncn) 
              goto 99                                     !--> ERROR
            endif
          endif
        enddo
# endif /* SEDIMENT_BIOLOGY  */

! Density anomaly.

        if (wrt(indxR)) then
          lvar=lenstr(vname(1,indxR))
          ierr=nf_inq_varid (ncid, vname(1,indxR)(1:lvar), vidR)
          if (ierr .ne. nf_noerr) then
            write(*,1) vname(1,indxR)(1:lvar), ncname(1:lncn)
            goto 99                                       !--> ERROR
          endif
        endif

! S-coordinate "omega" and "true" W- vertical velocities.

        if (wrt(indxO)) then
          lvar=lenstr(vname(1,indxO))
          ierr=nf_inq_varid (ncid, vname(1,indxO)(1:lvar), vidO)
          if (ierr .ne. nf_noerr) then
            write(*,1) vname(1,indxO)(1:lvar), ncname(1:lncn)
            goto 99                                       !--> ERROR
          endif
        endif

        if (wrt(indxW)) then
          lvar=lenstr(vname(1,indxW))
          ierr=nf_inq_varid (ncid, vname(1,indxW)(1:lvar), vidW)
          if (ierr .ne. nf_noerr) then
            write(*,1) vname(1,indxW)(1:lvar), ncname(1:lncn)
            goto 99                                       !--> ERROR
          endif
        endif

! Vertical viscosity, thermal conductivity, and salinity diffusion.

        if (wrt(indxAkv)) then
          lvar=lenstr(vname(1,indxAkv))
          ierr=nf_inq_varid (ncid, vname(1,indxAkv)(1:lvar), vidAkv)
          if (ierr .ne. nf_noerr) then
            write(*,1) vname(1,indxAkv)(1:lvar), ncname(1:lncn)
            goto 99                                       !--> ERROR
          endif
        endif

        if (wrt(indxAkt)) then
          lvar=lenstr(vname(1,indxAkt))
          ierr=nf_inq_varid (ncid,vname(1,indxAkt)(1:lvar), vidAkt)
          if (ierr .ne. nf_noerr) then
            write(*,1) vname(1,indxAkt)(1:lvar), ncname(1:lncn)
            goto 99                                       !--> ERROR
          endif
        endif
 
# ifdef SALINITY
        if (wrt(indxAks)) then
          lvar=lenstr(vname(1,indxAks))
          ierr=nf_inq_varid (ncid,vname(1,indxAks)(1:lvar), vidAks)
          if (ierr .ne. nf_noerr) then
            write(*,1) vname(1,indxAks)(1:lvar), ncname(1:lncn)
            goto 99                                       !--> ERROR
          endif
        endif
# endif
# ifdef LMD_KPP
! Depth of surface planetary boundary layer.

        if (wrt(indxHbls)) then
          lvar=lenstr(vname(1,indxHbls))
          ierr=nf_inq_varid (ncid,vname(1,indxHbls)(1:lvar), vidHbls)
          if (ierr .ne. nf_noerr) then
          write(*,1) vname(1,indxHbls)(1:lvar), ncname(1:lncn)
            goto 99                                       !--> ERROR
          endif
        endif
# endif
# ifdef LMD_BKPP
! Thickness of bottom boundary layer.

        if (wrt(indxHbbl)) then
          lvar=lenstr(vname(1,indxHbbl))
          ierr=nf_inq_varid (ncid,vname(1,indxHbbl)(1:lvar), vidHbbl)
          if (ierr .ne. nf_noerr) then
          write(*,1) vname(1,indxHbbl)(1:lvar), ncname(1:lncn)
            goto 99                                       !--> ERROR
          endif
        endif
# endif

# if (defined BIOLOGY_BEC || (defined BIOLOGY_NPZDOC && \
     defined OXYGEN && defined CARBON) ) && \
         (! defined OUTPUT_AVG_PHYS_BIO || ! defined AVRH)
! pH
        lvar=lenstr(vname(1,indxPH_rst))
        ierr=nf_inq_varid (ncid, vname(1,indxPH_rst)(1:lvar), vidPH)
        if (ierr .ne. nf_noerr) then
           write(*,1) vname(1,indxPH_rst)(1:lvar), ncname(1:lncn)
           goto 99              !--> ERROR
        endif
! pCO2:
        lvar=lenstr(vname(1,indxPCO2_rst))
        ierr=nf_inq_varid (ncid, vname(1,indxPCO2_rst)(1:lvar), 
     &       vidPCO2)
        if (ierr .ne. nf_noerr) then
           write(*,1) vname(1,indxPCO2_rst)(1:lvar), ncname(1:lncn)
           goto 99              !--> ERROR
        endif
! pCO2air:
        lvar=lenstr(vname(1,indxPCO2air_rst))
        ierr=nf_inq_varid (ncid, vname(1,indxPCO2air_rst)(1:lvar),
     &       vidPCO2air)
        if (ierr .ne. nf_noerr) then
           write(*,1) vname(1,indxPCO2air_rst)(1:lvar), ncname(1:lncn)
           goto 99              !--> ERROR
        endif
# endif /* BIOLOGY_BEC || (BIOLOGY_NPZDOC && OXYGEN && CARBON) */

# if (defined BIOLOGY_BEC || defined BIOLOGY_NPZDOC) && (! defined OUTPUT_AVG_PHYS_BIO || ! defined AVRH)
! PARinc:
        lvar=lenstr(vname(1,indxPARinc_rst))
        ierr=nf_inq_varid (ncid, vname(1,indxPARinc_rst)(1:lvar),
     &       vidPARinc)
        if (ierr .ne. nf_noerr) then
           write(*,1) vname(1,indxPARinc_rst)(1:lvar), ncname(1:lncn)
           goto 99              !--> ERROR
        endif
! PAR:
        lvar=lenstr(vname(1,indxPAR_rst))
        ierr=nf_inq_varid (ncid, vname(1,indxPAR_rst)(1:lvar), vidPAR)
        if (ierr .ne. nf_noerr) then
           write(*,1) vname(1,indxPAR_rst)(1:lvar), ncname(1:lncn)
           goto 99              !--> ERROR
        endif
#endif /* BIOLOGY_BEC || BIOLOGY_NPZDOC */

! the following output is only done for avg files, not his files
# if defined OUTPUT_SPECDOM && defined AVRH && !defined OUTPUT_AVG_PHYS_BIO
! frequency of small phytoplankton dominance at the surface
        lvar=lenstr(vname(1,indxFreqDomSP_sfc))
        ierr=nf_inq_varid (ncid, vname(1,indxFreqDomSP_sfc)(1:lvar), 
     &       vidFreqDomSP_sfc)
        if (ierr .ne. nf_noerr) then
           write(*,1) vname(1,indxFreqDomSP_sfc)(1:lvar), ncname(1:lncn)
           goto 99              !--> ERROR
        endif
! frequency of diatom dominance at the surface
        lvar=lenstr(vname(1,indxFreqDomDIAT_sfc))
        ierr=nf_inq_varid (ncid, vname(1,indxFreqDomDIAT_sfc)(1:lvar), 
     &       vidFreqDomDIAT_sfc)
        if (ierr .ne. nf_noerr) then
           write(*,1) vname(1,indxFreqDomDIAT_sfc)(1:lvar), 
     &          ncname(1:lncn)
           goto 99              !--> ERROR
        endif
! frequency of diazotroph dominance at the surface
        lvar=lenstr(vname(1,indxFreqDomDIAZ_sfc))
        ierr=nf_inq_varid (ncid, vname(1,indxFreqDomDIAZ_sfc)(1:lvar), 
     &       vidFreqDomDIAZ_sfc)
        if (ierr .ne. nf_noerr) then
           write(*,1) vname(1,indxFreqDomDIAZ_sfc)(1:lvar), 
     &          ncname(1:lncn)
           goto 99              !--> ERROR
        endif
! frequency of small phytoplankton dominance, integrated
        lvar=lenstr(vname(1,indxFreqDomSP_int))
        ierr=nf_inq_varid (ncid, vname(1,indxFreqDomSP_int)(1:lvar), 
     &       vidFreqDomSP_int)
        if (ierr .ne. nf_noerr) then
           write(*,1) vname(1,indxFreqDomSP_int)(1:lvar), ncname(1:lncn)
           goto 99              !--> ERROR
        endif
! frequency of diatom dominance, integrated
        lvar=lenstr(vname(1,indxFreqDomDIAT_int))
        ierr=nf_inq_varid (ncid, vname(1,indxFreqDomDIAT_int)(1:lvar), 
     &       vidFreqDomDIAT_int)
        if (ierr .ne. nf_noerr) then
           write(*,1) vname(1,indxFreqDomDIAT_int)(1:lvar), 
     &          ncname(1:lncn)
           goto 99              !--> ERROR
        endif
! frequency of diazotroph dominance, integrated
        lvar=lenstr(vname(1,indxFreqDomDIAZ_int))
        ierr=nf_inq_varid (ncid, vname(1,indxFreqDomDIAZ_int)(1:lvar), 
     &       vidFreqDomDIAZ_int)
        if (ierr .ne. nf_noerr) then
           write(*,1) vname(1,indxFreqDomDIAZ_int)(1:lvar), 
     &          ncname(1:lncn)
           goto 99              !--> ERROR
        endif
# endif /* OUTPUT_SPECDOM */
#endif /* SOLVE3D */
#ifdef MPI_SILENT_MODE
        if (mynode.eq.0) then
#endif
#ifdef AVRH
          write(*,'(7x,4A,I4,2x,A,I4)') 'def_avg :: opened existing ',
#else
          write(*,'(7x,4A,I4,2x,A,I4)') 'def_his :: opened existing ',
#endif
     &   'file ''', ncname(1:lncn), ''' from rec =', rec MYID
#ifdef MPI_SILENT_MODE
        endif
#endif
 
#if defined MPI & !defined PARALLEL_FILES
      else
        ierr=nf_open (ncname(1:lncn), nf_write, ncid)
        if (ierr .eq. nf_noerr) then
# ifndef MPI_SILENT_MODE
#  ifdef AVRH
          write(*,'(7x,5A,I4)')   'def_avg :: reopened existing ',
#  else
          write(*,'(7x,5A,I4)')   'def_his :: reopened existing ',
#  endif
     &            'file ''', ncname(1:lncn), '''.' MYID
# endif
        else
# ifdef AVRH
          write(*,'(/1x,6A,I4/)')  '### ERROR: def_avg :: Cannot ',
# else
          write(*,'(/1x,6A,I4/)')  '### ERROR: def_his :: Cannot ',
# endif
     &                       'open file ''', ncname(1:lncn), '''.',
     &                        nf_strerror(ierr) MYID
          goto 99                                   !--> ERROR
        endif
#endif
      endif             !<--  create_new_file
 
      ierr=nf_set_fill (ncid, nf_nofill, lvar)
      if (ierr .ne. nf_noerr) then

#ifdef AVRH
        write(*,'(1x,4A,I4)') '### ERROR: def_avg :: Cannot switch ',
#else
        write(*,'(1x,4A,I4)') '### ERROR: def_his :: Cannot switch ',
#endif
     &          'to ''nf_nofill'' mode.', nf_strerror(ierr)
     &           MYID
      endif
#ifdef AVRH
   1  format(/1x,'### ERROR: def_avg :: Cannot find variable ''',
#else
   1  format(/1x,'### ERROR: def_his :: Cannot find variable ''',
#endif
     &                           A, ''' in file ''', A, '''.'/)
#if (defined PUT_GRID_INTO_HISTORY && !defined AVRH)\
 || (defined PUT_GRID_INTO_AVERAGES && defined AVRH)

! Write grid variables.
!------ ---- ----------

        if (total_rec.le.1) call wrt_grid (ncid, ncname, lncn)
#endif
  99  return
      end
 
#undef ncid
#undef ncname
#undef rec_per_file
#undef wrt
#undef vidTime
#undef vidZ
#undef vidUb
#undef vidVb
#undef vidU
#undef vidV
#undef vidT
#undef vidR
#undef vidO
#undef vidW
#undef vidAkv
#undef vidAkt
#undef vidAks
#undef vidHbls
#undef vidHbbl
#if defined BULK_FLUX
# undef vidSustr_blk
# undef vidSvstr_blk
# undef vidShflx_net
# undef vidShflx_lat
# undef vidShflx_sen
# undef vidShflx_rad
# undef vidSwflx_emp
# undef vidShflx_wwk
#endif
#if defined BIOLOGY_NPZDOC || defined BIOLOGY_BEC
# undef vidPH
# undef vidPCO2
# undef vidPCO2air
# undef vidPARinc
# undef vidPAR
#  ifdef SEDIMENT_BIOLOGY
#   undef vidTsed
#  endif
#  if defined OUTPUT_SPECDOM && defined AVRH
#   undef vidFreqDomSP_sfc
#   undef vidFreqDomDIAT_sfc
#   undef vidFreqDomDIAZ_sfc
#   undef vidFreqDomSP_int
#   undef vidFreqDomDIAT_int
#   undef vidFreqDomDIAZ_int
#  endif /* OUTPUT_SPECDOM */
#endif /* BIOLOGY_NPZDOC || BIOLOGY_BEC */
#ifdef WKB_WWAVE
# undef vidWKB
#endif
#ifdef MRL_WCI
# undef vidSUP
# undef vidUST2D
# undef vidVST2D
# undef vidUST
# undef vidVST
# undef vidWST
# undef vidAkb
# undef vidAkw
# undef vidKVF
# undef vidCALP
# undef vidKAPS
#endif
 
#ifndef AVRH
# define AVRH
# ifdef AVERAGES
#  include "def_his.F"
# endif
#endif
