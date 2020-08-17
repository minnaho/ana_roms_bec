!
!
#define AVRH
#ifndef AVRH
# include "cppdefs.h"
!hf: do NOT use his output at this time
!WARN      subroutine def_bio_his (ncid, total_rec, ierr)
#else
# include "cppdefs.h"
#ifdef OUTPUT_AVG_PHYS_BIO
      subroutine def_bio_avg (ncid, total_rec, ierr)
!
! Create/open averages/history netCDF file. In the case when a new
! netCDF file is created, define all variables, their dimensions and
! attributes. In the case when a previously existing netCDF file is
! to be opened for addition of new data, verify that all dimensions
! of the file are consistent with the present model configuration
! and all necessary variables exist. Save netCDF IDs for all needed
! variables. Also determine size of the unlimited dimension.
!
! The difference between def_his and def_avg is as follows: they
! have different netCDF file name (hisname/avgname); netCDF file ID
! (passed as argument); time record index (hisindx/avgindx); array
! of switches which variables to write (wrthis/wrtavg); and different
! sets of netCDF variable IDs (hisTime...hisHbl/avgTime...avgHbl);
! and the first attribute of each variable, long_name, has prefix
! 'averaged'. Because most of the code is identical for both
! routines, the second one is generated from the first entirely
! by CPP.
!
      implicit none
      logical create_new_file
      integer ncid, total_rec, ierr, rec, lstr,lvar,lenstr, timedim
     &      , r2dgrd(3),  u2dgrd(3), v2dgrd(3), auxil(2), checkdims
     &      , b3dgrd(4)
#ifdef SOLVE3D
     &      , r3dgrd(4),  u3dgrd(4), v3dgrd(4), w3dgrd(4), itrc
#endif
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
# define ncname hisname
# define rec_per_file nrpfhis
# define wrt wrthis
# define vidTime hisTime
# define vidTstep hisTstep
# define vidZ hisZ
# define vidT hisT
# define vidPH hisPH
# define vidPCO2 hisPCO2
# define vidPCO2air hisPCO2air
# define vidPARinc hisPARinc
# define vidPAR hisPAR
# ifdef SEDIMENT_BIOLOGY
# define vidTsed hisTsed
# endif
# ifdef BIOLOGY
#  define vidHel hisHel
#  ifdef BIO_NPZD
#   define vidChC hisChC
#   ifdef OXYGEN
#     define vidU10 hisU10
#     define vidKvO2 hisKvO2
#     define vidO2sat hisO2sat
#   endif
#  elif defined BIO_N2P2Z2D2 && defined VAR_CHL_C
#   define vidChC1 hisChC1
#   define vidChC2 hisChC2
#  endif
# endif
# define vidSed hisSed
#else
# define ncname bio_avgname
# define rec_per_file nrpf_bio_avg
# define wrt wrtavg
# define vidTime avg_bio_time
# define vidTstep avg_bio_Tstep
# define vidZ avg_bio_Z
# define vidT avgT
#define vidPH avgPH
#define vidPCO2 avgPCO2
#define vidPCO2air avgPCO2air
#define vidPARinc avgPARinc
#define vidPAR avgPAR
# ifdef SEDIMENT_BIOLOGY
# define vidTsed avgTsed
# endif
# ifdef OUTPUT_SPECDOM
#  define vidFreqDomSP_sfc avgFreqDomSP_sfc
#  define vidFreqDomDIAT_sfc avgFreqDomDIAT_sfc
#  define vidFreqDomDIAZ_sfc avgFreqDomDIAZ_sfc
#  define vidFreqDomSP_int avgFreqDomSP_int
#  define vidFreqDomDIAT_int avgFreqDomDIAT_int
#  define vidFreqDomDIAZ_int avgFreqDomDIAZ_int
# endif /* OUTPUT_SPECDOM */
# ifdef BIOLOGY
#  define vidHel avgHel
#  ifdef BIO_NPZD
#   define vidChC avgChC
#   ifdef OXYGEN
#    define vidU10 avgU10
#    define vidKvO2 avgKvO2
#    define vidO2sat avgO2sat
#   endif
#  elif defined BIO_N2P2Z2D2 && defined VAR_CHL_C
#   define vidChC1 avgChC1
#   define vidChC2 avgChC2
#  endif
# endif
      character*70 text
#endif
!
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
!
      ierr=0
      lstr=lenstr(ncname)
      if (rec_per_file.gt.0) then
        lvar=total_rec-(1+mod(total_rec-1, rec_per_file))
        call insert_time_index (ncname, lstr, lvar, ierr)
        if (ierr .ne. 0) goto 99
      endif
!
! Decide whether to create a new file, or open existing one.
! Overall the whole code below is organized into 3-way switch,
!
! 10  if (create_new_file) then
!        .... create new file, save netCDF ids for all variables;
!     elseif (ncid.eq.-1) then
!        .... try to open existing file and check its dimensions
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
! which is designed to implement flexible opening policy: 
! if ldefhis=.true., it forces creation of a new file [if the
! file already exists, it will be overwritten]; on the other hand,
! ldefhis=.false., it is assumed that the file already exists and
! an attempt to open it is made; if the attempt is successful, the
! file is prepared for appending hew data; if it fails, a new file
! is created.
!
      create_new_file=.true.
      if (ncid.ne.-1) create_new_file=.false.
#if defined MPI && !defined PARALLEL_FILES
      if (mynode.gt.0) create_new_file=.false.
#endif
!
! Create new history/averages file:    Put global attributes
! ====== === ======= ======== =====    and define all variables.
!
  10  if (create_new_file) then
        ierr=nf_create(ncname(1:lstr), nf_clobber, ncid)
        if (ierr .ne. nf_noerr) then
          write(*,'(/3(1x,A)/)') 'ERROR in def_bio_avg:',
     &           'Cannot create netCDF file:', ncname(1:lstr)
          goto 99                                         !--> ERROR
        endif
        if (rec_per_file.eq.0) total_rec=0  
!
! Put global attributes.
! --- ------ -----------
!
        call put_global_atts (ncid, ierr)
!
! Define dimensions of staggered fields.
! ------ ---------- -- --------- -------
!
        ierr=nf_def_dim (ncid, 'xi_rho',   xi_rho,   r2dgrd(1))
        ierr=nf_def_dim (ncid, 'xi_u',     xi_u,     u2dgrd(1))
        ierr=nf_def_dim (ncid, 'eta_rho',  eta_rho,  r2dgrd(2))
        ierr=nf_def_dim (ncid, 'eta_v',    eta_v,    v2dgrd(2))
#ifdef SEDIMENT
        ierr=nf_def_dim (ncid, 's_b',      NLAY,     b3dgrd(3))
#endif
#ifdef SOLVE3D
        ierr=nf_def_dim (ncid, 's_rho',    N,        r3dgrd(3))
        ierr=nf_def_dim (ncid, 's_w',      N+1,      w3dgrd(3))
#endif
        ierr=nf_def_dim (ncid, 'time', nf_unlimited, timedim)
        ierr=nf_def_dim (ncid, 'auxil',    4,        auxil(1))
        
        auxil(2)=timedim

        r2dgrd(3)=timedim           ! Free surface

        u2dgrd(2)=r2dgrd(2)         ! 2D UBAR-type
        u2dgrd(3)=timedim

        v2dgrd(1)=r2dgrd(1)         ! 2D VBAR-type
        v2dgrd(3)=timedim
        
        b3dgrd(1)=r2dgrd(1)         !
        b3dgrd(2)=r2dgrd(2)         ! 3D BED-type
        b3dgrd(4)=timedim           !

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
!
! Define grid variables.
! ------ ---- ----------
!
        if (total_rec.le.1) call def_grid (ncid, r2dgrd)
#endif

!
! Define evolving model variables.
! ------ -------- ----- ----------
!
! Time step number and time record indices:
!
        ierr=nf_def_var (ncid, 'time_step', nf_int, 2, auxil,
     &                                                 vidTstep)
        ierr=nf_put_att_text (ncid, vidTstep, 'long_name', 48,
     &       'time step and record numbers from initialization')
!
! Time.
!
        lvar=lenstr(vname(1,indxTime))
! hf: always use double precision for the time variable
        ierr=nf_def_var (ncid, vname(1,indxTime)(1:lvar),
     &                            NF_DOUBLE, 1, timedim, vidTime)
#ifndef AVRH
        lvar=lenstr(vname(2,indxTime))
        ierr=nf_put_att_text (ncid, vidTime, 'long_name', lvar,
     &                                vname(2,indxTime)(1:lvar))
#else
        text='averaged '/ /vname(2,indxTime)
        lvar=lenstr(text)
        ierr=nf_put_att_text (ncid, vidTime, 'long_name', lvar,
     &                                             text(1:lvar))
#endif
        lvar=lenstr(vname(3,indxTime))
        ierr=nf_put_att_text (ncid, vidTime, 'units',  lvar,
     &                                vname(3,indxTime)(1:lvar))
!
! Free-surface.
!
        if (wrt(indxZ)) then
          lvar=lenstr(vname(1,indxZ))
          ierr=nf_def_var (ncid, vname(1,indxZ)(1:lvar),
     &                              NF_FOUT, 3, r2dgrd, vidZ)
#ifndef AVRH
          lvar=lenstr(vname(2,indxZ))
          ierr=nf_put_att_text (ncid, vidZ, 'long_name', lvar,
     &                                  vname(2,indxZ)(1:lvar))
#else
          text='averaged '/ /vname(2,indxZ)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidZ, 'long_name', lvar,
     &                                            text(1:lvar))
#endif
          lvar=lenstr(vname(3,indxZ))
          ierr=nf_put_att_text (ncid, vidZ, 'units',     lvar,
     &                                  vname(3,indxZ)(1:lvar))
#ifdef MASK_LAND_DATA
          ierr=nf_put_att_FOUT (ncid, vidZ, '_FillValue', NF_FOUT,
     &                                                 1, spv_set)
#endif
        endif
#ifdef SOLVE3D
!
! Tracer variables.
!
        do itrc=2+ntrc_salt,NT
          if (wrt(indxT+itrc-1)) then
            lvar=lenstr(vname(1,indxT+itrc-1))
            ierr=nf_def_var (ncid, vname(1,indxT+itrc-1)(1:lvar),
     &                             NF_FOUT, 4, r3dgrd, vidT(itrc))
# ifndef AVRH
            lvar=lenstr(vname(2,indxT+itrc-1))
            ierr=nf_put_att_text (ncid, vidT(itrc), 'long_name',
     &                         lvar, vname(2,indxT+itrc-1)(1:lvar))
# else
            text='averaged '/ /vname(2,indxT+itrc-1)
            lvar=lenstr(text)
            ierr=nf_put_att_text (ncid, vidT(itrc), 'long_name',
     &                                          lvar, text(1:lvar))
# endif
            lvar=lenstr(vname(3,indxT+itrc-1))
            ierr=nf_put_att_text (ncid, vidT(itrc), 'units', lvar,
     &                               vname(3,indxT+itrc-1)(1:lvar))
#ifdef MASK_LAND_DATA
            ierr=nf_put_att_FOUT (ncid, vidT(itrc), '_FillValue', 
     &           NF_FOUT, 1, spv_set)
#endif
          endif
        enddo
# ifdef SEDIMENT_BIOLOGY
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
#ifdef MASK_LAND_DATA
            ierr=nf_put_att_FOUT (ncid, vidTsed(itrc), '_FillValue', 
     &           NF_FOUT, 1, spv_set)
#endif
          endif
        enddo
# endif /* SEDIMENT_BIOLOGY */

# if defined BIOLOGY_BEC || (defined BIOLOGY_NPZDOC && \
     defined OXYGEN && defined CARBON)

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

# if defined BIOLOGY_BEC || defined BIOLOGY_NPZDOC
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
# if defined OUTPUT_SPECDOM
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


# ifdef BIOLOGY
!
! Depth of the euphotic layer.
!
        if (wrt(indxHel)) then
          lvar=lenstr(vname(1,indxHel))
          ierr=nf_def_var (ncid, vname(1,indxHel)(1:lvar),
     &                             NF_FOUT, 3, r2dgrd, vidHel)
#  ifndef AVRH
          lvar=lenstr(vname(2,indxHel))
          ierr=nf_put_att_text (ncid, vidHel, 'long_name', lvar,
     &                                  vname(2,indxHel)(1:lvar))
#  else
          text='averaged '/ /vname(2,indxHel)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidHel, 'long_name', lvar,
     &                                              text(1:lvar))
#  endif
          lvar=lenstr(vname(3,indxHel))
          ierr=nf_put_att_text (ncid, vidHel, 'units',     lvar,
     &                                  vname(3,indxHel)(1:lvar))
        endif
!
! Chlorophyll/Carbon ratio
!
#  ifdef BIO_NPZD
        if (wrt(indxChC)) then
          lvar=lenstr(vname(1,indxChC))
          ierr=nf_def_var (ncid, vname(1,indxChC)(1:lvar),
     &                             NF_FOUT, 4, r3dgrd, vidChC)
#   ifndef AVRH
          lvar=lenstr(vname(2,indxChC))
          ierr=nf_put_att_text (ncid, vidChC, 'long_name', lvar,
     &                                  vname(2,indxChC)(1:lvar))
#   else
          text='averaged '/ /vname(2,indxChC)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidChC, 'long_name', lvar,
     &                                              text(1:lvar))
#   endif
          lvar=lenstr(vname(3,indxChC))
          ierr=nf_put_att_text (ncid, vidChC, 'units',     lvar,
     &                                  vname(3,indxChC)(1:lvar))
        endif
#   ifdef OXYGEN
!
! Wind speed at 10 m
!
        if (wrt(indxU10)) then
          lvar=lenstr(vname(1,indxU10))
          ierr=nf_def_var (ncid, vname(1,indxU10)(1:lvar),
     &                             NF_FOUT, 3, r2dgrd, vidU10)
#  ifndef AVRH
          lvar=lenstr(vname(2,indxU10))
          ierr=nf_put_att_text (ncid, vidU10, 'long_name', lvar,
     &                                  vname(2,indxU10)(1:lvar))
#  else
          text='averaged '/ /vname(2,indxU10)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidU10, 'long_name', lvar,
     &                                              text(1:lvar))
#  endif
          lvar=lenstr(vname(3,indxU10))
          ierr=nf_put_att_text (ncid, vidU10, 'units',     lvar,
     &                                  vname(3,indxU10)(1:lvar))
        endif
!
! Gas exchange coefficient of O2
!
        if (wrt(indxKvO2)) then
          lvar=lenstr(vname(1,indxKvO2))
          ierr=nf_def_var (ncid, vname(1,indxKvO2)(1:lvar),
     &                             NF_FOUT, 3, r2dgrd, vidKvO2)
#  ifndef AVRH
          lvar=lenstr(vname(2,indxKvO2))
          ierr=nf_put_att_text (ncid, vidKvO2, 'long_name', lvar,
     &                                  vname(2,indxKvO2)(1:lvar))
#  else
          text='averaged '/ /vname(2,indxKvO2)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidKvO2, 'long_name', lvar,
     &                                              text(1:lvar))
#  endif
          lvar=lenstr(vname(3,indxKvO2))
          ierr=nf_put_att_text (ncid, vidKvO2, 'units',     lvar,
     &                                  vname(3,indxKvO2)(1:lvar))
        endif
!
! Saturation concentration of O2
!
        if (wrt(indxO2sat)) then
          lvar=lenstr(vname(1,indxO2sat))
          ierr=nf_def_var (ncid, vname(1,indxO2sat)(1:lvar),
     &                             NF_FOUT, 3, r2dgrd, vidO2sat)
#  ifndef AVRH
          lvar=lenstr(vname(2,indxO2sat))
          ierr=nf_put_att_text (ncid, vidO2sat, 'long_name', lvar,
     &                                  vname(2,indxO2sat)(1:lvar))
#  else
          text='averaged '/ /vname(2,indxO2sat)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidO2sat, 'long_name', lvar,
     &                                              text(1:lvar))
#  endif
          lvar=lenstr(vname(3,indxO2sat))
          ierr=nf_put_att_text (ncid, vidO2sat, 'units',     lvar,
     &                                  vname(3,indxO2sat)(1:lvar))
        endif
#   endif /* OXYGEN */
#  elif defined BIO_N2P2Z2D2 && defined VAR_CHL_C
        if (wrt(indxChC1)) then
          lvar=lenstr(vname(1,indxChC1))
          ierr=nf_def_var (ncid, vname(1,indxChC1)(1:lvar),
     &                             NF_FOUT, 4, r3dgrd, vidChC1)
#   ifndef AVRH
          lvar=lenstr(vname(2,indxChC1))
          ierr=nf_put_att_text (ncid, vidChC1, 'long_name', lvar,
     &                                  vname(2,indxChC1)(1:lvar))
#   else
          text='averaged '/ /vname(2,indxChC1)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidChC1, 'long_name', lvar,
     &                                              text(1:lvar))
#   endif
          lvar=lenstr(vname(3,indxChC1))
          ierr=nf_put_att_text (ncid, vidChC1, 'units',     lvar,
     &                                  vname(3,indxChC1)(1:lvar))
        endif

        if (wrt(indxChC2)) then
          lvar=lenstr(vname(1,indxChC2))
          ierr=nf_def_var (ncid, vname(1,indxChC2)(1:lvar),
     &                             NF_FOUT, 4, r3dgrd, vidChC2)
#   ifndef AVRH
          lvar=lenstr(vname(2,indxChC2))
          ierr=nf_put_att_text (ncid, vidChC2, 'long_name', lvar,
     &                                  vname(2,indxChC2)(1:lvar))
#   else
          text='averaged '/ /vname(2,indxChC2)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidChC2, 'long_name', lvar,
     &                                              text(1:lvar))
#   endif
          lvar=lenstr(vname(3,indxChC2))
          ierr=nf_put_att_text (ncid, vidChC2, 'units',     lvar,
     &                                  vname(3,indxChC2)(1:lvar))
        endif
#  endif
# endif
!
# ifdef SEDIMENT
!
! bed_thick  Sediment bed layer thickness (m)
! bed_poros  Bed porosity (-)
! bed_frac_sand/silt
!           volume fraction of size class ised in each bed layer (-)
!
        do itrc=1,NST+2
         indxWrk=indxSed+itrc-1
         if (wrt(indxWrk)) then
          lvar=lenstr(vname(1,indxWrk))
          ierr=nf_def_var (ncid, vname(1,indxWrk)(1:lvar),
     &                             NF_FOUT, 4, b3dgrd, vidSed(itrc))
#  ifndef AVRH
          lvar=lenstr(vname(2,indxWrk))
          ierr=nf_put_att_text (ncid, vidSed(itrc), 'long_name', lvar,
     &                                  vname(2,indxWrk)(1:lvar))
#  else
          text='averaged '/ /vname(2,indxWrk)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidSed(itrc), 'long_name', lvar,
     &                                              text(1:lvar))
#  endif
          lvar=lenstr(vname(3,indxWrk))
          ierr=nf_put_att_text (ncid, vidSed(itrc), 'units',     lvar,
     &                                  vname(3,indxWrk)(1:lvar))
         endif
        enddo
!
# endif /* SEDIMENT */
#endif /* SOLVE3D */
!
! Leave definition mode.
! ----- ---------- -----
!
        ierr=nf_enddef(ncid)
        write(*,'(6x,4A,1x,A,i4)') 'DEF_BIO_AVG - Created ',
     &                'new netCDF file ''', ncname(1:lstr), '''.'
     &                 MYID
!
! Open an existing file and prepare for appending data.
! ==== == ======== ==== === ======= === ========= =====
! Inquire about the dimensions and variables. Check for
! consistency with model dimensions. In the case when file 
! is rejected (whether it cannot be opened, or something
! is wrong with its dimensions) create a new file.
!
! After that verify that all necessary variables are already
! defined there and find their netCDF IDs.
!
      elseif (ncid.eq.-1) then
        ierr=nf_open (ncname(1:lstr), nf_write, ncid)
        if (ierr. eq. nf_noerr) then
          ierr=checkdims (ncid, ncname(1:lstr), rec) 
          if (ierr .eq. nf_noerr) then
            if (rec_per_file.eq.0) then
              ierr=rec+1 - total_rec
            else
              ierr=rec+1 - (1+mod(total_rec-1, rec_per_file))
            endif
            if (ierr.gt.0) then
              mpi_master_only write( *,
     &                 '(/1x,A,I5,1x,A/8x,3A,I5,/8x,A,I5,1x,A/)'
     &            ) 'DEF_BIO_AVG WARNING: Actual number of records',
     &               rec,  'in netCDF file',  '''',  ncname(1:lstr),
     &             ''' exceeds the record number from restart data',
     &             rec+1-ierr,'/', total_rec,', restart is assumed.'
              rec=rec-ierr
            elseif (rec_per_file.eq.0) then
              total_rec=rec+1           ! <-- set to the next record
#if defined MPI && !defined PARALLEL_FILES
              if (mynode.gt.0) total_rec=total_rec-1
#endif
            endif
            ierr=nf_noerr
          endif
        endif

        if (ierr. ne. nf_noerr) then
#if defined MPI && !defined PARALLEL_FILES
          if (mynode.eq.0) then
            create_new_file=.true.
            goto 10
          else
            write(*,'(/1x,4A,2x,A,I4/)') 'DEF_BIO_AVG ERROR: ',
     &                  'Cannot open file ''', ncname(1:lstr), '''.'
     &                   MYID
            goto 99                                       !--> ERROR
          endif
#else
          create_new_file=.true.
          goto 10
#endif
        endif
!
! Find netCDF IDs of evolving model variables:
! ---- ------ --- -- -------- ----- ----------
!
! Time step indices:
!
        ierr=nf_inq_varid (ncid, 'time_step', vidTstep)
        if (ierr .ne. nf_noerr) then
          write(*,1) 'time_step', ncname(1:lstr)
          goto 99                                         !--> ERROR
        endif
!
! Time.
!
        lvar=lenstr(vname(1,indxTime))
        ierr=nf_inq_varid (ncid,vname(1,indxTime)(1:lvar),vidTime)
        if (ierr .ne. nf_noerr) then
          write(*,1) vname(1,indxTime)(1:lvar), ncname(1:lstr)
          goto 99                                         !--> ERROR 
        endif
!
! Free-surface.
!
        if (wrt(indxZ)) then
          lvar=lenstr(vname(1,indxZ))
          ierr=nf_inq_varid (ncid, vname(1,indxZ)(1:lvar), vidZ)
          if (ierr .ne. nf_noerr) then
            write(*,1) vname(1,indxZ)(1:lvar), ncname(1:lstr)
            goto 99                                       !--> ERROR
          endif
        endif
#ifdef SOLVE3D
!
! Tracer variables.
!
        do itrc=2+ntrc_salt,NT
          if (wrt(indxT+itrc-1)) then
            lvar=lenstr(vname(1,indxT+itrc-1))
            ierr=nf_inq_varid (ncid, vname(1,indxT+itrc-1)(1:lvar),
     &                                                 vidT(itrc))
            if (ierr .ne. nf_noerr) then
              write(*,1) vname(1,indxT+itrc-1)(1:lvar),
     &                                       ncname(1:lstr) 
              goto 99                                     !--> ERROR
            endif
          endif
        enddo
# ifdef SEDIMENT_BIOLOGY
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
     &                                       ncname(1:lstr) 
              goto 99                                     !--> ERROR
            endif
          endif
        enddo
# endif /* SEDIMENT_BIOLOGY */

# if defined BIOLOGY_BEC || (defined BIOLOGY_NPZDOC && \
     defined OXYGEN && defined CARBON)

! pH
        lvar=lenstr(vname(1,indxPH_rst))
        ierr=nf_inq_varid (ncid, vname(1,indxPH_rst)(1:lvar), vidPH)
        if (ierr .ne. nf_noerr) then
           write(*,1) vname(1,indxPH_rst)(1:lvar), ncname(1:lstr)
           goto 99              !--> ERROR
        endif
! pCO2:
        lvar=lenstr(vname(1,indxPCO2_rst))
        ierr=nf_inq_varid (ncid, vname(1,indxPCO2_rst)(1:lvar), 
     &       vidPCO2)
        if (ierr .ne. nf_noerr) then
           write(*,1) vname(1,indxPCO2_rst)(1:lvar), ncname(1:lstr)
           goto 99              !--> ERROR
        endif
! pCO2air:
        lvar=lenstr(vname(1,indxPCO2air_rst))
        ierr=nf_inq_varid (ncid, vname(1,indxPCO2air_rst)(1:lvar),
     &       vidPCO2air)
        if (ierr .ne. nf_noerr) then
           write(*,1) vname(1,indxPCO2air_rst)(1:lvar), ncname(1:lstr)
           goto 99              !--> ERROR
        endif
# endif /* BIOLOGY_BEC || (BIOLOGY_NPZDOC && OXYGEN && CARBON) */
# if defined BIOLOGY_BEC || defined BIOLOGY_NPZDOC
! PARinc:
        lvar=lenstr(vname(1,indxPARinc_rst))
        ierr=nf_inq_varid (ncid, vname(1,indxPARinc_rst)(1:lvar),
     &       vidPARinc)
        if (ierr .ne. nf_noerr) then
           write(*,1) vname(1,indxPARinc_rst)(1:lvar), ncname(1:lstr)
           goto 99              !--> ERROR
        endif
! PAR:
        lvar=lenstr(vname(1,indxPAR_rst))
        ierr=nf_inq_varid (ncid, vname(1,indxPAR_rst)(1:lvar), vidPAR)
        if (ierr .ne. nf_noerr) then
           write(*,1) vname(1,indxPAR_rst)(1:lvar), ncname(1:lstr)
           goto 99              !--> ERROR
        endif
#endif /* BIOLOGY_BEC || BIOLOGY_NPZDOC */

! the following output is only done for avg files, not his files
# if defined OUTPUT_SPECDOM
! frequency of small phytoplankton dominance at the surface
        lvar=lenstr(vname(1,indxFreqDomSP_sfc))
        ierr=nf_inq_varid (ncid, vname(1,indxFreqDomSP_sfc)(1:lvar), 
     &       vidFreqDomSP_sfc)
        if (ierr .ne. nf_noerr) then
           write(*,1) vname(1,indxFreqDomSP_sfc)(1:lvar), ncname(1:lstr)
           goto 99              !--> ERROR
        endif
! frequency of diatom dominance at the surface
        lvar=lenstr(vname(1,indxFreqDomDIAT_sfc))
        ierr=nf_inq_varid (ncid, vname(1,indxFreqDomDIAT_sfc)(1:lvar), 
     &       vidFreqDomDIAT_sfc)
        if (ierr .ne. nf_noerr) then
           write(*,1) vname(1,indxFreqDomDIAT_sfc)(1:lvar), 
     &          ncname(1:lstr)
           goto 99              !--> ERROR
        endif
! frequency of diazotroph dominance at the surface
        lvar=lenstr(vname(1,indxFreqDomDIAZ_sfc))
        ierr=nf_inq_varid (ncid, vname(1,indxFreqDomDIAZ_sfc)(1:lvar), 
     &       vidFreqDomDIAZ_sfc)
        if (ierr .ne. nf_noerr) then
           write(*,1) vname(1,indxFreqDomDIAZ_sfc)(1:lvar), 
     &          ncname(1:lstr)
           goto 99              !--> ERROR
        endif
! frequency of small phytoplankton dominance, integrated
        lvar=lenstr(vname(1,indxFreqDomSP_int))
        ierr=nf_inq_varid (ncid, vname(1,indxFreqDomSP_int)(1:lvar), 
     &       vidFreqDomSP_int)
        if (ierr .ne. nf_noerr) then
           write(*,1) vname(1,indxFreqDomSP_int)(1:lvar), ncname(1:lstr)
           goto 99              !--> ERROR
        endif
! frequency of diatom dominance, integrated
        lvar=lenstr(vname(1,indxFreqDomDIAT_int))
        ierr=nf_inq_varid (ncid, vname(1,indxFreqDomDIAT_int)(1:lvar), 
     &       vidFreqDomDIAT_int)
        if (ierr .ne. nf_noerr) then
           write(*,1) vname(1,indxFreqDomDIAT_int)(1:lvar), 
     &          ncname(1:lstr)
           goto 99              !--> ERROR
        endif
! frequency of diazotroph dominance, integrated
        lvar=lenstr(vname(1,indxFreqDomDIAZ_int))
        ierr=nf_inq_varid (ncid, vname(1,indxFreqDomDIAZ_int)(1:lvar), 
     &       vidFreqDomDIAZ_int)
        if (ierr .ne. nf_noerr) then
           write(*,1) vname(1,indxFreqDomDIAZ_int)(1:lvar), 
     &          ncname(1:lstr)
           goto 99              !--> ERROR
        endif
# endif /* OUTPUT_SPECDOM */

# ifdef BIOLOGY
!
! Depth of the euphotic layer.
!
        if (wrt(indxHel)) then
          lvar=lenstr(vname(1,indxHel))
          ierr=nf_inq_varid (ncid,vname(1,indxHel)(1:lvar), vidHel)
          if (ierr .ne. nf_noerr) then
          write(*,1) vname(1,indxHel)(1:lvar), ncname(1:lstr)
            goto 99                                       !--> ERROR
          endif
        endif
!
! Chlorophyll/Carbon ratio.
!
#  ifdef BIO_NPZD
        if (wrt(indxChC)) then
          lvar=lenstr(vname(1,indxChC))
          ierr=nf_inq_varid (ncid,vname(1,indxChC)(1:lvar), vidChC)
          if (ierr .ne. nf_noerr) then
          write(*,1) vname(1,indxChC)(1:lvar), ncname(1:lstr)
            goto 99                                       !--> ERROR
          endif
        endif
#   ifdef OXYGEN
!
! Wind speed at 10 m
        if (wrt(indxU10)) then
          lvar=lenstr(vname(1,indxU10))
          ierr=nf_inq_varid (ncid,vname(1,indxU10)(1:lvar), vidU10)
          if (ierr .ne. nf_noerr) then
          write(*,1) vname(1,indxU10)(1:lvar), ncname(1:lstr)
            goto 99                                       !--> ERROR
          endif
        endif
!
! Gas exchange coefficient of O2
        if (wrt(indxKvO2)) then
          lvar=lenstr(vname(1,indxKvO2))
          ierr=nf_inq_varid (ncid,vname(1,indxKvO2)(1:lvar), vidChC)
          if (ierr .ne. nf_noerr) then
          write(*,1) vname(1,indxKvO2)(1:lvar), ncname(1:lstr)
            goto 99                                       !--> ERROR
          endif
        endif
!
! Saturation concentration of O2
        if (wrt(indxO2sat)) then
          lvar=lenstr(vname(1,indxO2sat))
          ierr=nf_inq_varid (ncid,vname(1,indxO2sat)(1:lvar), vidChC)
          if (ierr .ne. nf_noerr) then
          write(*,1) vname(1,indxO2sat)(1:lvar), ncname(1:lstr)
            goto 99                                       !--> ERROR
          endif
        endif
#   endif /* OXYGEN */
#  elif defined BIO_N2P2Z2D2 && defined VAR_CHL_C
        if (wrt(indxChC1)) then
          lvar=lenstr(vname(1,indxChC1))
          ierr=nf_inq_varid (ncid,vname(1,indxChC1)(1:lvar), vidChC1)
          if (ierr .ne. nf_noerr) then
          write(*,1) vname(1,indxChC1)(1:lvar), ncname(1:lstr)
            goto 99                                       !--> ERROR
          endif
        endif

        if (wrt(indxChC2)) then
          lvar=lenstr(vname(1,indxChC2))
          ierr=nf_inq_varid (ncid,vname(1,indxChC2)(1:lvar), vidChC2)
          if (ierr .ne. nf_noerr) then
          write(*,1) vname(1,indxChC2)(1:lvar), ncname(1:lstr)
            goto 99                                       !--> ERROR
          endif
        endif
#  endif
# endif
# ifdef SEDIMENT
!
! Sediment bed layer thickness, porosity, and volume fraction of 
! size class ised in sediment layer (2+NST r2dgrd variables)
!
        do itrc=1,NST+2
         indxWrk=indxSed+itrc-1
         if (wrt(indxWrk)) then
          lvar=lenstr(vname(1,indxWrk))
          ierr=nf_inq_varid (ncid,vname(1,indxWrk)(1:lvar), 
     &                                                vidSed(itrc))
          if (ierr .ne. nf_noerr) then
          write(*,1) vname(1,indxWrk)(1:lvar), ncname(1:lstr)
            goto 99                                       !--> ERROR
          endif
         endif
        enddo
# endif /*SEDIMENT */

#endif /* SOLVE3D */
        write(*,'(6x,2A,i4,1x,A,i4)') 'DEF_BIO_AVG -- Opened ',
     &                     'existing file  from record =', rec
     &                      MYID

#if defined MPI && !defined PARALLEL_FILES
      else
        ierr=nf_open (ncname(1:lstr), nf_write, ncid)
        if (ierr .ne. nf_noerr) then
          write(*,'(/1x,4A,2x,A,I4/)') 'DEF_BIO_AVG ERROR: ',
     &                'Cannot open file ''', ncname(1:lstr), '''.'
     &                 MYID
          goto 99                                         !--> ERROR
        endif
#endif
      endif             !<--  create_new_file  
      
      ierr=nf_set_fill (ncid, nf_nofill, lvar)
      if (ierr .ne. nf_noerr) then
        write(*,'(6x,2A,i4,1x,A,i4)') 'DEF_BIO_HIS/AVG ERROR: Cannot ', 
     &    'switch to ''nf_nofill'' more; netCDF error code =', ierr 
      endif
   1  format(/1x,'DEF_BIO_HIS/AVG ERROR: Cannot find variable ''',
     &                   A, ''' in netCDF file ''', A, '''.'/)
#if (defined PUT_GRID_INTO_HISTORY && !defined AVRH)\
 || (defined PUT_GRID_INTO_AVERAGES && defined AVRH)
!
! Write grid variables.
! ----- ---- ----------
!
        if (total_rec.le.1) call wrt_grid (ncid, ncname, lstr)
#endif
#if defined AVRH
!
! Set initial average time.
! --- ------- ------- -----
!
      if (nts_bio_avg.eq.1) then
        time_bio_avg=time-0.5*float(n_bio_avg)*dt
      else
        time_bio_avg=time-0.5*float(n_bio_avg)*dt+float(nts_bio_avg)*dt
      endif
#endif

  99  return
      end

#undef ncname
#undef rec_per_file
#undef wrt 
#undef vidTime
#undef vidTstep
#undef vidZ
#undef vidT
# ifdef SEDIMENT_BIOLOGY
# undef vidTsed
# endif
#ifdef BIOLOGY
# undef vidHel
# ifdef BIO_NPZD
#  undef vidChC
#  ifdef OXYGEN
#   undef vidU10
#   undef vidKvO2
#   undef vidO2sat
#  endif
# elif defined BIO_N2P2Z2D2 && defined VAR_CHL_C
#  undef vidChC1
#  undef vidChC2
# endif
#endif
#undef vidSed

#ifdef AVERAGES
# ifndef AVRH
!!#  define AVRH
!!#  include "def_his.F"
# endif
#endif
#else /* OUTPUT_AVG_PHYS_BIO */
      subroutine def_bio_avg_empty
      return
      end
#endif /* OUTPUT_AVG_PHYS_BIO */
#endif /* AVRH */
