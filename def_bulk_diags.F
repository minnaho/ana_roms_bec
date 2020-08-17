#ifndef AVRH
# include "cppdefs.h"
#endif
#ifdef BULK_DIAGS
#ifndef AVRH
      subroutine def_bulk_diags_his(ncid, total_rec, ierr)
#else
      subroutine def_bulk_diags_avg(ncid, total_rec, ierr)
#endif
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
!      USE scalars
!      USE ncscrum

      implicit none
      logical create_new_file
      integer ncid, total_rec, ierr, rec, lstr,lvar,lenstr, timedim
     &      , r2dgrd(3),  u2dgrd(3), v2dgrd(3),  auxil(2),  checkdims
#include "netcdf.inc"
#include "param.h"
#include "scalars.h"
#include "ncvars.h"
#ifdef BULK_DIAGS
# include "forces.h"
#endif

#ifdef MASK_LAND_DATA
# ifdef HIS_DOUBLE
      real*8 spv_set
      parameter (spv_set=1.D+33)
# else
      real*4 spv_set
      parameter (spv_set=1.E+33)
# endif
#endif

#ifndef AVRH
# define ncname bulk_diags_his_name
# define rec_per_file nrpf_bulk_diags_his
# define vidTime bulk_diags_hisTime
# define vidTstep bulk_diags_hisTstep
# define vidZ bulk_diags_hisZ
# define vidSustr_blk hisSustr_blk
# define vidSvstr_blk hisSvstr_blk
# define vidShflx_net hisShflx_net
# define vidShflx_lat hisShflx_lat
# define vidShflx_sen hisShflx_sen
# define vidShflx_rad hisShflx_rad
# define vidSwflx_emp hisSwflx_emp
# define vidShflx_lat hisShflx_lat
# define vidShflx_wwk hisShflx_wwk
# define vidSurf_u hisSurf_u
# define vidSurf_v hisSurf_v

#else /* AVRH */
# define ncname bulk_diags_avg_name
# define rec_per_file nrpf_bulk_diags_avg
# define vidTime bulk_diags_avgTime
# define vidTstep bulk_diags_avgTstep
# define vidZ bulk_diags_avgZ
# define vidSustr_blk avgSustr_blk
# define vidSvstr_blk avgSvstr_blk
# define vidShflx_net avgShflx_net
# define vidShflx_lat avgShflx_lat
# define vidShflx_sen avgShflx_sen
# define vidShflx_rad avgShflx_rad
# define vidSwflx_emp avgSwflx_emp
# define vidShflx_wwk avgShflx_wwk
# define vidSurf_u avgSurf_u
# define vidSurf_v avgSurf_v

#endif /* AVRH */
      character*60 text

! This check is not performed in def_his.F, but very useful:
! If the output frequency is not a positive number, do nothing at all
#ifndef AVRH
      if (n_bulk_diags_his < 1) return
#else
      if (n_bulk_diags_avg < 1) return
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
      lstr=lenstr(ncname)
      if (rec_per_file.gt.0) then
        ierr=0
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
! file is prepared for appending new data; if it fails, a new file
! is created.
!
#ifndef AVRH
      create_new_file = new_bulk_diags_his
#else
      create_new_file = new_bulk_diags_avg
#endif
      if (ncid.ne.-1) create_new_file=.false.
#if defined MPI && !defined PARALLEL_FILES
      if (mynode.gt.0) create_new_file=.false.
#endif
!
! Create new history/averages file:    Put global attributes
! ====== === ======= ======== =====    and define all variables.
!
  10  create_new_file_if: if (create_new_file) then
        ierr=nf_create(ncname(1:lstr), nf_clobber, ncid)
        if (ierr .ne. nf_noerr) then
          write(*,'(/3(1x,A)/)') 'ERROR in def_bdiags_his/avg:',
     &           'Cannot create netCDF file:', 
     &          ncname(1:lstr)
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
        ierr=nf_def_dim (ncid, 'xi_rho',   xi_rho,  r2dgrd(1))
        ierr=nf_def_dim (ncid, 'xi_u',     xi_u,     u2dgrd(1))
        ierr=nf_def_dim (ncid, 'eta_rho',  eta_rho,  r2dgrd(2))
        ierr=nf_def_dim (ncid, 'eta_v',    eta_v,    v2dgrd(2))
        ierr=nf_def_dim (ncid, 'time', nf_unlimited, timedim)
        ierr=nf_def_dim (ncid, 'auxil',    4,        auxil(1))
        auxil(2)=timedim

        r2dgrd(3)=timedim           ! Free surface

        u2dgrd(2)=r2dgrd(2)         ! 2D UBAR-type
        u2dgrd(3)=timedim

        v2dgrd(1)=r2dgrd(1)         ! 2D VBAR-type
        v2dgrd(3)=timedim

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
     &       vidTstep)
        ierr=nf_put_att_text (ncid, vidTstep, 'long_name', 48,
     &       'time step and record numbers from initialization')
!
! Time.
!
        lvar=lenstr(vname(1,indxTime))
        ierr=nf_def_var (ncid, vname(1,indxTime)(1:lvar),
     &       NF_DOUBLE, 1, timedim, vidTime)
#ifndef AVRH
        lvar=lenstr(vname(2,indxTime))
        ierr=nf_put_att_text (ncid, vidTime, 'long_name', 
     &       lvar, vname(2,indxTime)(1:lvar))
#else
        text='averaged '/ /vname(2,indxTime)
        lvar=lenstr(text)
        ierr=nf_put_att_text (ncid, vidTime, 'long_name', 
     &       lvar, text(1:lvar))
#endif
        lvar=lenstr(vname(3,indxTime))
        ierr=nf_put_att_text (ncid, vidTime, 'units',  lvar,
     &       vname(3,indxTime)(1:lvar))

! Free surface
        lvar=lenstr(vname(1,indxZ))
        ierr=nf_def_var (ncid, vname(1,indxZ)(1:lvar),
     &       NF_FOUT, 3, r2dgrd, vidZ)
#ifndef AVRH
          lvar=lenstr(vname(2,indxZ))
          ierr=nf_put_att_text (ncid, vidZ, 'long_name', lvar,
     &                                  vname(2,indxZ)(1:lvar))
#else
          text='averaged '/ /vname(2,indxZ)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidZ, 'long_name', lvar,
     &         text(1:lvar))
#endif
        lvar=lenstr(vname(3,indxZ))
        ierr=nf_put_att_text (ncid, vidZ, 'units',     lvar,
     &       vname(3,indxZ)(1:lvar))
#ifdef MASK_LAND_DATA
        ierr=nf_put_att_FOUT (ncid, vidZ, '_FillValue',
     &       NF_FOUT, 1, spv_set)
#endif



! u-stress
!
#ifdef WND_AT_RHO_POINTS
        ierr=nf_def_var (ncid, 'sustr',
     &   NF_FOUT,3, r2dgrd, vidSustr_blk)
#else
        ierr=nf_def_var (ncid, 'sustr',
     &   NF_FOUT,3, u2dgrd, vidSustr_blk)
#endif         
#ifdef AVRH
         text='averaged zonal wind stress'
#else
          text='zonal wind stress'
#endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidSustr_blk, 'long_name', lvar,
     &                                                   text(1:lvar))
          text='N m^-2'
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidSustr_blk, 'units',     lvar,
     &                                                   text(1:lvar))
#ifdef MASK_LAND_DATA 
          ierr=nf_put_att_FOUT (ncid, vidSustr_blk, '_FillValue',
     &                                   NF_FOUT, 1, spv_set)
#endif
!
! v-stress
!
# ifdef WND_AT_RHO_POINTS
        ierr=nf_def_var (ncid, 'svstr',
     &   NF_FOUT,3, r2dgrd, vidSvstr_blk)
# else
        ierr=nf_def_var (ncid, 'svstr',
     &   NF_FOUT,3, v2dgrd, vidSvstr_blk)
# endif
# ifdef AVRH
          text='averaged meridional wind stress'
# else
          text='meridional wind stress'
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidSvstr_blk, 'long_name', lvar,
     &                                                   text(1:lvar))
          text='N m^-2'
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidSvstr_blk, 'units',     lvar,
     &                                                   text(1:lvar))
# ifdef MASK_LAND_DATA 
          ierr=nf_put_att_FOUT (ncid, vidSvstr_blk, '_FillValue',
     &                                   NF_FOUT, 1, spv_set)
# endif
!
! Net SHFLX
!
        ierr=nf_def_var (ncid, 'netheat',
     &   NF_FOUT,3, r2dgrd, vidShflx_net)
# ifdef AVRH
          text='averaged net heat flux'
# else
          text='net heat flux'
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidShflx_net, 'long_name', lvar,
     &                                              text(1:lvar))
          text='W m^-2'
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidShflx_net, 'units',     lvar,
     &                                              text(1:lvar))
# ifdef MASK_LAND_DATA 
          ierr=nf_put_att_FOUT (ncid, vidShflx_net, '_FillValue',
     &                                   NF_FOUT, 1, spv_set)
# endif

!
! Latent heat flux
!
        ierr=nf_def_var (ncid, 'latent',
     &   NF_FOUT,3, r2dgrd, vidShflx_lat)
# ifdef AVRH
          text='averaged latent heat flux'
# else
          text='latent heat flux'
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidShflx_lat, 'long_name', lvar,
     &                                              text(1:lvar))
          text='W m^-2'
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidShflx_lat, 'units',     lvar,
     &                                              text(1:lvar))
# ifdef MASK_LAND_DATA 
          ierr=nf_put_att_FOUT (ncid, vidShflx_lat, '_FillValue',
     &                                   NF_FOUT, 1, spv_set)
# endif

!
! Sensible heat flux
!
        ierr=nf_def_var (ncid, 'sensible',
     &   NF_FOUT,3, r2dgrd, vidShflx_sen)
#  ifdef AVRH
          text='averaged sensible heat flux'
#  else
          text='sensible heat flux'
#  endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidShflx_sen, 'long_name', lvar,
     &                                              text(1:lvar))
          text='W m^-2'
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidShflx_sen, 'units',     lvar,
     &                                              text(1:lvar))
#ifdef MASK_LAND_DATA 
          ierr=nf_put_att_FOUT (ncid, vidShflx_sen, '_FillValue',
     &                                   NF_FOUT, 1, spv_set)
#endif



!
! Short wave SHFLX
!
        ierr=nf_def_var (ncid, 'swrad',
     &   NF_FOUT,3, r2dgrd, vidShflx_rad)
#  ifdef AVRH
          text='averaged net shortwave radiation'
#  else
          text='net shortwave radiation'
#  endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidShflx_rad, 'long_name', lvar,
     &                                              text(1:lvar))
          text='W m^-2'
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidShflx_rad, 'units',     lvar,
     &                                              text(1:lvar))
#ifdef MASK_LAND_DATA 
          ierr=nf_put_att_FOUT (ncid, vidShflx_rad, '_FillValue',
     &                                   NF_FOUT, 1, spv_set)
#endif
!
! Freshwater flux
!
        ierr=nf_def_var (ncid, 'freshwater',
     &   NF_FOUT,3, r2dgrd, vidSwflx_emp)
#  ifdef AVRH
          text='averaged surface freshwater flux (E-P)'
#  else
          text='surface freshwater flux (E-P)'
#  endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidSwflx_emp, 'long_name', lvar,
     &                                                   text(1:lvar))
          text='cm day^-1'
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidSwflx_emp, 'units',     lvar,
     &                                                   text(1:lvar))
# ifdef MASK_LAND_DATA 
          ierr=nf_put_att_FOUT (ncid, vidSwflx_emp, '_FillValue',
     &                                   NF_FOUT, 1, spv_set)
# endif

!
! Wind Work
!
        ierr=nf_def_var (ncid, 'windwork',
     &   NF_FOUT,3, r2dgrd, vidShflx_wwk)
# ifdef AVRH
          text='averaged wind work'
# else
          text='wind work'
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidShflx_wwk, 'long_name', lvar,
     &                                              text(1:lvar))
          text=''
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidShflx_wwk, 'units',     lvar,
     &                                              text(1:lvar))
# ifdef MASK_LAND_DATA 
          ierr=nf_put_att_FOUT (ncid, vidShflx_wwk, '_FillValue',
     &                                   NF_FOUT, 1, spv_set)
# endif


!
! Surface U and V
!
        ierr=nf_def_var (ncid, 'surf_u',
     &   NF_FOUT,3, r2dgrd, vidSurf_u)
# ifdef AVRH
          text='averaged surface u'
# else
          text='surface u'
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidSurf_u, 'long_name', lvar,
     &                                              text(1:lvar))
          text=''
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidSurf_u, 'units',     lvar,
     &                                              text(1:lvar))
# ifdef MASK_LAND_DATA 
          ierr=nf_put_att_FOUT (ncid, vidSurf_u, '_FillValue',
     &                                   NF_FOUT, 1, spv_set)
# endif

! Surface U and V
!
        ierr=nf_def_var (ncid, 'surf_v',
     &   NF_FOUT,3, r2dgrd, vidSurf_v)
# ifdef AVRH
          text='averaged surface v'
# else
          text='surface v'
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidSurf_v, 'long_name', lvar,
     &                                              text(1:lvar))
          text=''
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidSurf_v, 'units',     lvar,
     &                                              text(1:lvar))
# ifdef MASK_LAND_DATA 
          ierr=nf_put_att_FOUT (ncid, vidSurf_v, '_FillValue',
     &                                   NF_FOUT, 1, spv_set)
# endif



!
! Leave definition mode.
! ----- ---------- -----
!
        ierr=nf_enddef(ncid)
        mpi_master_only write(*,'(6x,4A,1x,A,i4)')
#ifdef AVRH
     &       'DEF_BDIAGS_AVG - Created ',
#else
     &       'DEF_BDIAGS_HIS - Created ',
#endif /* AVRH */
     &                'new netCDF file ''', 
     &       ncname(1:lstr), '''.'
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
      elseif (ncid.eq.-1) then create_new_file_if
        ierr=nf_open (ncname(1:lstr), nf_write, ncid)
        if (ierr. eq. nf_noerr) then
          mpi_master_only write(*,'(1x,2A,1x,A,I3)') 
     &          'Opened existing bdiags_his file ',
     &          ncname(1:lstr), 'ncid =', ncid
          ierr=checkdims (ncid, ncname(1:lstr), rec)
          if (ierr .eq. nf_noerr) then
            if (rec_per_file.eq.0) then
              ierr=rec+1 - total_rec
            else
              ierr=rec+1 - (1+mod(total_rec-1, rec_per_file))
            endif
            if (ierr.gt.0) then
              mpi_master_only write( *,
     &                 '(/1x,A,I5,1x,A/8x,3A,I5,/8x,A,I5,1x,A/)')
#ifdef AVRH
     &        'DEF_BDIAGS_AVG WARNING: Actual number of records',
#else
     &        'DEF_BDIAGS_HIS WARNING: Actual number of records',
#endif /* AVRH */
     &               rec,  'in netCDF file',  '''',  
     &              ncname(1:lstr),
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
            write(*,'(/1x,4A,2x,A,I4/)')
#ifdef AVRH
     &            'DEF_BDIAGS_AVG ERROR: ',
#else
     &            'DEF_BDIAGS_HIS ERROR: ',
#endif /* AVRH */
     &                  'Cannot open file ''', 
     &            ncname(1:lstr), '''.'
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
        ierr=nf_inq_varid (ncid,vname(1,indxTime)(1:lvar),
     &       vidTime)
        if (ierr .ne. nf_noerr) then
          write(*,1) vname(1,indxTime)(1:lvar), 
     &          ncname(1:lstr)
          goto 99                                         !--> ERROR 
        endif
! Free surface
        lvar=lenstr(vname(1,indxZ))
        ierr=nf_inq_varid (ncid,vname(1,indxZ)(1:lvar),
     &       vidZ)
        if (ierr .ne. nf_noerr) then
          write(*,1) vname(1,indxZ)(1:lvar), 
     &          ncname(1:lstr)
          goto 99                                         !--> ERROR 
        endif
!


        ierr=nf_inq_varid ( ncid,'sustr'     , vidSustr_blk )
        ierr=nf_inq_varid ( ncid,'svstr'     , vidSvstr_blk )
        ierr=nf_inq_varid ( ncid,'netheat'   , vidShflx_net )
        ierr=nf_inq_varid ( ncid,'swrad'     , vidShflx_rad )
        ierr=nf_inq_varid ( ncid,'freshwater', vidSwflx_emp )
        ierr=nf_inq_varid ( ncid,'windwork'  , vidShflx_wwk )


        mpi_master_only write(*,'(6x,2A,i4,1x,A,i4)')
#ifdef AVRH
     &       'DEF_BDIAGS_AVG -- Opened ',
#else
     &       'DEF_BDIAGS_HIS -- Opened ',
#endif /* AVRH */
     &                     'existing file  from record =', rec
     &                      MYID

#if defined MPI && !defined PARALLEL_FILES
      else
        ierr=nf_open (ncname(1:lstr), nf_write, ncid)
        if (ierr .ne. nf_noerr) then
#ifdef AVRH
          write(*,'(/1x,4A,2x,A,I4/)') 'DEF_BDIAGS_AVG ERROR: ',
#else
          write(*,'(/1x,4A,2x,A,I4/)') 'DEF_BDIAGS_HIS ERROR: ',
#endif /* AVRH */
     &                'Cannot open file ''', ncname(1:lstr), '''.'
     &                 MYID
          goto 99                                         !--> ERROR
        endif
#endif
      end if create_new_file_if
      
      ierr=nf_set_fill (ncid, nf_nofill, lvar)
      if (ierr .ne. nf_noerr) then
        write(*,'(6x,2A,i4,1x,A,i4)') 
#ifdef AVRH
     &        'DEF_BDIAGS_AVG ERROR: Cannot ', 
#else
     &        'DEF_BDIAGS_HIS ERROR: Cannot ', 
#endif /* AVRH */
     &    'switch to ''nf_nofill'' more; netCDF error code =', ierr 
      endif
#ifdef AVRH
   1  format(/1x,'DEF_BDIAGS_AVG ERROR: Cannot find variable ''',
#else
   1  format(/1x,'DEF_BDIAGS_HIS ERROR: Cannot find variable ''',
#endif /* AVRH */
     &                   A, ''' in netCDF file ''', A, '''.'/)
#if (defined PUT_GRID_INTO_HISTORY && !defined AVRH)\
 || (defined PUT_GRID_INTO_AVERAGES && defined AVRH)
!
! Write grid variables.
! ----- ---- ----------
!
        if (total_rec.le.1) call wrt_grid (ncid, 
     &      ncname, lstr)
#endif
  99  return
      end

#undef ncname
#undef rec_per_file
#undef vidTime
#undef vidTstep
#undef vidZ
#undef vidSustr_blk
#undef vidSvstr_blk
#undef vidShflx_net
#undef vidShflx_rad
#undef vidSwflx_emp
#undef vidShflx_wwk
#undef vidShflx_lat
#undef vidShflx_sen
#undef vidSurf_u
#undef vidSurf_v

# ifndef AVRH
#  ifdef AVERAGES
#   define AVRH
#   include "def_bulk_diags.F"
#  endif /* AVERAGES */
#  endif
#else /* BDIAGS_ANALYSIS */
      subroutine def_bulk_diags_empty()
      return
      end
#endif /* BDIAGS_ANALYSIS */
