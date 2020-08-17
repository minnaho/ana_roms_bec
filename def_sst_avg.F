!
!
# include "cppdefs.h"
#if defined SOLVE3D && defined OUTPUT_SST_AVG
      subroutine def_sst_avg (ncid, total_rec, ierr)
!
! Create/open SST averages netCDF file. In the case when a new
! netCDF file is created, define all variables, their dimensions and
! attributes. In the case when a previously existing netCDF file is
! to be opened for addition of new data, verify that all dimensions
! of the file are consistent with the present model configuration
! and all necessary variables exist. Save netCDF IDs for all needed
! variables. Also determine size of the unlimited dimension.
!
      implicit none
      logical create_new_file
      integer ncid, total_rec, ierr, rec, lstr,lvar,lenstr, timedim
     &      , r2dgrd(3),  auxil(2), checkdims
     &      , r3dgrd(4),  itrc
#include "param.h"
#include "scalars.h"
#include "ncvars.h"
#include "netcdf.inc"

      character*70 text
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
      lstr=lenstr(sst_avgname)
      if (nrpf_sst_avg.gt.0) then
        lvar=total_rec-(1+mod(total_rec-1, nrpf_sst_avg))
        call insert_time_index (sst_avgname, lstr, lvar, ierr)
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
!!HF      create_new_file=ldefhis

! hf: default - create new file
      create_new_file=.true.
      if (ncid.ne.-1) create_new_file=.false.
#if defined MPI & !defined PARALLEL_FILES
      if (mynode.gt.0) create_new_file=.false.
#endif
!
! Create new SST averages file:    Put global attributes
! ====== === === ======== =====    and define all variables.
!
  10  if (create_new_file) then
        ierr=nf_create(sst_avgname(1:lstr), nf_clobber, ncid)
        if (ierr .ne. nf_noerr) then
          write(*,'(/3(1x,A)/)') 'ERROR in def_sst_avg:',
     &           'Cannot create netCDF file:', sst_avgname(1:lstr)
          goto 99                                         !--> ERROR
        endif
        if (nrpf_sst_avg.eq.0) total_rec=0  
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
        ierr=nf_def_dim (ncid, 'eta_rho',  eta_rho,  r2dgrd(2))
        ierr=nf_def_dim (ncid, 'time', nf_unlimited, timedim)
        ierr=nf_def_dim (ncid, 'auxil',    4,        auxil(1))
        
        auxil(2)=timedim

        r2dgrd(3)=timedim           ! Free surface

#ifdef PUT_GRID_INTO_AVERAGES
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
     &       avg_sst_Tstep)
        ierr=nf_put_att_text (ncid, avg_sst_Tstep, 'long_name', 48,
     &       'time step and record numbers from initialization')
!
! Time.
!
        lvar=lenstr(vname(1,indxTime))
! hf: always use double precision for the time variable
        ierr=nf_def_var (ncid, vname(1,indxTime)(1:lvar),
     &                            NF_DOUBLE, 1, timedim, avg_sst_time)
        text='Averaged '/ /vname(2,indxTime)
        lvar=lenstr(text)
        ierr=nf_put_att_text (ncid, avg_sst_time, 'long_name', lvar,
     &                                             text(1:lvar))
        lvar=lenstr(vname(3,indxTime))
        ierr=nf_put_att_text (ncid, avg_sst_time, 'units',  lvar,
     &                                vname(3,indxTime)(1:lvar))
!
! SST
!
        text='SST'
        lvar=lenstr(text)
        ierr=nf_def_var (ncid, text(1:lvar),
     &       NF_FOUT, 3, r2dgrd, avgSST)
        text='Averaged SST'
        lvar=lenstr(text)
        ierr=nf_put_att_text (ncid, avgSST, 'long_name',
     &       lvar, text(1:lvar))
        lvar=lenstr(vname(3,indxT))
        ierr=nf_put_att_text (ncid, avgSST, 'units', lvar,
     &       vname(3,indxT)(1:lvar))
!
! Leave definition mode.
! ----- ---------- -----
!
        ierr=nf_enddef(ncid)
        write(*,'(6x,4A,1x,A,i4)') 'DEF_SST_AVG - Created ',
     &                'new netCDF file ''', sst_avgname(1:lstr), '''.'
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
        ierr=nf_open (sst_avgname(1:lstr), nf_write, ncid)
        if (ierr. eq. nf_noerr) then
          ierr=checkdims (ncid, sst_avgname(1:lstr), rec) 
          if (ierr .eq. nf_noerr) then
            if (nrpf_sst_avg.eq.0) then
              ierr=rec+1 - total_rec
            else
              ierr=rec+1 - (1+mod(total_rec-1, nrpf_sst_avg))
            endif
            if (ierr.gt.0) then
              mpi_master_only write( *,
     &                 '(/1x,A,I5,1x,A/8x,3A,I5,/8x,A,I5,1x,A/)'
     &            ) 'DEF_SST_AVG WARNING: Actual number of records',
     &               rec,  'in netCDF file',  '''', sst_avgname(1:lstr),
     &             ''' exceeds the record number from restart data',
     &             rec+1-ierr,'/', total_rec,', restart is assumed.'
              rec=rec-ierr
            elseif (nrpf_sst_avg.eq.0) then
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
            write(*,'(/1x,4A,2x,A,I4/)') 'DEF_SST_AVG ERROR: ',
     &            'Cannot open file ''', sst_avgname(1:lstr), '''.'
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
        ierr=nf_inq_varid (ncid, 'time_step', avg_sst_Tstep)
        if (ierr .ne. nf_noerr) then
          write(*,1) 'time_step', sst_avgname(1:lstr)
          goto 99                                         !--> ERROR
        endif
!
! Time.
!
        lvar=lenstr(vname(1,indxTime))
        ierr=nf_inq_varid (ncid,vname(1,indxTime)(1:lvar),avg_sst_time)
        if (ierr .ne. nf_noerr) then
          write(*,1) vname(1,indxTime)(1:lvar), sst_avgname(1:lstr)
          goto 99                                         !--> ERROR 
        endif
!
! SST
!
        text='SST'
        lvar=lenstr(text)
        ierr=nf_inq_varid (ncid, text(1:lvar), avgSST)
        if (ierr .ne. nf_noerr) then
           write(*,1) text(1:lvar), sst_avgname(1:lstr) 
           goto 99              !--> ERROR
        endif
        write(*,'(6x,2A,i4,1x,A,i4)') 'DEF_SST_AVG -- Opened ',
     &                     'existing file  from record =', rec
     &                      MYID

#if defined MPI & !defined PARALLEL_FILES
      else
        ierr=nf_open (sst_avgname(1:lstr), nf_write, ncid)
        if (ierr .ne. nf_noerr) then
          write(*,'(/1x,4A,2x,A,I4/)') 'DEF_SST_AVG ERROR: ',
     &                'Cannot open file ''', sst_avgname(1:lstr), '''.'
     &                 MYID
          goto 99                                         !--> ERROR
        endif
#endif
      endif             !<--  create_new_file  
      
      ierr=nf_set_fill (ncid, nf_nofill, lvar)
      if (ierr .ne. nf_noerr) then
        write(*,'(6x,2A,i4,1x,A,i4)') 'DEF_SST_AVG ERROR: Cannot ', 
     &    'switch to ''nf_nofill'' more; netCDF error code =', ierr 
      endif
   1  format(/1x,'DEF_SST_AVG ERROR: Cannot find variable ''',
     &                   A, ''' in netCDF file ''', A, '''.'/)
#ifdef PUT_GRID_INTO_AVERAGES
!
! Write grid variables.
! ----- ---- ----------
!
        if (total_rec.le.1) call wrt_grid (ncid, sst_avgname, lstr)
#endif
!
! Set initial average time.
! --- ------- ------- -----
!
      if (nts_sst_avg.eq.1) then
        time_sst_avg=time-0.5*float(n_sst_avg)*dt
      else
        time_sst_avg=time-0.5*float(n_sst_avg)*dt+float(nts_sst_avg)*dt
      endif

  99  return
      end

#else /* SOLVE3D && OUTPUT_SST_AVG */
      subroutine def_sst_avg_empty
      return
      end
#endif  /* SOLVE3D && OUTPUT_SST_AVG */
