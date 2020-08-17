#include "cppdefs.h"
 
      integer function nf_read_bry_EW (A, ncid, varid, record, type)
 
! Read side boundary from an input netCDF file.
!
! Arguments: A       real array of standard horizontal dimensions
!                             which is to be read or written.
!            ncid    NetCDF ID of in the file.
!            varid   variable ID of that variable in NetCDF file.
!            record  record number.
!            type    type of the grid (RHO-, U, V, W, PSI etc.)
 
      implicit none
#include "param.h"
      real A(0:Mm+1,N)
      integer ncid, varid, record, type,  i,j,k, kmax, shft,
     &                                    ierr, lstr, lenstr
      character(len=16) vname
#include "netcdf.inc"
#include "buffer.h"
#if defined MPI
# include "hidden_mpi_vars.h"
# ifdef HOFFMAN2
      include "mpif.h"
# else
# include "mpif.h"
# endif
      integer, parameter :: ibss=100000
      integer req(2), status(MPI_STATUS_SIZE,2)
      real bffr_rsv_S(N+1), bffr_snd_S(N+1),
     &     bffr_rsv_N(N+1), bffr_snd_N(N+1)
# endif
#include "compute_starts_counts.h"
 
#ifdef VERBOSE
      write(*,*) 'entered nf_read_bry_EW...' MYID
#endif
                              ! Note that "ierr" is initialised
      if (ierr.eq.0) then     ! inside "compute_starts_counts.h".
        start(1)=start(2)     ! Suppress XI- the fist horizontal
        count(1)=count(2)     ! dimension: after this step updated
        start(2)=start(3)     ! start/count(1) corresponds to ETA-
        count(2)=count(3)     ! netCDF dimension; index (2) to
        start(3)=start(4)     ! vertical dimension (if any) or time
        count(3)=count(4)     ! record; (3) (if any) to time record.
 
! Read array from the disk.
!===== ===== ==== === =====
 
        ierr=nf_get_vara_FTYPE (ncid, varid, start, count, buff)
        if (ierr .eq. nf_noerr) then
          kmax=count(2)
          do k=1,kmax
            shft=1-jmin + (k-1)*count(1)
            do j=jmin,jmax
              A(j,k)=buff(j+shft)
            enddo
          enddo
#ifdef NS_PERIODIC
# ifdef MPI
          if (NP_ETA.eq.1) then
# endif
            do k=1,kmax            ! note that this periodic BC code
              A(Mm+1,k)=A( 1,k)    ! is never called from inside Open
              A(   0,k)=A(Mm,k)    ! MP parllel region associated
            enddo                  ! with tiling.
# ifdef MPI
          else
# endif
#endif
#ifdef MPI
# ifdef VERBOSE
          write(*,*) 'starting nf_read_bry_EW mpi exchange part' MYID
# endif
            if (south_msg_exch) then
              call MPI_Irecv (bffr_rsv_S, kmax, MPI_DOUBLE_PRECISION,
     &                    p_S, ibss+1, ocean_grid_comm, req(1), ierr)
            endif
            if (north_msg_exch) then
              call MPI_Irecv (bffr_rsv_N, kmax, MPI_DOUBLE_PRECISION,
     &                    p_N, ibss+2, ocean_grid_comm, req(2), ierr)
            endif
            if (south_msg_exch) then
              do k=1,kmax
                bffr_snd_S(k)=A(jmin,k)
              enddo
              call MPI_Send  (bffr_snd_S, kmax, MPI_DOUBLE_PRECISION,
     &                     p_S, ibss+2, ocean_grid_comm,        ierr)
            endif
            if (north_msg_exch) then
              do k=1,kmax
                bffr_snd_N(k)=A(jmax,k)
              enddo
              call MPI_Send  (bffr_snd_N, kmax, MPI_DOUBLE_PRECISION,
     &                     p_N, ibss+1, ocean_grid_comm,        ierr)
            endif
            if (south_msg_exch) then
              call MPI_Wait (req(1),status(1,1),ierr)
              do k=1,kmax
                A(jmin-1,k)=bffr_rsv_S(k)
              enddo
            endif
            if (north_msg_exch) then
              call MPI_Wait (req(2),status(1,2),ierr)
              do k=1,kmax
                A(jmax+1,k)=bffr_rsv_N(k)
              enddo
            endif
# ifdef NS_PERIODIC
          endif ! <-- NP_ETA.eq.1
# endif
#endif
        else
          ierr=nf_inq_varname (ncid, varid, vname)
          if (ierr.eq.nf_noerr) then
            lstr=lenstr(vname)
            write(*,'(1x,2A,I6,1x,3A,2x,A,2x,A,I4)')  '### ERROR: ',
     &              'nf_read_bry_EW ::Cannot read rec =',   record,
     &              'of variable ''',  vname(1:lstr),  '''.',
     &               nf_strerror(ierr) MYID
          else
            write(*,'(1x,2A,I6,1x,A,I4,2x,A,2x,A,I4)') '### ERROR: ',
     &              'nf_read_bry_EW ::Cannot read rec =',   record,
     &              'of variable with netCDF ID =',  varid,
     &               nf_strerror(ierr) MYID
          endif
        endif
 
      else
        write(*,'(/1x,2A,I4/)') '### ERROR: nf_read_bry_EW :: ',
     &                          'illegal grid type', type
      endif
      nf_read_bry_EW=ierr
#ifdef VERBOSE
      write(*,*) 'return from nf_read_bry_EW'  MYID
#endif
      return
      end
 
 
 
      integer function nf_read_bry_NS (A, ncid, varid, record, type)
      implicit none
#include "param.h"
      real A(0:Lm+1,N)
      integer ncid, varid, record, type,  i, k, kmax, shft,
     &                                    ierr, lstr, lenstr
      character(len=16) vname
#include "netcdf.inc"
#include "buffer.h"
#if defined MPI
# include "hidden_mpi_vars.h"
# ifdef HOFFMAN2
      include "mpif.h"
# else
# include "mpif.h"
# endif
      integer, parameter :: ibss=100000
      integer req(2), status(MPI_STATUS_SIZE,2)
      real bffr_rsv_W(N+1), bffr_rsv_E(N+1),
     &     bffr_snd_W(N+1), bffr_snd_E(N+1)
#endif
#include "compute_starts_counts.h"
 
#ifdef VERBOSE
      write(*,*) 'entered nf_read_bry_NS...'   MYID
#endif
      if (ierr.eq.0) then     ! Suppress ETA-netCDF dimension. After
        start(2)=start(3)     ! this change start/count(1) corresponds
        count(2)=count(3)     ! to XI-dimensions (as before); index
        start(3)=start(4)     ! (2) to vertical dimension (if any) or
        count(3)=count(4)     ! time record; (3) (if any) to time
                              ! record.
! Read array from the disk.
!===== ===== ==== === =====
 
        ierr=nf_get_vara_FTYPE (ncid, varid, start, count, buff)
        if (ierr .eq. nf_noerr) then
          kmax=count(2)
          do k=1,kmax
            shft=1-imin + (k-1)*count(1)
            do i=imin,imax
              A(i,k)=buff(i+shft)
            enddo
          enddo
#ifdef EW_PERIODIC
# ifdef MPI
          if (NP_XI.eq.1) then
# endif
            do k=1,kmax            ! note that this periodic BC code
              A(Lm+1,k)=A( 1,k)    ! is never called from inside Open
              A(   0,k)=A(Lm,k)    ! MP parllel region associated
            enddo                  ! with tiling.
# ifdef MPI
          else
# endif
#endif
#ifdef MPI
# ifdef VERBOSE
          write(*,*) 'starting nf_read_bry_NS mpi exchange part' MYID
# endif
            if (west_msg_exch) then
              call MPI_Irecv (bffr_rsv_W, kmax, MPI_DOUBLE_PRECISION,
     &                    p_W, ibss+3, ocean_grid_comm, req(1), ierr)
            endif
            if (east_msg_exch) then
              call MPI_Irecv (bffr_rsv_E, kmax, MPI_DOUBLE_PRECISION,
     &                    p_E, ibss+4, ocean_grid_comm, req(2), ierr)
            endif
            if (west_msg_exch) then
              do k=1,kmax
                bffr_snd_W(k)=A(imin,k)
              enddo
              call MPI_Send  (bffr_snd_W, kmax, MPI_DOUBLE_PRECISION,
     &                      p_W, ibss+4, ocean_grid_comm,       ierr)
            endif
            if (east_msg_exch) then
              do k=1,kmax
                bffr_snd_E(k)=A(imax,k)
              enddo
              call MPI_Send  (bffr_snd_E, kmax, MPI_DOUBLE_PRECISION,
     &                      p_E, ibss+3, ocean_grid_comm,       ierr)
            endif
            if (west_msg_exch) then
              call MPI_Wait (req(1),status(1,1),ierr)
              do k=1,kmax
                A(imin-1,k)=bffr_rsv_W(k)
              enddo
            endif
            if (east_msg_exch) then
              call MPI_Wait (req(2),status(1,2),ierr)
              do k=1,kmax
                A(imax+1,k)=bffr_rsv_E(k)
              enddo
            endif
# ifdef EW_PERIODIC
          endif !<-- NP_XI.eq.1
# endif
#endif
        else
          ierr=nf_inq_varname (ncid, varid, vname)
          if (ierr.eq.nf_noerr) then
            lstr=lenstr(vname)
            write(*,'(1x,2A,I6,1x,3A,2x,A,2x,A,I4)')  '### ERROR: ',
     &              'nf_read_bry_NS ::Cannot read rec =',   record,
     &              'of variable ''',  vname(1:lstr),  '''.',
     &               nf_strerror(ierr) MYID
          else
            write(*,'(1x,2A,I6,1x,A,I4,2x,A,2x,A,I4)') '### ERROR: ',
     &              'nf_read_bry_NS ::Cannot read rec =',   record,
     &              'of variable with netCDF ID =',  varid,
     &               nf_strerror(ierr) MYID
          endif
        endif
 
      else
        write(*,'(/1x,2A,I4/)') '### ERROR: nf_read_bry_NS :: ',
     &                          'illegal grid type', type
      endif
      nf_read_bry_NS=ierr
#ifdef VERBOSE
      write(*,*) 'return from nf_read_bry_NS' MYID
#endif
      return
      end
