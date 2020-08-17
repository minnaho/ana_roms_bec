#include "cppdefs.h"
!#if defined PSOURCE || defined Q_PSOURCE && !defined ANA_PSOURCE
#ifndef ANA_PSOURCE  && defined PSOURCE

      subroutine get_psrcs (ierr)
!      subroutine get_psrcs (istr,iend,jstr,jend,ierr)

! Read in point sources of flow and tracer at the appropriate
! time from psource NetCDF files.  Written by Eileen Deng-Idica
! (05-03-2007) & modified by Yusuke Uchiyama (06-20-2011)
! modified by Faycal Kessouri (08-2018) : read point-sources of 
! pollution inputs and rivers from a NetCDF file
      implicit none
# include "param.h"
# include "forces.h"
# include "sources.h"
# include "scalars.h"
# include "netcdf.inc"
# include "ncvars.h"
# include "ocean3d.h"
# include "grid.h"
# ifdef MPI
#  include "hidden_mpi_vars.h"
# endif

      real cff, csum, cff1,cff2
      real Zcff, Scff
      integer ierr, lstr,lvar,lenstr,ifrc,Npas,nf_fread
      integer istr,iend,jstr,jend
      integer is, itrc, ncid, varid, i,j,k, it1,it2
#ifdef MPI
     &                        , i1,i2,j1,j2
#endif
      integer, parameter :: pdims=2
      integer            :: start(pdims), count(pdims)
      real,    dimension(:,:), allocatable :: D
      integer, dimension(:,:), allocatable :: C
      character(len=7)  snsrc, sk
      character(len=64) fmt
!      character*(max_name_size) psrc_file

! Initialization: Inquire about the contents of psource NetCDF file:
!================  variables and dimensions. Check for consistency.

      if (may_day_flag.ne.0) return      !-->  EXIT


      lstr=lenstr(psrc_file)
      ncid=ncidfpsrc
      ierr=nf_noerr

      if (iic.eq.ntstart) then

        if (ncid.eq.-1) then
          ierr=nf_open(psrc_file(1:lstr), nf_nowrite, ncid)
          if (ierr.eq.nf_noerr) then
            ncidfpsrc=ncid
          else
            write(*,'(/1x,4A,1x,2A,I4 /)') '### ERROR: get_psource',
     &         ' :: Cannot open netCDF file ''', psrc_file(1:lstr),
     &         '''.',nf_strerror(ierr) MYID
            goto 99                                      !--> ERROR
          endif
        endif

      write(*,*) 'source_file_name = ',psrc_file(1:lstr)

! time
        ierr=nf_inq_varid (ncid, 'psrc_time', psrc_tid)
        if (ierr .eq. nf_noerr) then
          call set_cycle (ncid, psrc_tid, ntpsrc, psrc_cycle,
     &                    psrc_ncycle, psrc_rec, ierr)
          itpsrc=1
!        do i=1:Stepsrc
        do i=1,ntpsrc
          psrc_time(i)=-1.E+20
        enddo
        else
          write(*,1) 'psrc_time', psrc_file(1:lstr)
        endif


! Nsrc: number of point sources
        ierr=nf_inq_dimid (ncid, 'Nsrc', varid)
        if (ierr .eq. nf_noerr) then
          ierr=nf_inq_dimlen (ncid, varid, Nsrc)
          if (ierr .ne. nf_noerr) then
            write(*,1) 'Nsrc', psrc_file(1:lstr)
            goto 99                                      !--> ERROR
          else
            if (Nsrc.gt.Msrc) then
              write(*,*) 'FATAL ERROR: get_psrcs :: Found Nsrc > Msrc'
              goto 99                                    !--> ERROR
            endif
          endif
        else
          write(*,1) 'Nsrc', psrc_file(1:lstr)
          goto 99                                        !--> ERROR
        endif

        write(snsrc,'(i7)') Nsrc
        if (.not.allocated(C)) allocate (C(Nsrc,NT))
        if (.not.allocated(D)) allocate (D(Nsrc,N))

             write(*,*) 'get_psrc:: Nsrc = ',Nsrc

! Npas: number of tracers
        ierr=nf_inq_dimid (ncid, 'Npas', varid)
        if (ierr .eq. nf_noerr) then
          ierr=nf_inq_dimlen (ncid, varid, Npas)
          if (ierr .ne. nf_noerr) then
            write(*,1) 'Npas', psrc_file(1:lstr)
            goto 99                                      !--> ERROR
          else
            if (Npas.gt.NT) then
              write(*,*) 'FATAL ERROR: get_psrcs :: Found Npas > NT'
              goto 99                                    !--> ERROR
            endif
          endif
        else
          write(*,1) 'Npas', psrc_file(1:lstr)
          goto 99                                        !--> ERROR
        endif

! Isrc & Jsrc (Nsrc): point source locations relative to the global grid
        ierr=nf_inq_varid (ncid, 'Isrc', varid)
        if (ierr .eq. nf_noerr) then
          ierr=nf_get_var_int (ncid, varid, Isrc)
          if (ierr .ne. nf_noerr) then
            write(*,1) 'Isrc', psrc_file(1:lstr)
            goto 99                                      !--> ERROR
          endif
        else
          write(*,1) 'Isrc id', psrc_file(1:lstr)
          goto 99                                        !--> ERROR
        endif

        ierr=nf_inq_varid (ncid, 'Jsrc', varid)
        if (ierr .eq. nf_noerr) then
          ierr=nf_get_var_int (ncid, varid, Jsrc)
          if (ierr .ne. nf_noerr) then
            write(*,1) 'Jsrc', psrc_file(1:lstr)
            goto 99                                      !--> ERROR
          endif
        else
          write(*,1) 'Jsrc id', psrc_file(1:lstr)
          goto 99                                        !--> ERROR
        endif

! Initially indices "Isrc,Jsrc" are assumed to be defined relative
! to the absolute physical grid. Translate them into would be local
! indices within the current MPI-subdomain and check whether they
! fall within. If not, set unrealistic negative value to ignore
! these particular point sources.
        do is=1,Nsrc
            do itrc=1,NT
              Lsrc(is,itrc)=.true.
            enddo
        enddo

        do is=1,Nsrc
!        is = 30
#  ifdef MPI

          if ( (Isrc(is)-iSW_corn).lt.iwest  .or. (Isrc(is)-iSW_corn).gt.ieast  .or.
     &         (Jsrc(is)-jSW_corn).lt.jsouth .or. (Jsrc(is)-jSW_corn).gt.jnorth ) then

            do itrc=1,NT
              Lsrc(is,itrc)=.false.
            enddo

          else


            write(*,'(2X,A,X,I4,X,A,X,I3,A,2(I4,A))')
     &         '--> PSOURCE No.',is,'is found in node',mynode,
     &         ' : (i,j)=(',Isrc(is),',',Jsrc(is),')'
          endif
# endif
        enddo ! is

! Dsrc (Nsrc): direction of mass flux (0:xi-, 1:eta-, 2:s-direction)
        ierr=nf_inq_varid (ncid, 'Dsrc', varid)
        if (ierr .eq. nf_noerr) then
          ierr=nf_get_var_int (ncid, varid, Dsrc)
          if (ierr .ne. nf_noerr) then
            write(*,1) 'Dsrc', psrc_file(1:lstr)
            goto 99                                      !--> ERROR
          endif
        else
          write(*,1) 'Dsrc id', psrc_file(1:lstr)
          goto 99                                        !--> ERROR
        endif

        Dsrc(1) = 2
        Dsrc(2) = 2

! Qshape (Nsrc,N)
! from the file
        ierr=nf_inq_varid (ncid, 'Qshape', varid)
        if (ierr .eq. nf_noerr) then
          ierr=nf_get_var_FTYPE (ncid, varid, D)

! recalcaculate it online
!        write(*,*) 'D1 = '  , D(1,:)

            do is=1,Nsrc
            Scff=0.05
        cff=0.0
        if (is.le.28) then ! HYPERION
                Scff=0.02
                Zcff = abs(z_r(Isrc(is),Jsrc(is),0)) ! this is bottom,
                                                     ! orignal is 37.0
        elseif (is.le.56 .and. is.ge.29) then ! JWPCP
                Scff=0.02
                Zcff = 37.0
        elseif (is.le.70 .and. is.ge.57) then ! OCSD
                Scff=0.02
                Zcff = 37.0
        elseif (is.le.96 .and. is.ge.71) then  !POINT LOMA
                Scff=0.02
                Zcff = 65.0
        elseif (is.le.110 .and. is.ge.72) then ! small POTW
                Scff=0.02
                Zcff = max(2., abs(z_r(Isrc(is),Jsrc(is),0)) -15.0)
        elseif (is.ge.111) then  ! rivers
                Scff=0.02
                Zcff = max(2., abs(z_r(Isrc(is),Jsrc(is),59)))
        endif
                      do k=1,N
                        D(is,k)=exp(-Scff*
     &                  ( (z_r(Isrc(is),Jsrc(is),k)+Zcff)**2) )
                        cff=cff+D(is,k)
                      enddo ! k
          cff=1.0/cff
          do k=1,N
            D(is,k)=D(is,k)*cff
          enddo
            enddo ! is


          if (ierr .eq. nf_noerr) then
            Qshape=0.D0
            do is=1,Nsrc
!             is =30
              cff=0.D0; csum=0.D0
              do k=1,N
                Qshape(is,k)=D(is,k)
                cff=cff+Qshape(is,k)
              enddo
              if (cff.ne.0.D0) then
                csum=1./cff
              else
                csum=0.D0
              endif
              if (abs(cff-1.0).gt.1.D-6) then
                write(*,*) 'WARNING : get_psrcs :: Qshape is not ',
     &                     'normalized.  Correction has been made.'
                do k=1,N            ! normalize Qshape
                  Qshape(is,k)=Qshape(is,k)*csum
                enddo
              endif
            enddo   ! is
          else
            write(*,1) 'Qshape', psrc_file(1:lstr)
            goto 99                                      !--> ERROR
          endif
        else
          write(*,1) 'Qshape', psrc_file(1:lstr)
          goto 99                                        !--> ERROR
        endif

!        write(*,*) 'D2_2 = ' , D(1,:)

! Qbar id
        ierr=nf_inq_varid (ncid, 'Qbar', QbarG_id)
        if (ierr .ne. nf_noerr) then
          write(*,1) 'Qbar', psrc_file(1:lstr)
          goto 99                                        !--> ERROR
        endif


! var id
          ierr=nf_inq_varid (ncid, 'temp',tempG_id)
          if (ierr .ne. nf_noerr) then
            write(*,1) 'temp'
          endif

          ierr=nf_inq_varid (ncid, 'salt',saltG_id)
          if (ierr .ne. nf_noerr) then
            write(*,1) 'salt'
          endif

          ierr=nf_inq_varid (ncid, 'PO4',PO4G_id)
          if (ierr .ne. nf_noerr) then
            write(*,1) 'PO4'
          endif

          ierr=nf_inq_varid (ncid, 'NO3',NO3G_id)
          if (ierr .ne. nf_noerr) then
            write(*,1) 'NO3'
          endif

          ierr=nf_inq_varid (ncid, 'NH4',NH4G_id)
          if (ierr .ne. nf_noerr) then
            write(*,1) 'NH4'
          endif

          ierr=nf_inq_varid (ncid, 'Fe',FEG_id)
          if (ierr .ne. nf_noerr) then
            write(*,1) 'FE'
          endif

          ierr=nf_inq_varid (ncid, 'O2',O2G_id)
          if (ierr .ne. nf_noerr) then
            write(*,1) 'O2'
          endif

          ierr=nf_inq_varid (ncid, 'DIC',DICG_id)
          if (ierr .ne. nf_noerr) then
            write(*,1) 'DIC'
          endif

          ierr=nf_inq_varid (ncid, 'Alk',ALKG_id)
          if (ierr .ne. nf_noerr) then
            write(*,1) 'ALK'
          endif

          ierr=nf_inq_varid (ncid, 'DOC',DOCG_id)
          if (ierr .ne. nf_noerr) then
            write(*,1) 'DOC'
          endif

          ierr=nf_inq_varid (ncid, 'DON',DONG_id)
          if (ierr .ne. nf_noerr) then
            write(*,1) 'DON'
          endif

          ierr=nf_inq_varid (ncid, 'DOP',DOPG_id)
          if (ierr .ne. nf_noerr) then
            write(*,1) 'DOP'
          endif

          ierr=nf_inq_varid (ncid, 'NO2',NO2G_id)
          if (ierr .ne. nf_noerr) then
            write(*,1) 'NO2'
          endif

!!  Print diagnosis

        if (allocated(C)) deallocate (C)
        if (allocated(D)) deallocate (D)

      endif             ! iic.eq.ntstart

  1   format(1x,'WARNING: get_psource :: Cannot find variable ''',
     &                          A, ''' in file ''', A, '''.')



        call advance_cycle (psrc_cycle, ntpsrc,
     &                         psrc_ncycle, psrc_rec, ierr)


            start(2)=1           ! sources
            count(2)=Nsrc        ! sources
            start(1)=1           ! time 
            count(1)=ntpsrc      ! time

! time

          ierr=nf_get_vara_FTYPE (ncid, psrc_tid, start(1),count(1)
     &        , psrc_time)

! Qbar
           ierr=nf_get_vara_FTYPE (ncid, QbarG_id, start,count,
     &                                           QbarG)
!! Tsrc
             ierr=nf_get_vara_FTYPE (ncid, tempG_id, start, count,
     &                                           tempG)

             ierr=nf_get_vara_FTYPE (ncid, saltG_id, start, count,
     &                                           saltG)

             ierr=nf_get_vara_FTYPE (ncid, PO4G_id, start, count,
     &                                           PO4G)

             ierr=nf_get_vara_FTYPE (ncid, NO3G_id, start, count,
     &                                           NO3G)

             ierr=nf_get_vara_FTYPE (ncid, NH4G_id, start, count,
     &                                           NH4G)

             ierr=nf_get_vara_FTYPE (ncid, FEG_id, start, count,
     &                                           FEG)

             ierr=nf_get_vara_FTYPE (ncid, O2G_id, start, count,
     &                                           O2G)

             ierr=nf_get_vara_FTYPE (ncid, DICG_id, start, count,
     &                                           DICG)

             ierr=nf_get_vara_FTYPE (ncid, ALKG_id, start, count,
     &                                           ALKG)

             ierr=nf_get_vara_FTYPE (ncid, DOCG_id, start, count,
     &                                           DOCG)

             ierr=nf_get_vara_FTYPE (ncid, DONG_id, start, count,
     &                                           DONG)

             ierr=nf_get_vara_FTYPE (ncid, DOPG_id, start, count,
     &                                           DOPG)

             ierr=nf_get_vara_FTYPE (ncid, NO2G_id, start, count,
     &                                           NO2G)
 

! Read data from the file: Check if model time is bounded by the past
!----- ---- ---- --- ----- and the future data times: if not, increment
! record and cycling indices, flip time index and read a new portion of
! the data. Repeat until model time falls between the two data times.

  2   format(/1x,'### ERROR: get_psource :: Cannot read variable ''',
     &                                A, ''' for rec =', I4/12x,A/)
      return

  99  may_day_flag=2
!      return
      end


!!!!!! LAST PART !!!!!!!!!!!!!!!!!
!!!!!! TIME INTERPOLATION  !!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine set_psource_tile (istr,iend,jstr,jend,ierr)

! Set-up psource data for current tile.

      implicit none
# include "param.h"
# include "sources.h"
# include "forces.h"
# include "scalars.h"
# include "ocean3d.h"
      integer istr,iend,jstr,jend, ierr, it1,it2
      integer is, itrc, k
      real cff, cff1,cff2,cff3,cff4

c**
c         integer itpsrc
c         real psrc_time(2), psrc_cycle
c         real Tsrc0(Msrc,NT), TsrcG(Msrc,2,NT)

c**

# include "compute_extended_bounds.h"


!         write(*,*) 'time model = ',(time+0.5*dt)*sec2day

        do while (psrc_time(itpsrc).lt.(time+0.5*dt)*sec2day)
        itpsrc = min(itpsrc+1,ntpsrc+1)
        enddo

!         write(*,*) 'itpsrc1 = ', itpsrc

      it1=itpsrc-1             ! Always set synchro_flag, regardless
      it2=itpsrc                ! whether the boundary is needed or not:
      cff=(time+0.5*dt)*sec2day          ! this is necessary to make sure that
      cff1=psrc_time(it2)-cff   ! running out of data records condition
      cff2=cff-psrc_time(it1)   ! is the same for all MPI nodes.

      if ((cff1.ge.0. .and. cff2.ge.0.) .and. (itpsrc.lt.ntpsrc+1)) then

        cff=1./(cff1+cff2)
        cff3=cff1*cff
        cff4=cff2*cff

!        write(*,*) 'cff_1 = ',cff3
!        write(*,*) 'cff_2 = ',cff4

        do is=1,Nsrc

          do k=1,N

            Qbar(is)=(cff3*QbarG(it1,is)+cff4*QbarG(it2,is))
            Qsrc(is,k)=(cff3*QbarG(it1,is)+cff4*QbarG(it2,is))*Qshape(is,k)
            Qsrctemp(k)=cff3*tempG(it1,is)+cff4*tempG(it2,is)
            Qsrcsalt(k)=cff3*saltG(it1,is)+cff4*saltG(it2,is)
            QsrcPO4(k)=cff3*PO4G(it1,is)+cff4*PO4G(it2,is)
            QsrcNO3(k)=cff3*NO3G(it1,is)+cff4*NO3G(it2,is)
            QsrcNH4(k)=cff3*NH4G(it1,is)+cff4*NH4G(it2,is)
            QsrcFE(k)=cff3*FEG(it1,is)+cff4*FEG(it2,is)
            QsrcO2(k)=cff3*O2G(it1,is)+cff4*O2G(it2,is)
            QsrcDIC(k)=cff3*DICG(it1,is)+cff4*DICG(it2,is)
            QsrcALK(k)=cff3*ALKG(it1,is)+cff4*ALKG(it2,is)
            QsrcDOC(k)=cff3*DOCG(it1,is)+cff4*DOCG(it2,is)
            QsrcDON(k)=cff3*DONG(it1,is)+cff4*DONG(it2,is)
            QsrcDOP(k)=cff3*DOPG(it1,is)+cff4*DOPG(it2,is)
            QsrcNO2(k)=cff3*NO2G(it1,is)+cff4*NO2G(it2,is)

                 Tsrc(is,k,itemp)=Qsrctemp(k)
                 Tsrc(is,k,isalt)=Qsrcsalt(k)
                 Tsrc(is,k,iPO4)=QsrcPO4(k)
                 Tsrc(is,k,iNO3)=QsrcNO3(k)
                 Tsrc(is,k,iNH4)=QsrcNH4(k)
                 Tsrc(is,k,iFE)=QsrcFE(k)
                 Tsrc(is,k,iO2)=QsrcO2(k)
                 Tsrc(is,k,iDIC)=QsrcDIC(k)
                 Tsrc(is,k,iALK)=QsrcALK(k)
                 Tsrc(is,k,iDOC)=QsrcDOC(k)
                 Tsrc(is,k,iDON)=QsrcDON(k)
                 Tsrc(is,k,iDOP)=QsrcDOP(k)
                 Tsrc(is,k,iNO2)=QsrcNO2(k)

          enddo   ! k
        enddo ! is


!  Unable to set-up Qbar and Tsrc0:
!  Complain about the error and signal to quit.
      else
        if (ZEROTH_TILE) then
            write(*,'(/1x,3A/3(1x,A,F16.10)/)')      '### WARNING: ',
     &    'set_psource_tile :: Current model time is outside bounds ',
     &    'of ''psrc_time''.', 'psrc_tstart =', psrc_time(it1)*sec2day,
     &                         'tdays =',       tdays,
     &                         'psrc_tend =',   psrc_time(it2)*sec2day
          may_day_flag=2
        ierr=ierr+1
        endif
!        ierr=ierr+1
      endif
      return
      end
#else
      subroutine get_psource_empty
      return
      end
#endif /* PSOURCE && !ANA_PSOURCE */
