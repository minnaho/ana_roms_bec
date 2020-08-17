       do itrc=1,NT
          Lsrc(:,itrc)=.FALSE.
          if (itrc.eq.itemp) then
            Lsrc(:,itrc)=.true.
            ierr=nf_inq_varid (ncid,'temp', TsrcG_id(itrc))
            if (ierr .ne. nf_noerr) then
              Lsrc(:,itrc)=.false.
              write(*,2) 'temp', frcname(1:lstr)
            endif
#ifdef SALINITY
          elseif (itrc.eq.isalt) then
            Lsrc(:,itrc)=.true.
            ierr=nf_inq_varid (ncid,'salt', TsrcG_id(itrc))
            if (ierr .ne. nf_noerr) then
              Lsrc(:,itrc)=.false.
              write(*,2) 'salt', frcname(1:lstr)
            endif
#endif /* SALINITY  */
#ifdef BIOLOGY
#ifdef BEC_LIANG
          elseif (itrc.eq.iNO3_) then
            Lsrc(itrc)=.true.
            ierr=nf_inq_varid (ncid,'NO3', TsrcG_id(itrc))
            if (ierr .ne. nf_noerr) then
              Lsrc(:,itrc)=.false.
              write(*,2) 'no3', frcname(1:lstr)
c     goto 99           !--> ERROR
            endif
# endif /* BEC_LIANG */
# endif /* BIOLOGY */
