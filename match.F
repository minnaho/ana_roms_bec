      program match

! Finds and reports unmatched curly { }-brackets in LaTeX file.

      implicit none
      integer, parameter :: max_wdth=128
      character(len=1) S(max_wdth)
      character(len=64) fname
      character(len=7) fmt
      character(len=4) mess
      logical math, math1, math2, mathbracket,  end_of_file
      integer nargs,iargc, wdth, i,iocheck, line, TSUM,TDIF,LSUM,LDIF
      character(len=1) backslash
      backslash=char(92)


      nargs=iargc()
      if (nargs.eq.1) then
        call getarg (nargs,fname)
      else
        write(*,'(/1x,2A/1x,A/24x,A/)') 'This is "match" - a tool ',
     &       'to check parity of curly { }-brackets in LaTeX file.',
     &               'It should be used as:',  'match filename.tex'
        stop
      endif



      TSUM=0              ; math=.false.
      TDIF=0              ; math1=.false.
      end_of_file=.false. ; math2=.false.

      write(*,*) 'Processing ',fname,'...'
      open(20,file=fname,status='old',form='formatted')

      line=0
      do while (.not.end_of_file)
        line=line+1            !<-- count lines in the file
        do i=1,max_wdth
           S(i)=' '   ! <-- reset string
        enddo
        write(fmt,'(A,I3,A)') '(', max_wdth, 'A1)'
        read(20,fmt,iostat=iocheck,end=2) (S(i),i=1,max_wdth)
        go to 3
   2    end_of_file=.true.

   3    i=max_wdth                            ! Find the actual width
        do while(i.gt.1 .and. S(i).eq.' ')    ! of the line as the
          i=i-1                               ! position of the last
        enddo                                 ! non-blank character
        if (S(i).eq.' ') then
          wdth=i-1
        else
          wdth=i
        endif                                 ! Check that some
        i=0                                   ! portion of the line
        do while (i.lt.wdth)                  ! may be commented out
          i=i+1                               ! by % and restrict
          if (S(i).eq.'%') then               ! width accordingly
            if (i.gt.1) then
              if (S(i-1).ne.backslash) wdth=i-1
            else
              wdth=i-1
            endif
          endif
        enddo
        if (wdth.gt.0) then
          LSUM=0
          LDIF=0
          do i=1,wdth

              if (S(i).eq.'$') then
                math1=.not.math1
                if (.not.math1 .and. i.ne.1) then
                  if (S(i-1).eq.'$') then
                    if (.not.math2) then
                       math=.true.
                       math2=.true.
                       write(*,'(4x,A)') 'MATH MODE...'
                    else
                       math2=.false.
                    endif
                  endif
                endif
              endif

            if (S(i).eq.'{') then           ! Sort out curly brackets
              mathbracket=.false.
              if (i.gt.1) then
                if (S(i-1).eq.backslash) then
                  write(*,*)'       \\{'
                  mathbracket=.true.
                endif
              endif
              if (.not.mathbracket) then
                LSUM=LSUM+1
                  LDIF=LDIF+1
                  TSUM=TSUM+1
                  TDIF=TDIF+1
                endif
            elseif (S(i).eq.'}') then
              mathbracket=.false.
              if (i.gt.1) then
                if (S(i-1).eq.backslash) then
                  write(*,*)'       \\}'
                  mathbracket=.true.
                endif
              endif
              if (.not.mathbracket) then
                LSUM=LSUM+1
                LDIF=LDIF-1
                TSUM=TSUM+1
                TDIF=TDIF-1
              endif
            endif

          enddo

          if (wdth.lt.max_wdth .and. S(wdth+1).eq.'%') then
            MESS='comm'
          else
            MESS='    '
          endif

          if (LSUM.ne.0) then
            if(math) then
               write(*,98) line,TSUM,TDIF,MESS
  98           format(10x,'line=',I4,' {+}=',I5,' {-}=',I2,2x,A4)
            else
               write(*,99) line,TSUM,TDIF,MESS
  99           format('line=',I4,10x,' {+}=',I5,' {-}=',I2,2x,A4)
            endif
          endif

          if (math .and. .not.math2) then
             math=.false.
             write(*,'(4x,A)') ' ...END OF MATH MODE.'
          endif

        endif  !<-- wdth.gt.0
      enddo

      stop
      end
