#if defined PSOURCE || defined Q_PSOURCE

! Nsrc       Number of point Sources/Sinks.
! Dsrc       Direction of point Sources/Sinks:  0 = along XI-;
!                                          1= along ETA-direction.
!                                          2 = along s-direction
! Isrc,Jsrc  i,j-grid location of point Sources/Sinks,
!              0 <= Isrc <= Lm+1;   0 =< Jsrc <= Mm+1.
! Lsrc       Logical switch indicating which tracer field to apply
!                                         the point Sources/Sinks.
! Qsrc       Mass transport profile (m3/s) of point Sources/Sinks.
! Qbar       Vertically integrated Qsrc (m3/s) of point
! QbarG      Latest two-time snapshots of vertically integrated
!              mass transport (m3/s) of point Sources/Sinks.
! Tsrc       Tracer (tracer units) point Sources/Sinks.
! TsrcG      Latest two-time snapshots of tracer (tracer units)
!              point Sources/Sinks.
! Qshape     Nondimensional shape function to distribute mass
!             ass point Sources/Sinks vertically.

      real Qbar(Msrc), Qsrc(Msrc,N),  Qshape(Msrc,N),
     &                                Tsrc(Msrc,N,NT),Tsrc0(Msrc,NT),
     &     Qtemp(Msrc),Qsalt(Msrc),QPO4(Msrc),QNO3(Msrc),
     &     QNH4(Msrc),QFE(Msrc),QO2(Msrc),QDIC(Msrc),QALK(Msrc),
     &     QDOC(Msrc),QDON(Msrc),QDOP(Msrc),QNO2(Msrc),
     &     Qsrctemp(N),Qsrcsalt(N),QsrcPO4(N),QsrcNO3(N),
     &     QsrcNH4(N),QsrcFE(N),QsrcO2(N),QsrcDIC(N),QsrcALK(N),
     &     QsrcDOC(N),QsrcDON(N),QsrcDOP(N),QsrcNO2(N)


# if defined PSOURCE
     &              , QbarG(Stepsrc,Msrc), TsrcG(Stepsrc,Msrc,NT), PO4G(Stepsrc,Msrc),
     &              tempG(Stepsrc,Msrc),saltG(Stepsrc,Msrc),
     &              NO3G(Stepsrc,Msrc),NH4G(Stepsrc,Msrc),FEG(Stepsrc,Msrc),O2G(Stepsrc,Msrc),
     &              DICG(Stepsrc,Msrc),ALKG(Stepsrc,Msrc),DOCG(Stepsrc,Msrc),DONG(Stepsrc,Msrc),
     &              DOPG(Stepsrc,Msrc),NO2G(Stepsrc,Msrc),
     &              psrc_Hz(N), psrc_pm(1),psrc_pn(1)

      real psrc_time(Stepsrc), psrc_cycle, var1d(Msrc)

      integer itpsrc, psrc_ncycle, psrc_rec, psrc_tid,
     &        QbarG_id, TsrcG_id(NT), ncidfpsrc, ntpsrc, 
     &        tempG_id,saltG_id,
     &        PO4G_id,NO3G_id,NH4G_id,FEG_id,O2G_id,
     &        DICG_id,ALKG_id,DOCG_id,DONG_id,DOPG_id,NO2G_id

      common /source_psrcdat/ psrc_time, psrc_cycle,
     &        itpsrc, psrc_ncycle, psrc_rec, psrc_tid,
     &        QbarG_id, TsrcG_id, ncidfpsrc, ntpsrc,
     &        tempG_id,saltG_id,PO4G_id,NO3G_id,NH4G_id,
     &        FEG_id,O2G_id,DICG_id,ALKG_id,DOCG_id,DONG_id,
     &        DOPG_id,NO2G_id


# endif
      integer Nsrc, Dsrc(Msrc), Isrc(Msrc), Jsrc(Msrc), indii(Msrc),indjj(Msrc),
     &        IIsrc(Msrc), JJsrc(Msrc),Mynnd(Msrc)
      logical Lsrc(Msrc,NT)
      common /sources_all/ Qbar, Qsrc, Qshape, Tsrc,
     &     Qtemp,Qsalt,QPO4,QNO3,
     &     QNH4,QFE,QO2,QDIC,QALK,
     &     QDOC,QDON,QDOP,QNO2,
     &     Qsrctemp,Qsrcsalt,QsrcPO4,QsrcNO3,
     &     QsrcNH4,QsrcFE,QsrcO2,QsrcDIC,QsrcALK,
     &     QsrcDOC,QsrcDON,QsrcDOP,QsrcNO2
# if defined PSOURCE
     &              , QbarG, TsrcG, PO4G,
     &              tempG,saltG,
     &              NO3G,NH4G,FEG,O2G,
     &              DICG,ALKG,DOCG,DONG,
     &              DOPG,NO2G,
# endif
     &                 Nsrc, Dsrc, Isrc, Jsrc, Lsrc
#endif

# if defined ANA_PSOURCE

!      real Qbar(Msrc), Qsrc(Msrc,N),  Qshape(Msrc,N),
!     &                                Tsrc(Msrc,N,NT),Tsrc0(Msrc,NT),
!     &     Qsrctemp(N),Qsrcsalt(N),QsrcPO4(N),QsrcNO3(N),
!     &     QsrcNH4(N),QsrcFE(N),QsrcO2(N),QsrcDIC(N),QsrcALK(N),
!     &     QsrcDOC(N),QsrcDON(N),QsrcDOP(N),QsrcNO2(N),
       real temp,salt,PO4,NO3,NH4,iron_FE,O2,DIC,ALK,DOC,DON,DOP,NO2

!      integer Nsrc, Dsrc(Msrc), Isrc(Msrc), Jsrc(Msrc), indii(Msrc),indjj(Msrc),
!     &        IIsrc(Msrc), JJsrc(Msrc),Mynnd(Msrc)
!      logical Lsrc(Msrc,NT)
!      common /sources_all/ Qbar, Qsrc, Qshape, Tsrc,
!     &     Qsrctemp,Qsrcsalt,QsrcPO4,QsrcNO3,
!     &     QsrcNH4,QsrcFE,QsrcO2,QsrcDIC,QsrcALK,
!     &     QsrcDOC,QsrcDON,QsrcDOP,QsrcNO2
# endif
