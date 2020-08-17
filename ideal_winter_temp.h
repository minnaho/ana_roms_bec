! include file for initial conditions in ana_initial.F
! idealized linearly stratified temp conditions 

           t(i,j,1,1,itemp)  = 11.6552899348752   ! bottom
           t(i,j,2,1,itemp)  = 11.689323697713162
           t(i,j,3,1,itemp)  = 11.725761373059267
           t(i,j,4,1,itemp)  = 11.764695819647693
           t(i,j,5,1,itemp)  = 11.806210463360992
           t(i,j,6,1,itemp)  = 11.850376858280255
           t(i,j,7,1,itemp)  = 11.89725210910698
           t(i,j,8,1,itemp)  = 11.946876207764266
           t(i,j,9,1,itemp)  = 11.999269351622347
!           t(i,j,10,1,itemp) = 12.054429325191242
!           t(i,j,11,1,itemp) = 12.154429325191242
!           t(i,j,12,1,itemp) = 12.254429325191242
!           t(i,j,13,1,itemp) = 12.354429325191242
!           t(i,j,14,1,itemp) = 12.454429325191242
!           t(i,j,15,1,itemp) = 12.654429325191242
!           t(i,j,16,1,itemp) = 12.754429325191242
!           t(i,j,17,1,itemp) = 12.854429325191242
!           t(i,j,18,1,itemp) = 12.954429325191242
!           t(i,j,19,1,itemp) = 13.054429325191242
!           t(i,j,20,1,itemp) = 13.154429325191242
           t(i,j,k,2,itemp)=t(i,j,k,1,isalt)

