#include "cppdefs.h"
#ifdef MPI
 
C$    subroutine master_num_threads (nthrds)
C$    implicit none
C$    integer nthrds, trd, omp_get_num_threads, omp_get_thread_num
C$    trd=omp_get_thread_num()
C$    if (trd.eq.0) nthrds=omp_get_num_threads()
C$    return
C$    end
 
 
 
      subroutine mpi_setup (ierr)
      implicit none
      integer ierr
# include "param.h"
# include "hidden_mpi_vars.h"
# include "ncvars.h"
# ifdef HOFFMAN2
      include "mpif.h"
# else
# include "mpif.h"
# endif
C$    integer nthrds
      integer nsize, i_W, i_E, j_S, j_N, off_XI, off_ETA
# ifndef MPI_SILENT_MODE
      character(len=25) greet ! <-- greeting message to print 
# endif
 
      ocean_grid_comm=MPI_COMM_WORLD
      call MPI_Comm_size (ocean_grid_comm, nsize,  ierr)
      call MPI_Comm_rank (ocean_grid_comm, mynode, ierr)

# ifndef MPI_SILENT_MODE
      write(greet,'(A,I5,1x,A,I5)') 'MPI rank', mynode, 'out of', nsize
      call system ('echo "'/ /greet/ /' running on `hostname -s`"')
# endif


      exc_call_count=0
                                         ! Determine number of threads
C$    if (mynode.eq.0) then              ! on MPI node rank=0, then
C$OMP PARALLEL SHARED(nthrds)            ! broadcast it, so all others
C$      call master_num_threads (nthrds) ! can set the same number of
C$OMP END PARALLEL                       ! threads as the rank=0 node.
C$    endif
C$    call MPI_Bcast(nthrds, 1, MPI_INTEGER, 0, ocean_grid_comm, ierr)
C$    if (mynode.gt.0) then
C$      call omp_set_num_threads (nthrds)
C$    endif
                                   ! Check whether the number of
      if (nsize.eq.NNODES) then    ! nodes specified by -np argument
        inode=mod(mynode,NP_XI)    ! of mpiexec matches the parameter
        jnode=mynode/NP_XI         ! settings in the code, and find
                                   ! indices inode,jnode identifying
        if (NP_XI.gt.1) then       ! location of the current subdomain
# ifdef EW_PERIODIC
          west_msg_exch=.true.      ! on the "processor grid".
          east_msg_exch=.true.      ! Determine whether the subdomain
# else
          west_msg_exch=inode.gt.0        ! has neighbours on the four
          east_msg_exch=inode.lt.NP_XI-1  ! sides around it and set the
# endif
        else
          west_msg_exch=.false.     ! corresponding logical flags: e.g.
          east_msg_exch=.false.     ! west_msg_exch.eqv=.true. means tha
        endif                      ! there is a neighbor on the west
# ifndef EW_PERIODIC
        west_exchng=west_msg_exch
        east_exchng=east_msg_exch
# endif
 
                                   ! side, so this subdomain must send
        if (NP_ETA.gt.1) then      ! and expects incoming  messages
# ifdef NS_PERIODIC
          south_msg_exch=.true.     ! to/from the neighbour. Note that
          north_msg_exch=.true.     ! periodic boundaries (if any) are
# else
          south_msg_exch=jnode.gt.0         ! treated exclussively via
          north_msg_exch=jnode.lt.NP_ETA-1  ! exchange of computational
# endif
        else                        ! margins, so that in this case
          south_msg_exch=.false.     ! communications take place even
          north_msg_exch=.false.     ! if the subdomain is located near
        endif                       ! the edge of the grid.
# ifndef NS_PERIODIC
        south_exchng=south_msg_exch
        north_exchng=north_msg_exch
# endif
 
 
        i_W=mod(inode-1+NP_XI,NP_XI)
        i_E=mod(inode+1       ,NP_XI)
        j_S=mod(jnode-1+NP_ETA,NP_ETA)
        j_N=mod(jnode+1       ,NP_ETA)
 
        p_W=i_W +NP_XI*jnode     ! Determine MPI-ranks of neighbors
        p_E=i_E +NP_XI*jnode     ! from the sides and corners, which
        p_S=inode+NP_XI*j_S      ! will later be used to designate
        p_N=inode+NP_XI*j_N      ! sources of incoming and targets for
                                 ! outgoing messages. Here they are
        p_NW=i_W+NP_XI*j_N       ! set as for double-periodic grid
        p_SW=i_W+NP_XI*j_S       ! regardless of the actual boundary
        p_NE=i_E+NP_XI*j_N       ! conditions. There is no ambiguity,
        p_SE=i_E+NP_XI*j_S       ! since west_msg_exch...etc logic
                                 ! blocks the annecessary messages.
 
 
        off_XI=NP_XI*Lm-LLm
        iSW_corn=inode*Lm-off_XI/2
        if (inode.eq.0) then
          iwest=1+off_XI/2
        else
          iwest=1
        endif
        if (inode.lt.NP_XI-1) then
          ieast=Lm
        else
          ieast=Lm -(off_XI+1)/2
        endif
 
        off_ETA=NP_ETA*Mm-MMm
        jSW_corn=jnode*Mm-off_ETA/2
        if (jnode.eq.0) then
          jsouth=1+off_ETA/2
        else
          jsouth=1
        endif
        if (jnode.lt.NP_ETA-1) then
          jnorth=Mm
        else
          jnorth=Mm -(off_ETA+1)/2
        endif
 
#ifdef VERBOSE
      write(*,'(A,7I5,1x,A,I4)') 'XI:', LLm, off_XI, iSW_corn, Lm,
     & ieast-iwest+1, iwest+iSW_corn,ieast+iSW_corn, 'node=', mynode
      write(*,'(A,7I5,1x,A,I4)') 'ETA:',MMm, off_ETA, jSW_corn, Mm,
     & jnorth-jsouth+1,jsouth+jSW_corn,jnorth+jSW_corn,'node=',mynode
#endif
 
 
 
# ifdef PARALLEL_FILES
         xi_rho=ieast-iwest+1
         if (EASTERN_MPI_EDGE) then
           xi_rho=xi_rho+1
         endif
         if (WESTERN_MPI_EDGE) then
           xi_rho=xi_rho+1
           xi_u=xi_rho-1
         else
           xi_u=xi_rho
         endif
 
         eta_rho=jnorth-jsouth+1
         if (NORTHERN_MPI_EDGE) then
           eta_rho=eta_rho+1
         endif
         if (SOUTHERN_MPI_EDGE) then
           eta_rho=eta_rho+1
           eta_v=eta_rho-1
         else
           eta_v=eta_rho
         endif
# endif
 
 
c#ifdef PARALLEL_FILES
c# ifndef EW_PERIODIC
c        xi_rho=Lm
c        xi_u=xi_rho
c        if (inode.eq.0) xi_rho=xi_rho+1
c        if (inode.eq.NP_XI-1) then
c          xi_rho=xi_rho+1
c          xi_u=xi_u+1
c        endif
c# endif
c# ifndef NS_PERIODIC
c        eta_rho=Mm
c        eta_v=eta_rho
c        if (jnode.eq.0) eta_rho=eta_rho+1
c        if (jnode.eq.NP_ETA-1) then
c          eta_rho=eta_rho+1
c          eta_v=eta_v+1
c        endif
c# endif
c#endif
 
        ierr=0
      else
        mpi_master_only write(*,'(/1x,A,I4,1x,A,I3,A/)')
     &   '### ERROR: mpi_setup: number of MPI-nodes should be',
     &                         NNODES, 'instead of', nsize, '.'
        ierr=99
      endif
      return
      end
#else
      subroutine MPI_Setup_empty
      end
#endif    /* MPI */
 
 
