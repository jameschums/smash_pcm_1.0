!----------------------------------
  subroutine para_init(mpi_comm1)
!----------------------------------
!
! Start MPI and set mpi_comm1=MPI_COMM_WORLD
!
#ifndef noMPI
#ifndef ILP64
      use mpi
      implicit none
      integer(selected_int_kind(9)) :: ierr
#else
      use mpi
      implicit none
      integer(selected_int_kind(18)) :: ierr 
#endif
      integer,intent(out) :: mpi_comm1
      call mpi_init(ierr)
      mpi_comm1= MPI_COMM_WORLD
#endif
      return
end
