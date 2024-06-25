!-------------------------------------
  subroutine para_waitall(nump,ireq)
!-------------------------------------
!
!  MPI_Waitall
!
#ifndef noMPI
#ifndef ILP64
      use mpi
      implicit none
      integer,intent(in) :: nump, ireq(nump)
      integer(selected_int_kind(9)) :: nump4, ireq4(nump), STATUS(MPI_STATUS_SIZE,nump), ierr
#else
      use mpi
      implicit none
      integer,intent(in) :: nump, ireq(nump)
      integer(selected_int_kind(18)) :: nump4, ireq4(nump), STATUS(MPI_STATUS_SIZE,nump), ierr
#endif
      integer :: ii
      nump4= nump
      do ii= 1,nump
        ireq4(ii)= ireq(ii)
      enddo
      call mpi_waitall(nump4,ireq4,STATUS,ierr)
#endif
      return
end
