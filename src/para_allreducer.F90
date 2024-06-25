!---------------------------------------------------------
  subroutine para_allreducer(sbuff,rbuff,num,mpi_commin)
!---------------------------------------------------------
!
! Accumulate real(8) values from all processes and 
! distributes the result back to all processes in mpi_comm
!
#ifndef noMPI
#ifndef ILP64
      use mpi
      implicit none
      integer(selected_int_kind(9)) :: num4, mpi_comm4, ierr
#else
      use mpi
      implicit none
      integer(selected_int_kind(18)) :: num4, mpi_comm4, ierr
#endif
      integer,intent(in) :: num, mpi_commin
      real(8),intent(in) :: sbuff(*)
      real(8),intent(out) :: rbuff(*)
!
      num4= num
      mpi_comm4= mpi_commin
      call mpi_allreduce(sbuff,rbuff,num4,mpi_real8,MPI_SUM,mpi_comm4,ierr)
#else
      implicit none
      integer,intent(in) :: num, mpi_commin
      real(8),intent(in) :: sbuff(*)
      real(8),intent(out) :: rbuff(*)
!
      call dcopy(num,sbuff,1,rbuff,1)
#endif
      return
end
