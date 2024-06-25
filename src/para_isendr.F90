!--------------------------------------------------------------
  subroutine para_isendr(buff,num,idest,ntag,mpi_commin,ireq)
!--------------------------------------------------------------
!
! Non-blocking MPI_Isend of real(8) data
!
#ifndef noMPI
#ifndef ILP64
      use mpi
      implicit none
      integer,intent(in) :: num, idest, ntag, mpi_commin
      integer,intent(out) :: ireq
      integer(selected_int_kind(9)) :: num4, idest4, ntag4, ireq4, mpi_comm4, ierr
#else
      use mpi
      implicit none
      integer,intent(in) :: num, idest, ntag, mpi_commin
      integer,intent(out) :: ireq
      integer(selected_int_kind(18)) :: num4, idest4, ntag4, ireq4, mpi_comm4, ierr
#endif
      real(8),intent(in) :: buff(*)
      num4= num
      idest4= idest
      ntag4= ntag
      mpi_comm4= mpi_commin
      call mpi_isend(buff,num4,MPI_DOUBLE_PRECISION,idest4,ntag4,mpi_comm4,ireq4,ierr)
      ireq= ireq4
#endif
      return
end
