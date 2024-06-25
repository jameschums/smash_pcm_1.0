!----------------------------------------------------------------
  subroutine para_irecvr(buff,num,isource,ntag,mpi_commin,ireq)
!----------------------------------------------------------------
!
! Non-blocking MPI_Irecv of real(8) data
!
#ifndef noMPI
#ifndef ILP64
      use mpi
      implicit none
      integer,intent(in) :: num, isource, ntag, mpi_commin
      integer,intent(out) :: ireq
      integer(selected_int_kind(9)) :: num4, isource4, ntag4, ireq4, mpi_comm4, ierr
#else
      use mpi
      implicit none
      integer,intent(in) :: num, isource, ntag, mpi_commin
      integer,intent(out) :: ireq
      integer(selected_int_kind(18)) :: num4, isource4, ntag4, ireq4, mpi_comm4, ierr
#endif
      real(8),intent(out) :: buff(*)
      num4= num
      isource4= isource
      ntag4= ntag
      mpi_comm4= mpi_commin
      call mpi_irecv(buff,num4,MPI_DOUBLE_PRECISION,isource4,ntag4,mpi_comm4,ireq4,ierr)
      ireq= ireq4
#endif
      return
end

