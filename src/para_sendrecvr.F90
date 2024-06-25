!----------------------------------------------------------------------------------------
  subroutine para_sendrecvr(sbuff,nums,idest,ntags,rbuff,numr,isource,ntagr,mpi_commin)
!----------------------------------------------------------------------------------------
!
! Send and receive real(8) data
!
#ifndef noMPI
#ifndef ILP64
      use mpi
      implicit none
      integer,intent(in) :: nums, idest, ntags, numr, isource, ntagr, mpi_commin
      integer(selected_int_kind(9)) :: nums4, idest4, ntags4, numr4, isource4, ntagr4
      integer(selected_int_kind(9)) :: mpi_comm4, ierr, STATUS(MPI_STATUS_SIZE)
#else
      use mpi
      implicit none
      integer,intent(in) :: nums, idest, ntags, numr, isource, ntagr, mpi_commin
      integer(selected_int_kind(18)) :: nums4, idest4, ntags4, numr4, isource4, ntagr4
      integer(selected_int_kind(18)) :: mpi_comm4, ierr, STATUS(MPI_STATUS_SIZE)
#endif
      real(8),intent(in) :: sbuff(*)
      real(8),intent(out) :: rbuff(*)
      nums4= nums
      idest4= idest
      ntags4= ntags
      numr4= numr
      isource4= isource
      ntagr4 = ntagr
      mpi_comm4= mpi_commin
      call mpi_sendrecv(sbuff,nums4,MPI_DOUBLE_PRECISION,idest4,ntags4, &
&                       rbuff,numr4,MPI_DOUBLE_PRECISION,isource4,ntagr4,mpi_comm4,STATUS,ierr)
#else
      implicit none
      integer,intent(in) :: nums, idest, ntags, numr, isource, ntagr, mpi_commin
      real(8),intent(in) :: sbuff(*)
      real(8),intent(out) :: rbuff(*)
      call dcopy(nums,sbuff,1,rbuff,1)
#endif
      return
end


