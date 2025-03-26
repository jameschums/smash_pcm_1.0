!----------------------------------------------------------------------------
  subroutine para_allgathervr(sbuff,num,rbuff,idisa,idisb,nproc,mpi_commin)
!----------------------------------------------------------------------------
!
! Gather data from all tasks and deliver the combined data to all tasks in mpi_comm
!
#ifndef noMPI
#ifndef ILP64
      use mpi
      implicit none
      integer,intent(in) :: num, nproc, idisa(nproc), idisb(nproc), mpi_commin
      integer(selected_int_kind(9)) :: num4, idisa4(nproc), idisb4(nproc), mpi_comm4, ierr
!      !write(*,'("james para_allgathervr line 13 use mpi")')
#else
      use mpi
      implicit none
      integer,intent(in) :: num, nproc, idisa(nproc), idisb(nproc), mpi_commin
      integer(selected_int_kind(18)) :: num4, idisa4(nproc), idisb4(nproc), mpi_comm4, ierr
!      write(*,'("james para_allgathervr line 19 use mpi")')
#endif
      real(8),intent(in) :: sbuff(*)
      real(8),intent(out):: rbuff(*)
      integer :: ii
!
      num4= num
      mpi_comm4= mpi_commin
      do ii= 1,nproc
        idisa4(ii)= idisa(ii)
        idisb4(ii)= idisb(ii)
      enddo
!      write(*,'("james para_allgathervr line 31 num4 mpi")')
      call mpi_allgatherv(sbuff,num4,mpi_real8,rbuff,idisa4,idisb4,mpi_real8,mpi_comm4,ierr)
#else
      implicit none
      integer,intent(in) :: num, nproc, idisa(nproc), idisb(nproc), mpi_commin
      real(8),intent(in) :: sbuff(*)
      real(8),intent(out):: rbuff(*)
!      write(*,'("james para_allgathervr line 37 use mpi")')
!
      call dcopy(num,sbuff,1,rbuff,1)
#endif
!
      return
end
