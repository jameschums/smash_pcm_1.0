!---------------------------------------------------------
  subroutine para_allreducei(sbuff,rbuff,num,mpi_commin)
!---------------------------------------------------------
!
! Accumulate integer values from all processes and 
! distributes the result back to all processes in mpi_comm
!
#ifndef noMPI
#ifndef ILP64
      use mpi
      use modparallel, only : checkintsize
      implicit none
      integer(selected_int_kind(9)) :: num4, mpi_comm4, ierr
!      write(*,'("james line 14 use mpi para_allreducei")')
#else
      use mpi
      use modparallel, only : checkintsize
      implicit none
      integer(selected_int_kind(18)) :: num4, mpi_comm4, ierr
!      write(*,'("james line 20 use mpif.h para_allreducei")')
#endif
      integer,intent(in) :: num, mpi_commin
      integer,intent(in) :: sbuff(*)
      integer,intent(out) :: rbuff(*)
      integer :: isize
!
      num4= num
      mpi_comm4= mpi_commin
!      write(*,'("james line 29 use num4 para_allreducei")')
!
      call checkintsize(isize)
      if(isize == 4) then
        call mpi_allreduce(sbuff,rbuff,num4,mpi_integer4,MPI_SUM,mpi_comm4,ierr)
      elseif(isize == 8) then
        call mpi_allreduce(sbuff,rbuff,num4,mpi_integer8,MPI_SUM,mpi_comm4,ierr)
      endif
#else
      implicit none
      integer,intent(in) :: num, mpi_commin
      integer,intent(in) :: sbuff(*)
      integer,intent(out) :: rbuff(*)
      integer ii
!      write(*,'("james line 43 use mpi_commin para_allreducei")')
!
      do ii= 1,num
        rbuff(ii)= sbuff(ii)
      enddo
#endif
      return
end
