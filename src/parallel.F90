! Copyright 2014-2017  Kazuya Ishimura
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.
!



!-----------------------------------------------
  subroutine para_comm_size(nproc,mpi_commin)
!-----------------------------------------------
!
! Return the number of processes in mpi_comm
!
#ifndef noMPI
#ifndef ILP64
      use mpi
      implicit none
      integer(selected_int_kind(9)) :: mpi_comm4, nproc4, ierr
#else
      use mpi
      implicit none
      integer(selected_int_kind(18)) :: mpi_comm4, nproc4, ierr
#endif
      integer,intent(in) :: mpi_commin
      integer,intent(out) :: nproc
!
      mpi_comm4= mpi_commin
      call mpi_comm_size(mpi_comm4,nproc4,ierr)
      nproc= nproc4
#else
      integer,intent(in) :: mpi_commin
      integer,intent(out) :: nproc
!
      nproc= 1
#endif
      return
end


!-----------------------------------------------
  subroutine para_comm_rank(myrank,mpi_commin)
!-----------------------------------------------
!
! Return the MPI rank in mpi_comm
!
#ifndef noMPI
#ifndef ILP64
      use mpi
      implicit none
      integer(selected_int_kind(9)) :: mpi_comm4, myrank4, ierr
#else
      use mpi
      implicit none
      integer(selected_int_kind(18)) :: mpi_comm4, myrank4, ierr
#endif
      integer,intent(in) :: mpi_commin
      integer,intent(out) :: myrank
!
      mpi_comm4= mpi_commin
      call mpi_comm_rank(mpi_comm4,myrank4,ierr)
      myrank= myrank4
#else
      integer,intent(in) :: mpi_commin
      integer,intent(out) :: myrank
!
      myrank= 0
#endif
      return
end


!---------------------------
  subroutine para_finalize
!---------------------------
!
! Finalize MPI
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
!
      call mpi_finalize(ierr)
#endif
      return
end


!------------------------
  subroutine para_abort
!------------------------
#ifndef noMPI
#ifndef ILP64
      use mpi
      implicit none
      integer(selected_int_kind(9)) :: icode, ierr
#else
      use mpi 
      implicit none
      integer(selected_int_kind(18)) :: icode, ierr
#endif
!
      icode=9
      call mpi_abort(MPI_COMM_WORLD,icode,ierr)
#endif
      return
end


!----------------------------------
  subroutine checkintsize4(isize)
!----------------------------------
      implicit none
      integer(selected_int_kind(9)),intent(out) :: isize
!
      isize= 4
end


!----------------------------------
  subroutine checkintsize8(isize)
!----------------------------------
      implicit none
      integer(selected_int_kind(18)),intent(out) :: isize
!
      isize= 8
end


!----------------------------------------------------
  subroutine para_bcastr(buff,num,irank,mpi_commin)
!----------------------------------------------------
!
! Broadcast real(8) data
!
#ifndef noMPI
#ifndef ILP64
      use mpi
      implicit none
      integer(selected_int_kind(9)) :: num4, irank4, mpi_comm4, ierr
#else
      use mpi
      implicit none
      integer(selected_int_kind(18)) :: num4, irank4, mpi_comm4, ierr
#endif
      integer,intent(in) :: num, irank, mpi_commin
      real(8),intent(inout) :: buff(*)
!
      num4= num
      irank4= irank
      mpi_comm4= mpi_commin
!
      call mpi_bcast(buff,num4,mpi_real8,irank4,mpi_comm4,ierr)
#endif
      return
end


!-----------------------------------------------------
  subroutine para_bcasti(ibuff,num,irank,mpi_commin)
!-----------------------------------------------------
!
! Broadcast integer data
!
#ifndef noMPI
#ifndef ILP64
      use mpi
      use modparallel, only : checkintsize
      implicit none
      integer(selected_int_kind(9)) :: num4, irank4, mpi_comm4, ierr
#else
      use mpi
      use modparallel, only : checkintsize
      implicit none
      integer(selected_int_kind(18)) :: num4, irank4, mpi_comm4, ierr
#endif
      integer,intent(in) :: num, irank, mpi_commin
      integer,intent(inout) :: ibuff(*)
      integer :: isize
!
      num4= num
      irank4= irank
      mpi_comm4= mpi_commin
!
      call checkintsize(isize)
      if(isize == 4) then
        call mpi_bcast(ibuff,num4,mpi_integer4,irank4,mpi_comm4,ierr)
      elseif(isize == 8) then
        call mpi_bcast(ibuff,num4,mpi_integer8,irank4,mpi_comm4,ierr)
      endif
#endif
      return
end


!----------------------------------------------------
  subroutine para_bcastc(buff,num,irank,mpi_commin)
!----------------------------------------------------
!
! Broadcast character data
!
#ifndef noMPI
#ifndef ILP64
      use mpi
      implicit none
      integer(selected_int_kind(9)) :: num4, irank4, mpi_comm4, ierr
#else
      use mpi
      implicit none
      integer(selected_int_kind(18)) :: num4, irank4, mpi_comm4, ierr
#endif
      integer,intent(in) :: num, irank, mpi_commin
      character(*),intent(inout) :: buff(*)
!
      num4= num
      irank4= irank
      mpi_comm4= mpi_commin
!
      call mpi_bcast(buff,num4,mpi_character,irank4,mpi_comm4,ierr)
#endif
      return
end


!----------------------------------------------------
  subroutine para_bcastl(buff,num,irank,mpi_commin)
!----------------------------------------------------
!
! Broadcast logical data
!
#ifndef noMPI
#ifndef ILP64
      use mpi
      implicit none
      integer,intent(in) :: num, irank, mpi_commin
      integer(selected_int_kind(9)) :: num4, irank4, mpi_comm4, ierr, itmp(num)
#else
      use mpi
      implicit none
      integer,intent(in) :: num, irank, mpi_commin
      integer(selected_int_kind(18)) :: num4, irank4, mpi_comm4, ierr
      integer(selected_int_kind(9)) :: itmp(num)
#endif
      integer :: ii, myrank
      logical,intent(inout) :: buff(*)
!
      call para_comm_rank(myrank,mpi_commin)
!
      if(irank == myrank) then
        do ii= 1,num
          if(buff(ii)) then
            itmp(ii)= 1
          else
            itmp(ii)= 0
          endif
        enddo
      endif
!
      num4= num
      irank4= irank
      mpi_comm4= mpi_commin
!
      call mpi_bcast(itmp,num4,mpi_integer4,irank4,mpi_comm4,ierr)
!
      do ii= 1,num
        if(itmp(ii) == 1) then
          buff(ii)= .true.
        else
          buff(ii)= .false.
        endif
      enddo
#endif
      return
end

