!-------------------------
  subroutine setparallel
!-------------------------
!
!  Initialize MPI execution environment
!
      use modparallel, only : master, nproc1, nproc2, myrank1, myrank2, &
&                             mpi_comm1, mpi_comm2
      implicit none
!
! Initialize variables for parallelization
!
      master = .true.
!
! Start MPI parallelization and set mpi_comm1=MPI_COMM_WORLD
!
      write(*,'("james - one day I will work out what the hell this does.. mp2 smash")')
      call para_init(mpi_comm1)
      call para_comm_size(nproc1,mpi_comm1)
      call para_comm_rank(myrank1,mpi_comm1)
!
      nproc2= nproc1
      myrank2= myrank1
      mpi_comm2= mpi_comm1
!
      if(nproc1.gt.1) then
        master =(myrank1 == 0)
      endif
!
      return
end


