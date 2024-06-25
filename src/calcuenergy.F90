!----------------------------------------------------------------------------
  subroutine calcuenergy(nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
!----------------------------------------------------------------------------
!
! Driver of open-shell energy calculation
!
! Parallel information
!   nproc1, myrank1, mpi_comm1 : MPI_COMM_WORLD (all nodes)
!   nproc2, myrank2, mpi_comm2 : new communicator for matrix operations
!                               (default: MPI_COMM_WORLD)
!
      use modparallel, only : master
      use modiofile, only : check
      use modbasis, only : nao, nshell
      use modenergy, only : enuc
      use modmolecule, only : nmo
      use modjob, only : method
      use moddft, only : idftex, idftcor
      use modguess, only : guess
      use modprint, only : iprint
      use modscf, only : dconv
      use modthresh, only : cutint2, threshover
      use modprop, only : octupole
      implicit none
      integer,intent(in) :: nproc1, nproc2, myrank1, myrank2, mpi_comm1, mpi_comm2
      integer :: nao2, nao3, nshell3
      real(8), allocatable :: h1mtrx(:), smtrx(:), tmtrx(:), cmoa(:), cmob(:), ortho(:)
      real(8), allocatable :: dmtrxa(:), dmtrxb(:), xint(:), energymoa(:), energymob(:)
      real(8), allocatable :: overinv(:), work(:)
      real(8) :: savedconv, savecutint2
!
      nao2= nao*nao
      nao3=(nao*(nao+1))/2
      nshell3=(nshell*(nshell+1))/2
!
! Set arrays 1
!
      call memset(nao3*5+nao2*3+nshell3+nao*2)
      allocate(h1mtrx(nao3),smtrx(nao3),tmtrx(nao3),cmoa(nao2),cmob(nao2),ortho(nao2),&
&              dmtrxa(nao3),dmtrxb(nao3),xint(nshell3),energymoa(nao),energymob(nao))
!
! Calculate nuclear repulsion energy
!
      call nucenergy
      if(master) then
        write(*,'(" Nuclear repulsion energy =",f15.8," a.u.",/)') enuc
      endif
!
! Set arrays 2
!
      call memset(nao2*2)
      allocate(overinv(nao2),work(nao2))
!
! Calculate overlap and 1-electron integrals
!
      call oneei(h1mtrx,smtrx,tmtrx,work,nproc2,myrank2,mpi_comm2)
!
! Calculate canonicalization and inverse overlap matrices
!
      call fullmtrx(smtrx,work,nao)
      call mtrxcanoninv(ortho,overinv,work,nao,nmo,threshover,nproc2,myrank2,mpi_comm2)
!
! Calculate initial MOs
!
      call guessmo(cmoa,cmob,overinv,h1mtrx,ortho, &
&                  nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
!
! Unset arrays 2
!
      deallocate(overinv,work)
      call memunset(nao2*2)
      call tstamp(1)
!
! Start SCF
!
      if(method == 'HARTREE-FOCK') then
        call calcuhf(h1mtrx,cmoa,cmob,ortho,smtrx,dmtrxa,dmtrxb,xint,energymoa,energymob, &
&                    nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
        call writeeigenvalue(energymoa,energymob,2)
        call tstamp(1)
      elseif((idftex >= 1).or.(idftcor >= 1)) then
        if(guess == 'HUCKEL') then
          savedconv= dconv
          savecutint2= cutint2
          dconv= max(dconv,1.0D-2)
          cutint2= max(cutint2,1.0D-9)
          call calcuhf(h1mtrx,cmoa,cmob,ortho,smtrx,dmtrxa,dmtrxb,xint,energymoa,energymob, &
&                      nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
          dconv= savedconv
          cutint2= savecutint2
          call tstamp(1)
        endif
        call calcudft(h1mtrx,cmoa,cmob,ortho,smtrx,dmtrxa,dmtrxb,xint,energymoa,energymob, &
&                     nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
        call writeeigenvalue(energymoa,energymob,2)
        call tstamp(1)
!     elseif(method == 'MP2') then
!       call calcuhf(h1mtrx,cmoa,ortho,smtrx,xint,energymoa)
!       call tstamp(1)
!       call calcump2(cmoa,energymoa,xint)
!       call tstamp(1)
      else
        if(master) then
          write(*,'(" Error! This program does not support method= ",a16,".")')method
          call iabort
        endif
      endif
!
! Print MOs
!
      if(master.and.(iprint >= 2)) then
        write(*,'("  -------------------------")')
        write(*,'("    Alpha MO coefficients")')
        write(*,'("  -------------------------")')
        call writeeigenvector(cmoa,energymoa)
        write(*,'("  ------------------------")')
        write(*,'("    Beta MO coefficients")')
        write(*,'("  ------------------------")')
        call writeeigenvector(cmob,energymob)
      endif
!
! Calculate Mulliken charge
!
      call calcumulliken(dmtrxa,dmtrxb,smtrx)
!
! Calculate dipole, quadrupole, and octupole moments
!
      if(octupole) then
        call memset(nao3*29)
        allocate(work(nao3*29))
        call calcuoctupole(work,work(nao3*3+1),work(nao3*9+1),work(nao3*19+1),dmtrxa,dmtrxb, &
&                          nproc1,myrank1,mpi_comm1)
        deallocate(work)
        call memunset(nao3*29)
      else
!
! Calculate dipole moment
!
        call memset(nao3*6)
        allocate(work(nao3*6))
        call calcudipole(work,work(nao3*3+1),dmtrxa,dmtrxb,nproc1,myrank1,mpi_comm1)
        deallocate(work)
        call memunset(nao3*6)
      endif
!
! Write checkpoint file
!
      if(master.and.(check /= '')) call writecheck(cmoa,cmob,dmtrxa,dmtrxb,energymoa,energymob)
!
! Unset arrays 1
!
      deallocate(h1mtrx,smtrx,tmtrx,cmoa,cmob,ortho, &
&                dmtrxa,dmtrxb,xint,energymoa,energymob)
      call memunset(nao3*5+nao2*3+nshell3+nao*2)
      return
end


