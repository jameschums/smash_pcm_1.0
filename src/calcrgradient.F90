!------------------------------------------------------------------------------
  subroutine calcrgradient(nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
!------------------------------------------------------------------------------
!
! Driver of energy gradient calculation
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
      use modmolecule, only : nmo, natom
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
      real(8), allocatable :: h1mtrx(:), smtrx(:), tmtrx(:), cmo(:), ortho(:), dmtrx(:)
      real(8), allocatable :: xint(:), energymo(:)
      real(8), allocatable :: overinv(:), work(:)
      real(8), allocatable :: egrad(:)
      real(8) :: egradmax, egradrms
      real(8) :: savedconv, savecutint2
!
      nao2= nao*nao
      nao3=(nao*(nao+1))/2
      nshell3=(nshell*(nshell+1))/2
!
! Set arrays 1
!
      call memset(nao3*4+nao2*2+nshell3+nao)
      allocate(h1mtrx(nao3),smtrx(nao3),tmtrx(nao3),cmo(nao2),ortho(nao2),dmtrx(nao3),&
&              xint(nshell3),energymo(nao))
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
      call guessmo(cmo,cmo,overinv,h1mtrx,ortho, &
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
      if((method == 'HARTREE-FOCK').or.(method == 'MP2')) then
        call calcrhf(h1mtrx,cmo,ortho,smtrx,dmtrx,xint,energymo, &
&                    nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
        call writeeigenvalue(energymo,energymo,1)
        call tstamp(1)
      elseif((idftex >= 1).or.(idftcor >= 1)) then
        if(guess == 'HUCKEL') then
          savedconv= dconv
          savecutint2= cutint2
          dconv= max(dconv,1.0D-2)
          cutint2= max(cutint2,1.0D-9)
          call calcrhf(h1mtrx,cmo,ortho,smtrx,dmtrx,xint,energymo, &
&                      nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
          dconv= savedconv
          cutint2= savecutint2
          call tstamp(1)
        endif
        call calcrdft(h1mtrx,cmo,ortho,smtrx,dmtrx,xint,energymo, &
&                     nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
        call writeeigenvalue(energymo,energymo,1)
        call tstamp(1)
      else
        if(master) then
          write(*,'(" Error! This program does not support method= ",a16,".")')method
          call iabort
        endif
      endif
!
! Set arrays 3
!
      call memset(natom*3)
      allocate(egrad(natom*3))
!
! Calculate energy gradient
!
      if(method == 'HARTREE-FOCK') then
        call calcgradrhf(cmo,energymo,xint,egrad,nproc1,myrank1,mpi_comm1)
        call tstamp(1)
      elseif((idftex >= 1).or.(idftcor >= 1)) then
        call calcgradrdft(cmo,energymo,xint,egrad,nproc1,myrank1,mpi_comm1)
        call tstamp(1)
      elseif(method == 'MP2') then
        call calcgradrmp2(cmo,energymo,xint,egrad,nproc1,myrank1,mpi_comm1)
        call tstamp(1)
      else
        if(master) then
          write(*,'(" Error! This program does not support method= ",a16,".")')method
          call iabort
        endif
      endif
!
! Calculate maximum and root mean square gradient values
!
      call calcmaxgrad(egradmax,egradrms,egrad,natom*3)
      if(master) write(*,'(" Maximum gradient =",f13.8,"  RMS gradient =",f13.8,/)') &
&                      egradmax,egradrms
!
! Unset arrays 3
!
      deallocate(egrad)
      call memunset(natom*3)
!
! Print MOs
!
      if(master.and.(iprint >= 2)) then
        write(*,'("  -------------------")')
        write(*,'("    MO coefficients")')
        write(*,'("  -------------------")')
        call writeeigenvector(cmo,energymo)
      endif
!
! Calculate Mulliken charge
!
      call calcrmulliken(dmtrx,smtrx)
!
! Calculate dipole, quadrupole, and octupole moments
!
      if(octupole) then
        call memset(nao3*29)
        allocate(work(nao3*29))
        call calcroctupole(work,work(nao3*3+1),work(nao3*9+1),work(nao3*19+1),dmtrx, &
&                          nproc1,myrank1,mpi_comm1)
        deallocate(work)
        call memunset(nao3*29)
      else
!
! Calculate dipole moment
!
        call memset(nao3*6)
        allocate(work(nao3*6))
        call calcrdipole(work,work(nao3*3+1),dmtrx,nproc1,myrank1,mpi_comm1)
        deallocate(work)
        call memunset(nao3*6)
      endif
!
! Write checkpoint file
!
      if(master.and.(check /= '')) call writecheck(cmo,cmo,dmtrx,dmtrx,energymo,energymo)
!
! Unset arrays 1
!
      deallocate(h1mtrx,smtrx,tmtrx,cmo,ortho,dmtrx, &
&                xint,energymo)
      call memunset(nao3*4+nao2*2+nshell3+nao)
      call tstamp(1)
      return
end


