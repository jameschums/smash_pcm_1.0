!----------------------------------------------------------------------------------------
  subroutine calcrgeometry(converged,nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
!----------------------------------------------------------------------------------------
!
! Driver of geometry optimization calculation
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
      use modmolecule, only : nmo, natom, coord, coordold
      use modopt, only : nopt, optconv, cartesian
      use modwarn, only : nwarn
      use modjob, only : method
      use moddft, only : idftex, idftcor
      use modguess, only : guess
      use modprint, only : iprint
      use modscf, only : dconv
      use modthresh, only : cutint2, threshover
      use modprop, only : octupole
      implicit none
      integer,intent(in) :: nproc1, nproc2, myrank1, myrank2, mpi_comm1, mpi_comm2
      integer,allocatable :: iredun(:)
      integer :: nao2, nao3, nshell3, natom3, ii, iopt
      integer :: isizered, numbond, numangle, numtorsion, numredun, maxredun
      real(8), parameter :: third=0.3333333333333333D+00
      real(8), allocatable :: h1mtrx(:), smtrx(:), tmtrx(:), cmo(:), ortho(:), dmtrx(:)
      real(8), allocatable :: xint(:), energymo(:)
      real(8), allocatable :: egrad(:), egradold(:), ehess(:)
      real(8), allocatable :: overinv(:), work(:,:)
      real(8), allocatable :: workv(:), coordredun(:), egradredun(:)
      real(8) :: egradmax, egradrms
      real(8) :: savedconv, savecutint2
      logical,intent(out) :: converged
      logical :: exceed
!
      nao2= nao*nao
      nao3=(nao*(nao+1))/2
      nshell3=(nshell*(nshell+1))/2
      natom3= natom*3
      converged=.false.
!
! Calculate redundant coordinate
!
      if(.not.cartesian) then
        isizered= natom*4*10
        call memset(isizered)
        allocate(iredun(isizered))
        do ii= 1,10
          call setredundantcoord(iredun,isizered,numbond,numangle,numtorsion,exceed)
          if(.not.exceed) exit
          call memunset(isizered)
          deallocate(iredun)
          isizered= isizered*2
          call memset(isizered)
          allocate(iredun(isizered))
          if(ii == 10) then
            write(*,'(" Error! The array size for redundant coordinate is too large.")')
            call iabort
          endif
        enddo
        numredun= numbond+numangle+numtorsion
        maxredun= max(numredun,natom3)
      endif
!
! Set arrays for energy
!
      call memset(nao3*4+nao2*2+nshell3+nao)
      allocate(h1mtrx(nao3),smtrx(nao3),tmtrx(nao3),cmo(nao2),ortho(nao2),dmtrx(nao3), &
&              xint(nshell3),energymo(nao))
!
! Set arrays for energy gradient and geometry optimization
!
      if(cartesian) then
        call memset(natom3*2+natom3*(natom3+1)/2)
        allocate(egrad(natom3),egradold(natom3),ehess(natom3*(natom3+1)/2))
      else
        call memset(natom3+numredun*4+numredun*(numredun+1)/2)
        allocate(egrad(natom3),coordredun(numredun*2),egradredun(numredun*2), &
&                ehess(numredun*(numredun+1)/2))
      endif
!
! Start geometry optimization cycle
!
      do iopt= 1,nopt
!
! Print geometry
!
        if(iopt >= 2) call writegeom
!
! Calculate nuclear repulsion energy
!
        call nucenergy
        if(master) then
          write(*,'(" Nuclear repulsion energy =",f15.8," a.u.",/)') enuc
        endif
!
! Set work arrays 1
!
        call memset(nao2*2)
        allocate(overinv(nao2),work(nao,nao))
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
        if(iopt == 1) then
          call guessmo(cmo,cmo,overinv,h1mtrx,ortho, &
&                      nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
        endif
!
! Unset work arrays 1
!
        deallocate(overinv,work)
        call memunset(nao2*2)
        call tstamp(1)
!
! Calculate energy
!
        if((method == 'HARTREE-FOCK').or.(method == 'MP2')) then
          call calcrhf(h1mtrx,cmo,ortho,smtrx,dmtrx,xint,energymo, &
&                      nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
          if(iopt == 1) call writeeigenvalue(energymo,energymo,1)
          call tstamp(1)
        elseif((idftex >= 1).or.(idftcor >= 1)) then
          if((iopt == 1).and.(guess == 'HUCKEL')) then
            savedconv= dconv
            savecutint2= cutint2
            dconv= max(dconv,1.0D-2)
            cutint2= max(cutint2,1.0D-9)
            call calcrhf(h1mtrx,cmo,ortho,smtrx,dmtrx,xint,energymo, &
&                        nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
            dconv= savedconv
            cutint2= savecutint2
            call tstamp(1)
          endif
          call calcrdft(h1mtrx,cmo,ortho,smtrx,dmtrx,xint,energymo, &
&                       nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
          if(iopt == 1) call writeeigenvalue(energymo,energymo,1)
          call tstamp(1)
        else
          if(master) then
            write(*,'(" Error! This program does not support method= ",a16,".")')method
            call iabort
          endif
        endif
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
        call calcmaxgrad(egradmax,egradrms,egrad,natom3)
        if(master) write(*,'(" Optimization Cycle",i4,4x,"Maximum gradient =",f11.6,4x, &
&                            "RMS gradient =",f11.6,/)') iopt,egradmax,egradrms
!
! Check convergence
!
        if((egradmax <= optconv).and.(egradrms <= optconv*third)) then
          if(master) write(*,'(" Geometry converged.",/)')
          converged=.true.
          exit
        endif
!
! Write checkpoint file
!
        if(master.and.(check /= '')) call writecheck(cmo,cmo,dmtrx,dmtrx,energymo,energymo)
!
! Set work arrays 2
!
        if(cartesian) then
          call memset(natom3*3)
          allocate(workv(natom3*3))
        else
          call memset(maxredun*maxredun*4+maxredun*3)
          allocate(work(maxredun*maxredun,4),workv(maxredun*3))
        endif
!
! Calculate new coordinate
!
        if(cartesian) then
          call calcnewcoord(coord,coordold,egrad,egradold,ehess,workv,natom3,iopt, &
&                           nproc2,myrank2,mpi_comm2)
        else
          call calcnewcoordred(coord,coordold,coordredun,egrad,egradredun,ehess,work(1,1), &
&                              work(1,2),work(1,3),work(1,4),workv,iopt,iredun,isizered, &
&                              maxredun,numbond,numangle,numtorsion,numredun, &
&                              nproc2,myrank2,mpi_comm2)
        endif
!
! Unset work arrays 2
!
        if(cartesian) then
          deallocate(workv)
          call memunset(natom3*3)
        else
          deallocate(work,workv)
          call memunset(maxredun*maxredun*4+maxredun*3)
        endif
!
! Set guess MO calculation flag from Huckel to projection
!
        call setnextopt(coordold,natom,iopt)
!
        if((iopt == nopt).and.master) then
          write(*,'("Warning! Geometry did not converge.")')
          nwarn= nwarn+1
          exit
        endif
        call tstamp(1)
      enddo
!
! End of optimization cycle 
!
      call writeeigenvalue(energymo,energymo,1)
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
        allocate(work(nao3,29))
        call calcroctupole(work,work(1,4),work(1,10),work(1,20),dmtrx, &
&                          nproc1,myrank1,mpi_comm1)
        deallocate(work)
        call memunset(nao3*29)
      else
!
! Calculate dipole moment
!
        call memset(nao3*6)
        allocate(work(nao3,6))
        call calcrdipole(work,work(1,4),dmtrx,nproc1,myrank1,mpi_comm1)
        deallocate(work)
        call memunset(nao3*6)
      endif
!
! Write optimized geometry
!
      if(master.and.converged) then
        write(*,'(" ==========================")')
        write(*,'("     Optimized Geometry")')
        write(*,'(" ==========================")')
        call writegeom
      endif
!
! Write checkpoint file
!
      if(master.and.(check /= '')) call writecheck(cmo,cmo,dmtrx,dmtrx,energymo,energymo)
!
! Unset arrays for energy gradient and geometry optimization
!
      if(cartesian) then
        deallocate(egrad,egradold,ehess)
        call memunset(natom3*2+natom3*(natom3+1)/2)
      else
        deallocate(egrad,coordredun,egradredun, &
&                  ehess)
        call memunset(natom3+numredun*4+numredun*(numredun+1)/2)
      endif
!
! Unset arrays for energy
!
      deallocate(h1mtrx,smtrx,tmtrx,cmo,ortho,dmtrx, &
&                xint,energymo)
      call memunset(nao3*4+nao2*2+nshell3+nao)
!
! Unset array for redundant coordinate
!
      if(.not.cartesian) then
        deallocate(iredun)
        call memunset(isizered)
      endif
!
      call tstamp(1)
      return
end


