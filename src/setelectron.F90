!-------------------------
  subroutine setelectron
!-------------------------
!
! Set number of electrons
!
      use modparallel, only : master
      use modmolecule, only : numatomic, neleca, nelecb, natom, multi, charge
      use modecp, only : flagecp, izcore
      use modjob, only : scftype
      use modwarn, only : nwarn
      implicit none
      integer :: nume, ii
!
! Calculate total number of electrons
!
      nume= -nint(charge)
      do ii= 1,natom
        if(numatomic(ii) > 0) nume= nume+numatomic(ii)
      enddo
!
! Subtract electrons of core potentials
!
      if(flagecp) then
        do ii= 1,natom
          nume= nume-izcore(ii)
        enddo
      endif

!
! Calculate numbers of alpha and beta electrons
!
      if((scftype == 'RHF').and.(multi /= 1)) then
        if(master) write(*,'(" Warning! SCFtype changes from RHF to UHF.")')
        scftype = 'UHF'
        nwarn= nwarn+1
      endif
!
      neleca=(nume+multi-1)/2
      nelecb=(nume-multi+1)/2
      if((neleca+nelecb)/= nume) then
        if(master) write(*,'(" Error! Spin multiplicity is ",i2, &
&                               ", but number of elctrons is ",i5,". so you are a twit..abort")')multi, nume
        call iabort
      endif
!
      return
end

