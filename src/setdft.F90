!--------------------
  subroutine setdft
!--------------------
!
! Set functional information
! Adjust the numbe of DFT grids when heavy elements are included
!
      use modparallel, only : master
      use moddft, only : idftex, idftcor, nrad, nleb, hfexchange, bqrad
      use modatom, only : atomrad
      use modmolecule, only : natom, numatomic
      use modjob, only : method
      use modwarn, only : nwarn
      use modbasis, only : nao
      implicit none
      integer :: ii, maxelem
!
      do ii= 1,9
        atomrad(-ii)= bqrad(ii)
      enddo
!
      select case(method)
        case('B3LYP')
          idftex = 1
          idftcor= 1
          hfexchange= 0.2D+00
        case('B3LYP5')
          idftex = 1
          idftcor= 2
          hfexchange= 0.2D+00
        case('HARTREE-FOCK','MP2')
        case default
          if(master) then
            write(*,'(" Error! This program does not support method= ",a16,".")') method
            call iabort
          endif
      endselect
!
      if((idftex >= 1).or.(idftcor >= 1)) then
        maxelem= maxval(numatomic(1:natom))
        if(((maxelem >= 55).or.(nao >= 2000)).and.((nrad == 96).and.(nleb == 302))) then
          nwarn= nwarn+1
          if(master) write(*,'(" Warning! The number of DFT grids may not be enough.")')
        endif
      endif
!
      return
end

