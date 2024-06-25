!---------------------------------------------
  subroutine setnextopt(coordold,natom,iopt)
!---------------------------------------------
!
! Set parameters for next optimization step
!
      use modbasis, only : ex, coeff, nshell, nao, nprim, locprim, locbf, &
&                          locatom, mprim, mbf, mtype, spher
      use modguess, only : ex_g, coeff_g, nshell_g, nao_g, nmo_g, nprim_g, locprim_g, locbf_g, &
&                          locatom_g, mprim_g, mbf_g, mtype_g, spher_g, coord_g, guess
      use modmolecule, only : nmo
      implicit none
      integer,intent(in) :: natom, iopt
      integer :: iprim, ishell, iatom
      real(8),intent(in) :: coordold(3,natom)
!
! Set MO projection as initial MO calculation
!
      guess= 'UPDATE'
!
! Copy coordinate and energy gradient
!
      do iatom= 1,natom
        coord_g(1,iatom)= coordold(1,iatom)
        coord_g(2,iatom)= coordold(2,iatom)
        coord_g(3,iatom)= coordold(3,iatom)
      enddo
!
! Copy basis set information
!
      nmo_g= nmo
!
      if(iopt == 1) then
        nao_g= nao
        nprim_g= nprim
        nshell_g= nshell
!
        do iprim= 1,nprim
          ex_g(iprim)= ex(iprim)
          coeff_g(iprim)= coeff(iprim)
        enddo
!
        do ishell= 1,nshell
          locprim_g(ishell)= locprim(ishell)
          locbf_g(ishell)  = locbf(ishell)
          locatom_g(ishell)= locatom(ishell)
          mprim_g(ishell)  = mprim(ishell)
          mbf_g(ishell)    = mbf(ishell)
          mtype_g(ishell)  = mtype(ishell)
        enddo
!
        spher_g= spher
      endif
!
      return
end


