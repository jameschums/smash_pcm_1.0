!-------------------------------------------------
  subroutine calchelfey(egrad,fulldmtrx,ish,jsh)
!-------------------------------------------------
!
! Driver of Helmann-Feynman gradient term
!
! In : fulldmtrx (density matrix)
!    : ish, jsh (shell indices)
! Inout : egrad (energy gradient value)
!
      use modparam, only : mxprsh
      use modmolecule, only : natom, coord, znuc
      use modbasis, only : locatom, locprim, locbf, mprim, mbf, mtype, ex, coeff, nao
      use modthresh, only : threshex
      implicit none
      integer,intent(in) :: ish, jsh
      integer :: nangij(2), nprimij(2), nbfij(2), locbfij(2), iatom, jatom
      integer :: iloc, jloc, iprim, jprim, i
      real(8),intent(in) :: fulldmtrx(nao,nao)
      real(8),intent(inout) :: egrad(3,natom)
      real(8) :: exij(mxprsh,2), cij(mxprsh,2), coordij(3,2)
      logical :: iandj
!
!      write(*,'("james calchelfey.F90 Driver of Helmann-Feynman gradient term")')
!
      iandj=(ish == jsh)
      nangij(1)= mtype(ish)
      nangij(2)= mtype(jsh)
      nprimij(1)= mprim(ish)
      nprimij(2)= mprim(jsh)
      nbfij(1)  = mbf(ish)
      nbfij(2)  = mbf(jsh)
      locbfij(1)= locbf(ish)
      locbfij(2)= locbf(jsh)
      iatom = locatom(ish)
      iloc  = locprim(ish)
      jatom = locatom(jsh)
      jloc  = locprim(jsh)
      do i= 1,3
        coordij(i,1)= coord(i,iatom)
        coordij(i,2)= coord(i,jatom)
      enddo
      do iprim= 1,nprimij(1)
        exij(iprim,1)= ex(iloc+iprim)
        cij(iprim,1) = coeff(iloc+iprim)
      enddo
      do jprim= 1,nprimij(2)
        exij(jprim,2)= ex(jloc+jprim)
        cij(jprim,2) = coeff(jloc+jprim)
      enddo
!
      if((nangij(1) <= 2).and.(nangij(2) <= 2)) then
        call int1cgmd(egrad,fulldmtrx,exij,cij,coordij,coord,znuc,natom,nao, &
&                     nprimij,nangij,nbfij,locbfij,mxprsh,threshex,iandj)
      elseif((nangij(1) <= 5).and.(nangij(2) <= 5)) then
        call int1grys(egrad,fulldmtrx,exij,cij,coordij,coord,znuc,natom,nao, &
&                     nprimij,nangij,nbfij,locbfij,mxprsh,threshex,iandj)
      else
        write(*,'(" Error! This program supports up to h function in helfey")')
        call iabort
      endif
!
      return
end


