!--------------------------------------------------------------------
  subroutine gradoneei(egrad,egrad1,fulldmtrx,ewdmtrx,nproc,myrank)
!--------------------------------------------------------------------
!
! Driver of derivatives for one-electron and overlap integrals
!
      use modparallel, only : master
      use modbasis, only : nshell, nao, mtype
      use modmolecule, only : natom
      use modecp, only : flagecp
      implicit none
      integer,intent(in) :: nproc, myrank
      integer :: ish, jsh, len1, i, maxfunc(0:6), maxbasis
      real(8),parameter :: zero=0.0D+00, two=2.0D+00
      real(8),intent(in) :: fulldmtrx(nao*nao), ewdmtrx(nao*(nao+1)/2)
      real(8),intent(out) :: egrad1(3*natom)
      real(8),intent(inout) :: egrad(3*natom)
      data maxfunc/1,3,6,10,15,21,28/
!
      maxbasis= maxval(mtype(1:nshell))
      if(maxbasis > 6) then
        if(master) write(*,'(" Error! This program supports up to h function in gradoneei")')
        call iabort
      endif
      len1= maxfunc(maxbasis+1)
      egrad1(:)= zero
!
!$OMP parallel reduction(+:egrad1)
      do ish= nshell-myrank,1,-nproc
!$OMP do
        do jsh= 1,ish
          call calcdoverlap(egrad1,ewdmtrx,ish,jsh)
          call calchelfey(egrad1,fulldmtrx,ish,jsh)
        enddo
!$OMP enddo
      enddo
!
      do ish= myrank+1,nshell,nproc
!$OMP do
        do jsh= 1,nshell
          call calcdkinetic(egrad1,fulldmtrx,ish,jsh)
          call calcdcoulomb(egrad1,fulldmtrx,ish,jsh,len1)
        enddo
!$OMP enddo
      enddo
!$OMP end parallel
!
      do i= 1,3*natom
        egrad(i)= egrad(i)+egrad1(i)*two
      enddo
!
! Add ECP derivative terms
!
      if(flagecp) then
        egrad1(:)= zero
        call gradoneeiecp(egrad1,fulldmtrx,nproc,myrank)
        do i= 1,3*natom
          egrad(i)= egrad(i)+egrad1(i)*two
        enddo
      endif
! 
      return
end

