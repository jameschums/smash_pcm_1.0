!------------------------------------------------------------------------------------
  subroutine calcnewcoord(coord,coordold,egrad,egradold,ehess,displc,natom3,iopt, &
&                         nproc,myrank,mpi_comm)
!------------------------------------------------------------------------------------
!
! Calculate new Cartesian coordinate with gradient and hessian
!
      use modparallel, only : master
      use modprint, only : iprint
      use modunit, only : toang
      use modmolecule, only : numatomic
      implicit none
      integer,intent(in) :: natom3, iopt, nproc, myrank, mpi_comm
      integer :: i, j, ii
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, third=0.3333333333333333D+00
      real(8),intent(inout) :: egrad(natom3), egradold(natom3), ehess(natom3*(natom3+1)/2)
      real(8),intent(inout) :: coord(natom3), coordold(natom3), displc(natom3*3)
!
      real(8) :: work(natom3,natom3),eigen(natom3),work2(natom3,natom3),work3(natom3,natom3)
!
      character(len=3) :: table(-9:112)= &
&     (/'Bq9','Bq8','Bq7','Bq6','Bq5','Bq4','Bq3','Bq2','Bq ','X  ',&
&       'H  ','He ','Li ','Be ','B  ','C  ','N  ','O  ','F  ','Ne ','Na ','Mg ','Al ','Si ','P  ',&
&       'S  ','Cl ','Ar ','K  ','Ca ','Sc ','Ti ','V  ','Cr ','Mn ','Fe ','Co ','Ni ','Cu ','Zn ',&
&       'Ga ','Ge ','As ','Se ','Br ','Kr ','Rb ','Sr ','Y  ','Zr ','Nb ','Mo ','Tc ','Ru ','Rh ',&
&       'Pd ','Ag ','Cd ','In ','Sn ','Sb ','Te ','I  ','Xe ','Cs ','Ba ','La ','Ce ','Pr ','Nd ',&
&       'Pm ','Sm ','Eu ','Gd ','Tb ','Dy ','Ho ','Er ','Tm ','Yb ','Lu ','Hf ','Ta ','W  ','Re ',&
&       'Os ','Ir ','Pt ','Au ','Hg ','Tl ','Pb ','Bi ','Po ','At ','Rn ','Fr ','Ra ','Ac ','Th ',&
&       'Pa ','U  ','Np ','Pu ','Am ','Cm ','Bk ','Cf ','Es ','Fm ','Md ','No ','Lr ','Rf ','Db ',&
&       'Sg ','Bh ','Hs ','Mt ','Uun','Uuu','Uub'/)
!
      if(iopt == 1) then
!
! Set initial hessian
!
        ehess(:)= zero
        do i= 1,natom3
          ehess(i*(i+1)/2)= third
        enddo 
      else
!
! Updata hessian matrix
!
        call hessianbfgs(ehess,coord,coordold,egrad,egradold,displc,natom3)
      endif
!
      do i=1,natom3
        ii= i*(i-1)/2
        do j=1,i
          work(j,i)=ehess(ii+j)
        enddo
      enddo
      call diag('V','U',natom3,work,natom3,eigen,nproc,myrank,mpi_comm)
!
      do i=1,natom3
        eigen(i)= one/eigen(i)
      enddo
      do i=1,natom3
        do j=1,natom3
          work2(j,i)=work(j,i)*eigen(i)
        enddo
      enddo
      call dgemm('N','T',natom3,natom3,natom3,one,work,natom3,work2,natom3,zero,work3,natom3)
!
! Copy old coordinate and gradient
!
      coordold(:)= coord(:)
      egradold(:)= egrad(:)
!
      do i=1,natom3
        do j=1,natom3
          coord(i)=coord(i)-work3(i,j)*egrad(j)
        enddo
      enddo
!
! Print delta xyz
!
      if(master.and.(iprint >= 2)) then
        write(*,'(" ----------------------------------------------------")')
        write(*,'("          Delta xyz (Angstrom)")')
        write(*,'("  Atom            X             Y             Z")')
        write(*,'(" ----------------------------------------------------")')
        do i= 1,natom3/3
          write(*,'(3x,a3,3x,3f14.7)')table(numatomic(i)), &
&              ((coord((i-1)*3+j)-coordold((i-1)*3+j))*toang,j=1,3)
        enddo
        write(*,'(" ----------------------------------------------------")')
      endif
!
      return
end
