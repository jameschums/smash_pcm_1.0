!-----------------------
 subroutine writegeom 
!-----------------------
!
! Write molecular geometry
!
      use modparallel, only : master
      use modmolecule, only : numatomic, natom, coord
      use modunit, only : toang
      use modenergy, only : escf, emp2, escsmp2
      use modjob, only : method
      use moddft, only : idftex, idftcor
      implicit none
      integer :: i, j ,jamesflag=0     
      real(8) :: jamese=0.0D+00
      logical :: file_exists, exist
      character(len=9) :: jamesfirst = 'first_run'
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
!!james
!      write(*,'("james line 305 writegeom.F90 jamesflag = ",i3," so far",/,&
! &               "hope it never prints zero mp2 ",f14.7," hartrees escf+escsmp2",/)')jamesflag,(escf+escsmp2)
      if (escf /= 0) then 
        write(*,'("energies scf(dft i think)",i3,2x,f14.7," a.u. escf")')jamesflag,(escf)
      endif
      if (emp2 /=0) then
        write(*,'("energies mp2 ",i3,2x,f14.7,"a.u. emp2")')jamesflag,(emp2)
      endif 
      if (escsmp2 /=0) then 
        write(*,'("energies mp2 ",i3,2x,f14.7," a.u. escsmp2")')jamesflag,(escsmp2)
      endif
 
      inquire(file="molden.xyz", exist=exist)
      if (exist) then
        open(13, file="molden.xyz", status="old", position="append", action="write")
      else
        open(13, file="molden.xyz", status="new", action="write")
      endif
      if (jamesflag >= 1 .and. (method == 'MP2')) then
        write(13, '(I3)')natom
        write(13, '(f14.7)')(escsmp2)
        do i= 1,natom
          write(13,'(3x,a3,3x,3f14.7)')table(numatomic(i)),(coord(j,i)*toang,j=1,3)
        end do
      elseif (jamesflag >= 1 .and. escf /=0 .and. ((idftex >= 1).or.(idftcor >= 1))) then
        write(13, '(I3)')natom
        write(13, '(f14.7)')(escf)
        do i= 1,natom
          write(13,'(3x,a3,3x,3f14.7)')table(numatomic(i)),(coord(j,i)*toang,j=1,3)
        end do
      endif
      close(13)
      jamesflag = ( jamesflag +1 )
! ----   
      if(master) then
        write(*,'("james writegeom.F90 line 62 = ",i3," so far, hope it never prints zero ")')jamesflag
        write(*,'(" ----------------------------------------------------")')
        write(*,'("          Molecular Geometry (Angstrom)")')
        write(*,'("  Atom            X             Y             Z")')
        write(*,'(" jgeomi----------------------------------------------------")')
        write(*, '(I3)')natom
        write(*, '(f14.7)')(escf+emp2)
        do i= 1,natom
          write(*,'(3x,a3,3x,3f14.7)')table(numatomic(i)),(coord(j,i)*toang,j=1,3)
        end do
        write(*,'(" jgeomj----------------------------------------------------",/)')
      endif
      return
end

