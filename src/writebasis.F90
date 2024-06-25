!------------------------
  subroutine writebasis
!------------------------
!
! Write basis functions
!
      use modparallel, only : master
      use modmolecule, only : numatomic
      use modbasis, only : nshell, ex, coeffinp, locprim, locatom, mprim, mtype
      use modprint, only : iprint
      implicit none
      integer :: iatom, ishell, iloc, iprim, jatomcheck(-9:112)=0
      logical :: second
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
! Write basis functions
!
      write(*,'("james writebasis.F90 needs iprint4 for printing")')
      second=.false.
      if(master.and.(iprint >= 4)) then
        iatom= 0
        write(*,'(" -------------")')
        write(*,'("   Basis set")')
        write(*,'(" -------------")')
        do ishell= 1, nshell
          if(iatom /= locatom(ishell)) then
            if(jatomcheck(numatomic(locatom(ishell))) /= 0)cycle
            jatomcheck(numatomic(locatom(ishell)))= 1
            if(second) write(*,'("  ****")')
            second=.true.
            write(*,'(2x,a3)')table(numatomic(locatom(ishell)))
            iatom= locatom(ishell)
          endif
          if(mtype(ishell) == 0) write(*,'(4x,"S",i3)') mprim(ishell)
          if(mtype(ishell) == 1) write(*,'(4x,"P",i3)') mprim(ishell)
          if(mtype(ishell) == 2) write(*,'(4x,"D",i3)') mprim(ishell)
          if(mtype(ishell) == 3) write(*,'(4x,"F",i3)') mprim(ishell)
          if(mtype(ishell) == 4) write(*,'(4x,"G",i3)') mprim(ishell)
          if(mtype(ishell) == 5) write(*,'(4x,"H",i3)') mprim(ishell)
          if(mtype(ishell) == 6) write(*,'(4x,"I",i3)') mprim(ishell)
!
          if(mtype(ishell) >  6) then
            write(*,'(" Error! The subroutine writebasis supports up to i functions.")')
            call iabort
          endif 
!
          iloc= locprim(ishell)
          do iprim= 1,mprim(ishell)
             write(*,'(3x,f16.7,1x,f15.8)') ex(iloc+iprim), coeffinp(iloc+iprim)
          enddo
        enddo
        write(*,'("  ****")')
      endif
      return
end


