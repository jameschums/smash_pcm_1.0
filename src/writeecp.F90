!----------------------
  subroutine writeecp
!----------------------
!
! Write ECP functions
!
      use modparallel, only : master
      use modmolecule, only : numatomic, natom
      use modecp, only : maxangecp, mtypeecp, locecp, mprimecp, execp, coeffecp, izcore
      use modprint, only : iprint
      implicit none
      integer :: iatom, ll, jprim, jloc, k, nprim, jatomcheck(-9:112)=0
      character(len=7) :: tblecp
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
! Write ECP functions
!
      if(master.and.(iprint >= 1)) then
        write(*,'(" ----------------")')
        write(*,'("   ECP function")')
        write(*,'(" ----------------")')
!
        do iatom= 1,natom
          if(maxangecp(iatom) /= -1)then
            if(jatomcheck(numatomic(iatom)) == 0) then
              jatomcheck(numatomic(iatom))= iatom
              write(*,'(2x,a3)')table(numatomic(iatom))
              tblecp=table(numatomic(iatom))
              ll= len_trim(tblecp)
              tblecp= tblecp(1:ll)//'-ECP'
              write(*,'(2x,a7,2x,i3,2x,i3)')tblecp, maxangecp(iatom), izcore(iatom)
              nprim= mprimecp(0,iatom)
              if(maxangecp(iatom) == 0) write(*,'(2x,"s-ul potential")')
              if(maxangecp(iatom) == 1) write(*,'(2x,"p-ul potential")')
              if(maxangecp(iatom) == 2) write(*,'(2x,"d-ul potential")')
              if(maxangecp(iatom) == 3) write(*,'(2x,"f-ul potential")')
              if(maxangecp(iatom) == 4) write(*,'(2x,"g-ul potential")')
              write(*,'(3x,i2)')nprim
              do jprim= 1,nprim
                jloc= jprim+locecp(0,iatom)
                write(*,'(2x,i1,2f16.7)')mtypeecp(jloc), execp(jloc), coeffecp(jloc)
              enddo
              do k= 1,maxangecp(iatom)
                nprim= mprimecp(k,iatom)
                if(k == 1) write(*,'(2x,"s-ul potential")')
                if(k == 2) write(*,'(2x,"p-ul potential")')
                if(k == 3) write(*,'(2x,"d-ul potential")')
                if(k == 4) write(*,'(2x,"f-ul potential")')
                write(*,'(3x,i2)')nprim
                do jprim= 1,nprim
                  jloc= jprim+locecp(k,iatom)
                  write(*,'(2x,i1,2f16.7)')mtypeecp(jloc), execp(jloc), coeffecp(jloc)
                enddo
              enddo
            endif
          endif
        enddo
!
        write(*,*)
      endif
!
      return
end

