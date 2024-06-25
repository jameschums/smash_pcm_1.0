!-----------------------------------------------------------------------------------
  subroutine calcbmatrix(coordredun,coord,bmat,iredun,numbond,numangle,numtorsion)
!-----------------------------------------------------------------------------------
!
! Calculate transformation matrix(B-matrix) from Cartesian to internal coordinate
!
      use modmolecule, only : natom
      implicit none
      integer,intent(in) :: numbond, numangle, numtorsion, iredun(4,numbond+numangle+numtorsion)
      integer :: iatom, jatom, katom, latom, ii, jj, kk, ll, mm, ibond, iangle, itorsion
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, two=2.0D+00, four=4.0D+00
      real(8),intent(in) :: coord(3,natom)
      real(8),intent(out) :: coordredun(numbond+numangle+numtorsion)
      real(8),intent(out) :: bmat(numbond+numangle+numtorsion,natom*3)
      real(8) :: rrij, rrji, rrjk, rji, rjk, rrkl, rij, rkl, pi, pi2
      real(8) :: dotj, dotk, dotjk, dotjk3, cp1(3), cp2(3), cp3(3), sinj, sink, b1, b2, b3
      real(8) :: xyzji(3), xyzjk(3), xyzij(3), xyzkl(3)
      real(8) :: unitji(3), unitjk(3), unitij(3), unitkl(3)
!
      pi= atan(one)*four
      pi2=pi*two
      bmat(:,:)= zero
!
! Bond strech
!
      do ibond= 1,numbond
        iatom= iredun(1,ibond)
        jatom= iredun(2,ibond)
        rrij= (coord(1,jatom)-coord(1,iatom))**2+(coord(2,jatom)-coord(2,iatom))**2 &
&            +(coord(3,jatom)-coord(3,iatom))**2
        rij= sqrt(rrij)
        coordredun(ibond)= rij
!
        ii=(iatom-1)*3
        jj=(jatom-1)*3
        do mm= 1,3
          bmat(ibond,ii+mm)=-(coord(mm,jatom)-coord(mm,iatom))/rij
          bmat(ibond,jj+mm)= (coord(mm,jatom)-coord(mm,iatom))/rij
        enddo
      enddo
!
! Bond angle
!
      do iangle= numbond+1,numbond+numangle
        iatom= iredun(1,iangle)
        jatom= iredun(2,iangle)
        katom= iredun(3,iangle)
        rrji= zero
        rrjk= zero
        do mm= 1,3
          xyzji(mm)= coord(mm,iatom)-coord(mm,jatom)
          xyzjk(mm)= coord(mm,katom)-coord(mm,jatom)
          rrji= rrji+xyzji(mm)*xyzji(mm)
          rrjk= rrjk+xyzjk(mm)*xyzjk(mm)
        enddo
        rji= sqrt(rrji)
        rjk= sqrt(rrjk)
        dotj= zero
        do mm= 1,3
          unitji(mm)= xyzji(mm)/rji
          unitjk(mm)= xyzjk(mm)/rjk
          dotj= dotj+unitji(mm)*unitjk(mm)
        enddo
        if(abs(dotj) >= one) then
          write(*,'(" Error! During calculation of bond angles in calcbmatrix.")')
          write(*,'(" Use Cartesian coordinate. The input is")')
          write(*,'("   opt cartesian=.true.",/)')
          call iabort
        endif
        coordredun(iangle)= acos(dotj)
!
        sinj= sqrt(one-dotj*dotj)
        ii=(iatom-1)*3
        jj=(jatom-1)*3
        kk=(katom-1)*3
        do mm= 1,3
          b1=(dotj*unitji(mm)-unitjk(mm))/(rji*sinj)
          b2=(dotj*unitjk(mm)-unitji(mm))/(rjk*sinj)
          bmat(iangle,ii+mm)= b1
          bmat(iangle,kk+mm)= b2
          bmat(iangle,jj+mm)=-(b1+b2)
        enddo
      enddo
!
! Set dihedral angle
!
      do itorsion= numbond+numangle+1,numbond+numangle+numtorsion
        iatom= iredun(1,itorsion)
        jatom= iredun(2,itorsion)
        katom= iredun(3,itorsion)
        latom= iredun(4,itorsion)
        rrij= zero
        rrjk= zero
        rrkl= zero
        do mm= 1,3
          xyzij(mm)= coord(mm,jatom)-coord(mm,iatom)
          xyzjk(mm)= coord(mm,katom)-coord(mm,jatom)
          xyzkl(mm)= coord(mm,latom)-coord(mm,katom)
          rrij= rrij+xyzij(mm)*xyzij(mm)
          rrjk= rrjk+xyzjk(mm)*xyzjk(mm)
          rrkl= rrkl+xyzkl(mm)*xyzkl(mm)
        enddo
        rij= sqrt(rrij)
        rjk= sqrt(rrjk)
        rkl= sqrt(rrkl)
        dotj= zero
        dotk= zero
        do mm= 1,3
          unitij(mm)= xyzij(mm)/rij
          unitjk(mm)= xyzjk(mm)/rjk
          unitkl(mm)= xyzkl(mm)/rkl
          dotj= dotj-unitij(mm)*unitjk(mm)
          dotk= dotk-unitjk(mm)*unitkl(mm)
        enddo
        cp1(1)= unitij(2)*unitjk(3)-unitij(3)*unitjk(2)
        cp1(2)= unitij(3)*unitjk(1)-unitij(1)*unitjk(3)
        cp1(3)= unitij(1)*unitjk(2)-unitij(2)*unitjk(1)
        cp2(1)= unitjk(2)*unitkl(3)-unitjk(3)*unitkl(2)
        cp2(2)= unitjk(3)*unitkl(1)-unitjk(1)*unitkl(3)
        cp2(3)= unitjk(1)*unitkl(2)-unitjk(2)*unitkl(1)
        cp3(1)= cp1(2)*cp2(3)-cp1(3)*cp2(2)
        cp3(2)= cp1(3)*cp2(1)-cp1(1)*cp2(3)
        cp3(3)= cp1(1)*cp2(2)-cp1(2)*cp2(1)
        if((abs(dotj) >= one).or.(abs(dotk) >= one)) then
          write(*,'(" Error! During calculation of torsion angles in calcbmatrix.")')
          call iabort
        endif
!
        sinj= sqrt(one-dotj*dotj)
        sink= sqrt(one-dotk*dotk)
        dotjk= zero
        dotjk3=zero
        do mm= 1,3
          dotjk= dotjk+cp1(mm)*cp2(mm)/(sinj*sink)
          dotjk3=dotjk3+unitjk(mm)*cp3(mm)
        enddo
        if(abs(dotjk) > one) dotjk= sign(one,dotjk)
        coordredun(itorsion)= acos(dotjk)
        if(dotjk3 < zero) coordredun(itorsion)=-coordredun(itorsion)
!
        if(coordredun(itorsion) > pi) coordredun(itorsion)= coordredun(itorsion)-pi2
        if(coordredun(itorsion) <-pi) coordredun(itorsion)= coordredun(itorsion)+pi2
!
        ii=(iatom-1)*3
        jj=(jatom-1)*3
        kk=(katom-1)*3
        ll=(latom-1)*3
        do mm= 1,3
          b1=-cp1(mm)/(rij*sinj*sinj)
          b2= cp1(mm)*(rjk-rij*dotj)/(rjk*rij*sinj*sinj)-dotk*cp2(mm)/(rjk*sink*sink)
          b3= cp2(mm)/(rkl*sink*sink)
          bmat(itorsion,ii+mm)= b1
          bmat(itorsion,jj+mm)= b2
          bmat(itorsion,ll+mm)= b3
          bmat(itorsion,kk+mm)=-(b1+b2+b3)
        enddo
      enddo
!
      return
end


