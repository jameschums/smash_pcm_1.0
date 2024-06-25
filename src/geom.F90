! Copyright 2014-2017  Kazuya Ishimura
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.
!
!-----------------------
  subroutine nucenergy
!-----------------------
!
! Calculate nuclear replusion energy
!
      use modparallel, only : master
      use modmolecule, only : natom, coord, znuc
      use modthresh, only : threshatom
      use modenergy, only : enuc
      use modwarn, only : nwarn
      implicit none
      integer :: iatom, jatom
      real(8),parameter :: zero=0.0D+00
      real(8) :: xyz(3), rr, chrgij
!
      enuc= zero

      do iatom= 2,natom
        do jatom= 1,iatom-1
          xyz(1)= coord(1,iatom)-coord(1,jatom)
          xyz(2)= coord(2,iatom)-coord(2,jatom)
          xyz(3)= coord(3,iatom)-coord(3,jatom)
          rr= sqrt(xyz(1)*xyz(1)+xyz(2)*xyz(2)+xyz(3)*xyz(3))
          chrgij= znuc(iatom)*znuc(jatom)
          if(rr /= zero) then
            enuc= enuc+chrgij/rr
            if((rr <= threshatom).and.master) then
              write(*,'("Warning! Distance of Atoms",i4," and",i4," is short!")') iatom, jatom
              nwarn= nwarn+1
            endif       
          else
            if((chrgij /= zero).and.master) then
              write(*,'("Error! Atoms",i4," and",i4," are the same position!")') iatom, jatom
              call iabort
            endif
          endif
        enddo
      enddo
      return
end


!---------------------------------------------
  subroutine nucgradient(egrad,nproc,myrank)
!---------------------------------------------
!
! Calculate gradinet of nuclear replusion energy
!
      use modmolecule, only : natom, coord, znuc
      implicit none
      integer,intent(in) :: nproc, myrank
      integer :: iatom, jatom, i
      real(8),intent(inout) :: egrad(3,natom)
      real(8) :: xyz(3), rr, chrgij
!
      do iatom= 2+myrank,natom,nproc
        do jatom= 1,iatom-1
          xyz(1)= coord(1,iatom)-coord(1,jatom)
          xyz(2)= coord(2,iatom)-coord(2,jatom)
          xyz(3)= coord(3,iatom)-coord(3,jatom)
          rr= xyz(1)*xyz(1)+xyz(2)*xyz(2)+xyz(3)*xyz(3)
          rr= rr*sqrt(rr)
          chrgij= znuc(iatom)*znuc(jatom)
          do i= 1,3
            egrad(i,iatom)= egrad(i,iatom)-xyz(i)*chrgij/rr
            egrad(i,jatom)= egrad(i,jatom)+xyz(i)*chrgij/rr
          enddo
        enddo
      enddo
      return
end


!-----------------------------------------------------------------------------------
  subroutine setredundantcoord(iredun,isizered,numbond,numangle,numtorsion,exceed)
!-----------------------------------------------------------------------------------
!
! Set redundant internal coordinate
! Covalent radii (H - Cn): P. Pyykko, M. Atsumi, Chem. Eur. J., 186 (2009) 15.
!
      use modparallel, only : master
      use modmolecule, only : coord, natom, numatomic
      use modunit, only : tobohr
      implicit none
      integer,parameter :: maxconnect=13
      integer,intent(in) :: isizered
      integer,intent(out) :: iredun(4,isizered/4), numbond, numangle, numtorsion
      integer :: numredun, maxsize, ijpair(natom,maxconnect), iatom, jatom, katom, icount
      integer :: ibond, kpair, iangle, jangle, numb, iblock(natom), ib, jb, ijatom(2) 
      real(8),parameter :: zero=0.0D+00
      real(8) :: radii(112), rrij, thresh, rrijmin
      logical,intent(out) :: exceed
      data radii/ &
&     0.32D+00, 0.46D+00, 1.33D+00, 1.02D+00, 0.85D+00, 0.75D+00, 0.71D+00, 0.63D+00, 0.64D+00, &
&     0.67D+00, 1.55D+00, 1.39D+00, 1.26D+00, 1.16D+00, 1.11D+00, 1.03D+00, 0.99D+00, 0.96D+00, &
&     1.96D+00, 1.71D+00, 1.48D+00, 1.36D+00, 1.34D+00, 1.22D+00, 1.19D+00, 1.16D+00, 1.11D+00, &
&     1.10D+00, 1.12D+00, 1.18D+00, 1.24D+00, 1.21D+00, 1.21D+00, 1.16D+00, 1.14D+00, 1.17D+00, &
&     2.10D+00, 1.85D+00, 1.63D+00, 1.54D+00, 1.47D+00, 1.38D+00, 1.28D+00, 1.25D+00, 1.25D+00, &
&     1.20D+00, 1.28D+00, 1.36D+00, 1.42D+00, 1.40D+00, 1.40D+00, 1.36D+00, 1.33D+00, 1.31D+00, &
&     2.32D+00, 1.96D+00, 1.80D+00, 1.63D+00, 1.76D+00, 1.74D+00, 1.73D+00, 1.72D+00, 1.68D+00, &
&     1.69D+00, 1.68D+00, 1.67D+00, 1.66D+00, 1.65D+00, 1.64D+00, 1.70D+00, 1.62D+00, 1.52D+00, &
&     1.46D+00, 1.37D+00, 1.31D+00, 1.29D+00, 1.22D+00, 1.23D+00, 1.24D+00, 1.33D+00, 1.44D+00, &
&     1.44D+00, 1.51D+00, 1.45D+00, 1.47D+00, 1.42D+00, 2.23D+00, 2.01D+00, 1.86D+00, 1.75D+00, &
&     1.69D+00, 1.70D+00, 1.71D+00, 1.72D+00, 1.66D+00, 1.66D+00, 1.68D+00, 1.68D+00, 1.65D+00, &
&     1.67D+00, 1.73D+00, 1.76D+00, 1.61D+00, 1.57D+00, 1.49D+00, 1.43D+00, 1.41D+00, 1.34D+00, &
&     1.29D+00, 1.28D+00, 1.21D+00, 1.22D+00/
!
      ijpair(:,:)= 0
      numredun= 0
      maxsize= isizered/4
      exceed=.false.
!
! Set bond strech
!
      numbond= 0
      do iatom= 1,natom
        if(numatomic(iatom) > 112) then
          if(master) &
&         write(*,'(" Error! This program supports up to Cn in Subroutine setredundantcoord.")')
          call iabort
        elseif(numatomic(iatom) < 1) then
          if(master) &
&         write(*,'(" Error! This program does not support dummy and ghost atoms ", &
&                   "in Subroutine setredundantcoord currently.")')
          call iabort
        endif
        icount= 0
        do jatom= 1,natom
          if(jatom == iatom) cycle
          thresh=(1.25D+00*tobohr*(radii(numatomic(iatom))+radii(numatomic(jatom))))**2
          rrij= (coord(1,jatom)-coord(1,iatom))**2+(coord(2,jatom)-coord(2,iatom))**2 &
&              +(coord(3,jatom)-coord(3,iatom))**2
          if(rrij <= thresh) then
            icount= icount+1
            if(icount > maxconnect) then
              if(master) write(*,'(" Error! There are too many atoms near Atom",i4,".")')iatom
              call iabort
            endif
            ijpair(iatom,icount)= jatom
            if(jatom > iatom) then
              numredun= numredun+1
              if(numredun > maxsize) then
                exceed=.true.
                return
              endif
              numbond= numbond+1
              iredun(1,numbond)= iatom
              iredun(2,numbond)= jatom
            endif
          endif
        enddo
      enddo
!
! Check the number of molecular fragments
!
      numb= 0
      iblock(:)= 0
 mblock:do iatom= 1,natom
          if(iblock(iatom) /= 0) cycle
          numb= numb+1
          iblock(iatom)= numb
          call checkbond(iatom,numb,iblock,ijpair,natom,maxconnect)
          do jatom= iatom+1,natom
            if(iblock(jatom) == 0) cycle mblock
          enddo
          exit mblock
        enddo mblock
!
     if(master.and.(numb /= 1)) then
       write(*,'(" There are",i3," moleclar blocks in redundant coordinate.")') numb
     endif
!
! Calculate length between molecular fragments
!
      do ib= 1, numb
        do jb= ib+1,numb
          rrijmin= zero
          do iatom= 1,natom
            if(iblock(iatom) /= ib) cycle
            do jatom= iatom+1,natom
              if(iblock(jatom) /= jb) cycle
                rrij= (coord(1,jatom)-coord(1,iatom))**2+(coord(2,jatom)-coord(2,iatom))**2 &
&                    +(coord(3,jatom)-coord(3,iatom))**2
              if((rrijmin == zero).or.(rrij < rrijmin)) then
                rrijmin= rrij
                ijatom(1)= min(iatom,jatom)
                ijatom(2)= max(iatom,jatom)
              endif
            enddo
          enddo
          do icount= 1,maxconnect
            if(ijpair(ijatom(1),icount) == 0) then
              ijpair(ijatom(1),icount)= ijatom(2)
              exit
            endif
            if(icount == maxconnect) then
              if(master) write(*,'(" Error! There are too many atoms near Atom",i4,".")')ijatom(1)
              call iabort
            endif
          enddo
          do icount= 1,maxconnect
            if(ijpair(ijatom(2),icount) == 0) then
              ijpair(ijatom(2),icount)= ijatom(1)
              exit
            endif
            if(icount == maxconnect) then
              if(master) write(*,'(" Error! There are too many atoms near Atom",i4,".")')ijatom(2)
              call iabort
            endif
          enddo
          numredun= numredun+1
          if(numredun > maxsize) then
            exceed=.true.
            return
          endif
          numbond= numbond+1
          iredun(1,numbond)= ijatom(1)
          iredun(2,numbond)= ijatom(2)
        enddo
      enddo
             
!
! Set bond angle
!
      numangle= 0
      do ibond= 1,numbond
        iatom= iredun(1,ibond)
        jatom= iredun(2,ibond)
        do kpair= 1,maxconnect
          katom= ijpair(jatom,kpair)
          if(katom == 0) exit
          if(katom > iatom) then
            numredun= numredun+1
            if(numredun > maxsize) then
              exceed=.true.
              return
            endif
            numangle= numangle+1
            iredun(1,numbond+numangle)= iatom
            iredun(2,numbond+numangle)= jatom
            iredun(3,numbond+numangle)= katom
          endif
        enddo
        iatom= iredun(2,ibond)
        jatom= iredun(1,ibond)
        do kpair= 1,maxconnect
          katom= ijpair(jatom,kpair)
          if(katom == 0) exit
          if(katom > iatom) then
            numredun= numredun+1
            if(numredun > maxsize) then
              exceed=.true.
              return
            endif
            numangle= numangle+1
            iredun(1,numbond+numangle)= iatom
            iredun(2,numbond+numangle)= jatom
            iredun(3,numbond+numangle)= katom
          endif
        enddo
      enddo
!
!
! Set dihedral angle
!
      numtorsion= 0
      do iangle= numbond+1,numbond+numangle
        iatom= iredun(1,iangle)
        jatom= iredun(2,iangle)
        katom= iredun(3,iangle)
        do jangle= iangle+1,numbond+numangle
          if(iredun(1,jangle) == jatom) then
            if(iredun(2,jangle) == katom) then
              if(iredun(3,jangle) /= iatom) then
                numredun= numredun+1
                if(numredun > maxsize) then
                  exceed=.true.
                  return
                endif
                numtorsion= numtorsion+1
                iredun(1,numbond+numangle+numtorsion)= iatom
                iredun(2,numbond+numangle+numtorsion)= jatom
                iredun(3,numbond+numangle+numtorsion)= katom
                iredun(4,numbond+numangle+numtorsion)= iredun(3,jangle)
              endif
            endif
          endif
          if(iredun(3,jangle) == jatom) then
            if(iredun(2,jangle) == katom) then
              if(iredun(1,jangle) /= iatom) then
                numredun= numredun+1
                if(numredun > maxsize) then
                  exceed=.true.
                  return
                endif
                numtorsion= numtorsion+1
                iredun(1,numbond+numangle+numtorsion)= iatom
                iredun(2,numbond+numangle+numtorsion)= jatom
                iredun(3,numbond+numangle+numtorsion)= katom
                iredun(4,numbond+numangle+numtorsion)= iredun(1,jangle)
              endif
            endif
          endif
        enddo
        iatom= iredun(3,iangle)
        jatom= iredun(2,iangle)
        katom= iredun(1,iangle)
        do jangle= iangle+1,numbond+numangle
          if(iredun(1,jangle) == jatom) then
            if(iredun(2,jangle) == katom) then
              if(iredun(3,jangle) /= iatom) then
                numredun= numredun+1
                if(numredun > maxsize) then
                  exceed=.true.
                  return
                endif
                numtorsion= numtorsion+1
                iredun(1,numbond+numangle+numtorsion)= iatom
                iredun(2,numbond+numangle+numtorsion)= jatom
                iredun(3,numbond+numangle+numtorsion)= katom
                iredun(4,numbond+numangle+numtorsion)= iredun(3,jangle)
              endif
            endif
          endif
          if(iredun(3,jangle) == jatom) then
            if(iredun(2,jangle) == katom) then
              if(iredun(1,jangle) /= iatom) then
                numredun= numredun+1
                if(numredun > maxsize) then
                  exceed=.true.
                  return
                endif
                numtorsion= numtorsion+1
                iredun(1,numbond+numangle+numtorsion)= iatom
                iredun(2,numbond+numangle+numtorsion)= jatom
                iredun(3,numbond+numangle+numtorsion)= katom
                iredun(4,numbond+numangle+numtorsion)= iredun(1,jangle)
              endif
            endif
          endif
        enddo
      enddo
!
      if((natom >= 4).and.(numtorsion == 0)) then
        do iangle= numbond+1,numbond+numangle
          iatom= iredun(1,iangle)
          jatom= iredun(2,iangle)
          katom= iredun(3,iangle)
          do jangle= iangle+1,numbond+numangle
            if(iredun(1,jangle) == katom) then
              if(iredun(2,jangle) == jatom) then
                numredun= numredun+1
                if(numredun > maxsize) then
                  exceed=.true.
                  return
                endif
                numtorsion= numtorsion+1
                iredun(1,numbond+numangle+numtorsion)= iatom
                iredun(2,numbond+numangle+numtorsion)= jatom
                iredun(3,numbond+numangle+numtorsion)= katom
                iredun(4,numbond+numangle+numtorsion)= iredun(3,jangle)
              endif
            endif
          enddo
        enddo
      endif
!
      return
end


!----------------------------------------------------------------------------
  recursive subroutine checkbond(iatom,numb,iblock,ijpair,natom,maxconnect)
!----------------------------------------------------------------------------
!
! Check the number of molecular fragments
!
      implicit none
      integer,intent(in) :: iatom, numb, natom, maxconnect, ijpair(natom,maxconnect)
      integer,intent(inout) :: iblock(natom)
      integer :: icount, jatom
!
      do icount= 1,maxconnect
        jatom= ijpair(iatom,icount)
        if(jatom == 0) exit
        if(iblock(jatom) /= 0) cycle
        iblock(jatom)= numb
        call checkbond(jatom,numb,iblock,ijpair,natom,maxconnect)
      enddo
!
      return
end
