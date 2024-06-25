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


!-------------------------
  subroutine opendatfile
!-------------------------
!
! Open temporary file
!
      use modparallel, only : master
      use modiofile, only : input
      implicit none
      integer(selected_int_kind(9)) :: getpid, iprocess
      character(len=30) :: filename
!
      if(master) then
        iprocess= getpid()
        write(filename,*) iprocess
        filename= "input.dat"//adjustl(filename)
        open(unit=input,file=filename,status='replace')
      endif
!
      return
end


!---------------------------
  subroutine opencheckfile
!---------------------------
!
! Open checkpoint file
!
      use modparallel, only : master
      use modiofile, only : icheck, check
      implicit none
!
      if(master) then
        open(unit=icheck,file=check,form='unformatted',status='unknown')
      endif
!
      return
end


!--------------------
  subroutine iabort
!--------------------
!
! Abort the calculation
!
      use modparallel, only : master
      implicit none
!
      if(master) then
        write(*,'("Pants!! :( Calculation finished abnormally.")')
        call tstamp(2)
      else
        call sleep(5)
      endif  
      call para_abort
      call exit
end
