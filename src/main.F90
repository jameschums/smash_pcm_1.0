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
!---------------
  program main
!---------------
!
! This is the main driver of Scalable Molecular Analysis Solver 
! for High performance computing systems (SMASH).
!
      use modparallel, only : master, nproc1, nproc2, myrank1, myrank2, mpi_comm1, mpi_comm2
      use modwarn, only : nwarn
      use modmemory, only : memusedmax
      use modjob, only : runtype, scftype
      use modiofile, only : input, icheck, check, version
      implicit none
      logical :: converged
!
      call setparallel
      version='2.3.0'
!
      if(master) then
        write(*,&
&           '(" *******************************************",/,&
&             "    Scalable Molecular Analysis Solver for",/,&
&             "      High performance computing systems",/,&
&             "            SMASH Version ",a10/,&
&             "           written by K. ISHIMURA",/,&
&             "           edited by J. Robinson Nov fed39",/,&
&             "           Gfortran lapack mpi fedora 39 linux AMD lenovo ",/,&
&             "----------- JJR mp2 bits only SMSH 2.3.0--------",/)') version
      endif
      call tstamp(0)
      call parallelinfo
!
! Read input file and set details
!
      call setdetails(mpi_comm1)
!
! Start calculations
!
      if(scftype == 'RHF') then
        if(runtype == 'ENERGY') then
          call calcrenergy(nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
        elseif(runtype == 'GRADIENT') then
          call calcrgradient(nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
        elseif(runtype == 'OPTIMIZE') then
          call calcrgeometry(converged,nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
        else
          if(master) then
            write(*,'(" Error! This program does not support runtype= ",a16,".")')runtype
            call iabort
          endif
        endif
      elseif(scftype == 'UHF') then
        if(runtype == 'ENERGY') then
          call calcuenergy(nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
        elseif(runtype == 'GRADIENT') then
          call calcugradient(nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
        elseif(runtype == 'OPTIMIZE') then
          call calcugeometry(converged,nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
        endif
      else
        if(master) write(*,'(" Error! SCFtype=",a16," is not supported.")')scftype
        call iabort
      endif
!
! Close input.dat and checkpoint files
!
      if(master) close(unit=input,status='DELETE')
      if(master.and.(check /= '')) close(unit=icheck)
!
      call para_finalize
      call memcheck
      call tstamp(2)
      if(master) then
        write(*,'(" Used memory :",1x,i6," MB")')memusedmax/125000
        if((runtype =='OPTIMIZE').and.(.not.converged))then
          write(*,'(/," ============================================================")')
          write(*,'("  Geometry optimization did not finish with",i3," warning(s)!")')nwarn
          write(*,'(" ============================================================")')
        else
          write(*,'("smash2.3 Lenovo IdeaPad 5 Pro 16ACH AMD Ryzen 9 5900HX with Radeon Graphics x16")')
	      write(*,'("smash2.3 NVIDIA Corporation GA107M [GeForce RTX 3050 Mobile] (rev a1) ")')
	      write(*,'("smash2.3 Your calculation finished with",i3," warning(s) nvfotran AMD Lenovo Fedora39.")')nwarn
        endif
      endif
end program main



