!--------------------------
  subroutine parallelinfo
!--------------------------
!
! Write hostname of master node and numbers of processes and threads
!
      use modparallel, only : master, nproc1
!$    use omp_lib
      implicit none
      integer :: nthread, istat, hostnm, llen, len_trim
      character(len=64) :: hostname
!     write(*,'("james parallelinfo.F90 line 12 nthread")')
      nthread=1
!$OMP parallel
!$OMP master
!$    nthread= omp_get_num_threads()
!$OMP end master
!$OMP end parallel
!
!      write(*,'("james like i know what the hell this does.. parallelinfo.F90")')
      if(master) then
        istat= hostnm(hostname)
        llen= len_trim(hostname)
        write(*,'(" Master node is ",a)')hostname(1:llen)
!
        write(*,'(" Number of processes =",i6  )')nproc1
        write(*,'(" Number of threads   =",i6,/)')nthread
      endif
      return
end
 
