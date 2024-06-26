!---------------------------
  subroutine tstamp(indext)
!---------------------------
!
! Print CPU and elapsed times
! indext = 0 : at the beginning of the program
! indext = 1 : at the end of a step
! indext = 2 : at the end of the program
!
      use modtclock, only : cpu0, cpu1, iwall0, iwall1
      use modparallel, only : master
      implicit none
      integer :: indext, iwall2, iwrate, iwmax, iday, ihour, imin
      real(8) :: cpu2, wall0, wall1, sec
      character(len=24) :: tdate
!
      write(*,'("james tstamp.F90")')
      
      if(.not.master) return
      if(indext == 0) then
        call cpu_time(cpu0)
        call system_clock(iwall0,iwrate,iwmax)
        call fdate(tdate)
        write(*,'(" The job started at ",a)')tdate
        cpu1= cpu0
        iwall1= iwall0
!
      elseif(indext == 1) then
        call cpu_time(cpu2)
        call system_clock(iwall2,iwrate,iwmax)
        if(iwall2 < iwall1) then
          iwall0= iwall0-iwmax
          iwall1= iwall1-iwmax
        endif
        wall0= dble(iwall2-iwall0)/dble(iwrate)
        wall1= dble(iwall2-iwall1)/dble(iwrate)
        call fdate(tdate)
        write(*,'(1x,"Step CPU :",f10.1,", Total CPU :",f10.1,&
&             " of Master node")') cpu2-cpu1, cpu2-cpu0
        write(*,'(1x,"Step Wall :",f9.1", Total Wall :",f9.1,&
&             " at ",a24,/)') wall1,wall0,tdate
        cpu1 = cpu2
        iwall1= iwall2
!
      elseif(indext == 2) then
        call cpu_time(cpu2)
        call system_clock(iwall2,iwrate,iwmax)
        if(iwall2 < iwall1) then
          iwall0= iwall0-iwmax
          iwall1= iwall1-iwmax
        endif
        wall0= dble(iwall2-iwall0)/dble(iwrate)
        call fdate(tdate)
        iday =(iwall2-iwall0)/iwrate/86400
        ihour= mod((iwall2-iwall0)/iwrate,86400)/3600
        imin = mod((iwall2-iwall0)/iwrate,3600)/60
        sec  = dble(iwall2-iwall0)/dble(iwrate)-dble(86400*iday+3600*ihour+60*imin)
        write(*,'(1x,"Total CPU time :",f11.1," seconds")') cpu2
        write(*,'(1x,"Total Wall time:",f11.1," seconds")') wall0
        write(*,'(17x,"(",i2," days",i3," hours",i3," minutes",f5.1," seconds)")')&
&                   iday, ihour, imin, sec
        write(*,'(" The job finished at ",a)')tdate
      endif
      return
end

