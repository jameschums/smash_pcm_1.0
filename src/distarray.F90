!------------------------------------------------------------------------
  subroutine distarray(idis,nmo,nao,nao3,nocca,nvira,noccb,nvirb,nproc)
!------------------------------------------------------------------------
!
! Distribute arrays
!
      implicit none
      integer,intent(in) :: nmo, nao, nao3, nocca, nvira, noccb, nvirb, nproc
      integer,intent(out) :: idis(nproc,14)
      integer :: isize1, isize2, isize3, isize4, i, istart, iend
!
      isize1= (nmo -1)/nproc+1
      isize2= (nao3-1)/nproc+1
      isize3= (nvira-1)/nproc+1
      isize4= (nvirb-1)/nproc+1
      do i=1,nproc
        istart=isize1*(i-1)
        if(istart >= nmo) then
          idis(i,1)= 0
          idis(i,2)= 1
          idis(i,3)= 0
          idis(i,4)= 1
        else
          iend  =isize1*i
          if(iend > nmo) iend=nmo
          idis(i,1)= iend-istart
          idis(i,2)= istart
          idis(i,3)=(iend-istart)*nao
          idis(i,4)= istart*nao
        endif
        istart=isize2*(i-1)
        if(istart >= nao3) then
          idis(i,5)= 0
          idis(i,6)= 1
        else
          iend  =isize2*i
          if(iend > nao3) iend=nao3
          idis(i,5)= iend-istart
          idis(i,6)= istart
        endif
        istart=isize3*(i-1)
        if(istart >= nvira) then
          idis(i, 7)= 0
          idis(i, 8)= 1
          idis(i, 9)= 0
          idis(i,10)= 1
        else
          iend  =isize3*i
          if(iend > nvira) iend=nvira
          idis(i, 7)= iend-istart
          idis(i, 8)= istart
          idis(i, 9)=(iend-istart)*nocca
          idis(i,10)= istart*nocca
        endif
        istart=isize4*(i-1)
        if(istart >= nvirb) then
          idis(i,11)= 0
          idis(i,12)= 1
          idis(i,13)= 0
          idis(i,14)= 1
        else
          iend  =isize4*i
          if(iend > nvirb) iend=nvirb
          idis(i,11)= iend-istart
          idis(i,12)= istart
          idis(i,13)=(iend-istart)*noccb
          idis(i,14)= istart*noccb
        endif
      enddo
!
      return
end

