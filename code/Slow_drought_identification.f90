subroutine slowdrought(n,a,b,drt)
!======================================================================
! References:
! 1. Yuan, X., L. Wang, P. Wu, P. Ji, J. Sheffield, and M. Zhang, 2019: 
! Anthropogenic shift towards higher risk of flash drought over China. 
! Nature Communications, 10, 4661, 
! https://doi.org/10.1038/s41467-019-12692-7
! 2. Yuan, X., Y. Wang, P. Ji, P. Wu, J. Sheffield, J. Otkin, 2023: A 
! global transition to flash droughts under climate change. Science, 
! https://doi.org/10.1126/science.abn6301
!======================================================================
! n, number of pentads for a given year
! a(n), pentad soil moisture percentile [0-100] 
! b(1), number of slow droughts; 
! b(2), mean duration of slow droughts
! b(3), mean severity of slow droughts 
! b(4), mean speed of slow droughts
! drt, 1-drought onset, 2-drought recovery, 0-no drought
!======================================================================
implicit none
real,parameter::thresh=40.   ! drought start threshold [percentile]
real,parameter::thresh1=20.  ! drought end threshold [percentile]
real,parameter::speed=5.     ! drought speed threshold [percentile]
integer,parameter::td=4      ! drought duration threshold [pentads]
integer n,i,j,k
real a(n),b(4),tmin,drt(n)
integer si,cnt,flag1
integer flag   ! flag=1, begin drought identification
integer flag2  ! flag2=1, onset period; flag2=0, recovery period
real odur(n)   ! total duration of each flash drought
real osev(n)   ! total severity of each flash drought
real ospd(n)   ! onset speed of each flash drought
integer cnt1,cnt2
real odur1(n)  ! total duration of each drought
real osev1(n)  ! total severity of each drought
real ospd1(n)  ! onset speed of each drought
real drt1(n)   


do i=1,n
   if(a(i)<0.)then
      print*,'error, stop'
      print*,a
      stop
   endif
end do

cnt=0; flag=0; si=0
odur1=0; osev1=0.; ospd1=0.
tmin=99.
flag1=0; flag2=0
drt1=0
do j=2,n
   if(a(j)<thresh)then
      if(flag==0)then                ! drought start
         flag1=1
         if(j>1)then
            if(a(j-1)<thresh) flag1=0
         endif
         if(flag1==1)then
            flag=1
            flag2=1                  ! onset period
            cnt=cnt+1
            si=j                     ! si is drought start time
            osev1(cnt)=thresh-a(j)
            odur1(cnt)=1
            ospd1(cnt)=thresh-a(j)
            drt1(j)=1
            if(a(j)<tmin) tmin=a(j)
         endif
      else
         if(flag2==1)then            ! onset period
            if(tmin<=thresh1)then        ! enter thresh1 (e.g., <20%)
               if(a(j)>=a(j-1))then      ! stop onset period
                  flag2=0
                  if(a(j)<=thresh1)then  ! enter recovery period
                     osev1(cnt)=osev1(cnt)+(thresh-a(j))
                     odur1(cnt)=odur1(cnt)+1
                     drt1(j)=2
                     if(a(j)<tmin) tmin=a(j)
                  else                   ! no recovery period
                     flag=0
                     tmin=99.
                     if(odur1(cnt)<td)then
                        osev1(cnt)=0.
                        odur1(cnt)=0
                        ospd1(cnt)=0.
                        drt1(si:j-1)=0
                        cnt=cnt-1
                     endif
                  endif
               else                      ! continue onset period
                  osev1(cnt)=osev1(cnt)+(thresh-a(j))
                  odur1(cnt)=odur1(cnt)+1
                  ospd1(cnt)=(thresh-a(j))/real(j-si+1)
                  drt1(j)=1
                  if(a(j)<tmin) tmin=a(j)
               endif
            else                        ! before entering thresh1, continue onset period
                osev1(cnt)=osev1(cnt)+(thresh-a(j))
                odur1(cnt)=odur1(cnt)+1
                ospd1(cnt)=(thresh-a(j))/real(j-si+1)
                drt1(j)=1
                if(a(j)<tmin) tmin=a(j)
            endif
         else  ! flag2=0, recovery period
            if(a(j)<=thresh1)then
               osev1(cnt)=osev1(cnt)+(thresh-a(j))
               odur1(cnt)=odur1(cnt)+1
               drt1(j)=2
            else                  ! drought dismiss
               if(odur1(cnt)<td)then
                  osev1(cnt)=0.
                  odur1(cnt)=0
                  ospd1(cnt)=0.
                  drt1(si:j-1)=0
                  cnt=cnt-1
               endif
               flag=0
               tmin=99.
            endif
         endif
      endif
      if(j==n .and. flag==1)then
         if(odur1(cnt)<td .or. tmin>thresh1)then
            osev1(cnt)=0.
            odur1(cnt)=0
            ospd1(cnt)=0.
            drt1(si:j)=0
            cnt=cnt-1
         endif
      endif
   else                           ! a(j)>=thresh
      if(flag==1)then             ! drought dismiss
         flag=0
         if(odur1(cnt)<td .or. tmin>thresh1)then
            drt1(si:j)=0
            osev1(cnt)=0.
            odur1(cnt)=0
            ospd1(cnt)=0.
            cnt=cnt-1
         endif
         tmin=99.
      endif
   endif
end do
cnt1=cnt

cnt=0; flag=0; si=0
odur=0; osev=0.; ospd=0.
tmin=99.
flag1=0; flag2=0
drt=0

do j=2,n
   if(a(j)<thresh)then
      if(flag==0)then                ! drought start
         flag1=1
         if(j>1)then
            if(a(j-1)<thresh) flag1=0
         endif
         if(flag1==1)then
            flag=1
            flag2=1                  ! onset period
            cnt=cnt+1
            si=j                     ! si is drought start time
            osev(cnt)=thresh-a(j)
            odur(cnt)=1
            ospd(cnt)=thresh-a(j)
            drt(j)=1
            if(a(j)<tmin) tmin=a(j)
         endif
      else
         if(flag2==1)then            ! onset period
            if(thresh-a(j)>=(j-si+1)*speed)then
               if(tmin<=thresh1)then        ! enter thresh1 (e.g., <20%)
                  if(a(j)>=a(j-1))then      ! stop onset period
                     flag2=0
                     if(a(j)<=thresh1)then   ! enter recovery period
                        osev(cnt)=osev(cnt)+(thresh-a(j))
                        odur(cnt)=odur(cnt)+1
                        drt(j)=2
                        if(a(j)<tmin) tmin=a(j)
                     else                   ! no recovery period
                        flag=0
                        tmin=99.
                        if(odur(cnt)<td)then
                           osev(cnt)=0.
                           odur(cnt)=0
                           ospd(cnt)=0.
                           drt(si:j-1)=0
                           cnt=cnt-1
                        endif
                     endif
                  else                      ! continue onset period
                     osev(cnt)=osev(cnt)+(thresh-a(j))
                     odur(cnt)=odur(cnt)+1
                     ospd(cnt)=(thresh-a(j))/real(j-si+1)
                     drt(j)=1
                     if(a(j)<tmin) tmin=a(j)
                  endif
               else                         ! before entering thresh1, continue onset period
                   osev(cnt)=osev(cnt)+(thresh-a(j))
                   odur(cnt)=odur(cnt)+1
                   ospd(cnt)=(thresh-a(j))/real(j-si+1)
                   drt(j)=1
                   if(a(j)<tmin) tmin=a(j)
               endif
            else               ! speed does not meet      
               flag2=0         ! so enter recovery period
               if(tmin<=thresh1)then
                  if(a(j)<=thresh1)then
                     osev(cnt) = osev(cnt)+(thresh-a(j))
                     odur(cnt) = odur(cnt)+1
                     drt(j) = 2
                  else            !drought dismiss
                     flag=0
                     tmin=99.
                     if(odur(cnt)<td)then
                        osev(cnt)=0.
                        odur(cnt)=0
                        ospd(cnt)=0.
                        drt(si:j-1)=0
                        cnt=cnt-1
                     endif
                  endif
               else            ! tmin>thresh1, does not meet drought criterion 
                  drt(si:j-1)=0
                  osev(cnt)=0.
                  odur(cnt)=0
                  ospd(cnt)=0.
                  cnt=cnt-1
                  flag=0
                  tmin=99.
               end if
            endif
         else  ! flag2=0, recovery period
            if(a(j)<=thresh1)then
               osev(cnt)=osev(cnt)+(thresh-a(j))
               odur(cnt)=odur(cnt)+1
               drt(j)=2
            else                  ! drought dismiss
               if(odur(cnt)<td)then
                  osev(cnt)=0.
                  odur(cnt)=0
                  ospd(cnt)=0.
                  drt(si:j-1)=0
                  cnt=cnt-1
               endif
               flag=0
               tmin=99.
            endif
         endif
      endif
      if(j==n .and. flag==1)then
         if(odur(cnt)<td .or. tmin>thresh1)then
            osev(cnt)=0.
            odur(cnt)=0
            ospd(cnt)=0.
            drt(si:j)=0
            cnt=cnt-1
         endif
      endif
   else                           ! a(j)>=thresh
      if(flag==1)then             ! drought dismiss
         flag=0
         if(odur(cnt)<td .or. tmin>thresh1)then
            drt(si:j)=0
            osev(cnt)=0.
            odur(cnt)=0
            ospd(cnt)=0.
            cnt=cnt-1
         endif
         tmin=99.
      endif
   endif
end do

where(drt>0) drt1=0
drt=drt1

ospd1=0.; cnt2=0; si=0
do i=2,n 
   if(drt(i)==1. .and. drt(i-1)==0.)then
      si=i
      cnt2=cnt2+1
   endif
   if(drt(i-1)==1. .and. abs(drt(i)-1)>0.5) ospd1(cnt2)=(thresh-a(i-1))/(i-si)
   if(i==n .and. drt(i)==1.) ospd1(cnt2)=(thresh-a(i))/(i-si+1)
end do

if(cnt1>0)then
   b(1)=real(cnt1-cnt)
   if(b(1)>0)then
      b(2)=real(sum(odur1(1:cnt1))-sum(odur(1:cnt)))/b(1)
      b(3)=(sum(osev1(1:cnt1))-sum(osev(1:cnt)))/b(1)
      b(4)=sum(ospd1(1:cnt2))/b(1)
    else
      b=0 
    endif
else
   b=0
endif
end

