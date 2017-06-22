program TraFFT


!-------------------------------------------------------------
!
!   TraFFT
!
!
!      
!
!
!                                 
!                                    2010.08. FUJI Nobuaki
!                                    2010.10. FUJI Nobuaki
!------------------------------------------------------------- 

  implicit none
  integer :: iPS = 3 ! 1:SH, 2:PSV, 3:ALL
  character(120) :: outputDir,psvmodel,modelname,stationsinf,parentDir
  integer :: itranslat
  character(120) :: coutfile
  real(kind(0d0)) :: tlen,r0min,r0max,r0delta,r0lat,r0lon

  real(kind(0d0)), allocatable :: stla(:),stlo(:)
  real(kind(0d0)), allocatable :: r_(:),r0(:),theta(:),phi(:)
  integer, allocatable :: updown(:) 
  integer :: r_n,r0_n,ir0,imt,theta_n,nsta
  integer :: i,j,jj,imin,imax,np

  !for FFT, imin should be 0
  real(kind(0d0)) :: omegai
  real(kind(0d0)), parameter :: pi=3.1415926535897932d0
  integer,parameter :: maxradsamples = 5
  ! FOR FOURIER TRANSFORM & BUTTERWORTH FILTERING
  
 
  real(kind(0d0)) :: mt(6) 

  integer :: iwindowStart, iwindowEnd
  integer, parameter:: ibwfilt = 0 ! if =0, we don't perform filtering
  real(kind(0d0)) :: samplingHz = 1.d1
  integer :: lsmooth,np0,np1
  real(kind(0d0)),parameter :: start = 0.d0 ! Time window
  real(kind(0d0)),parameter :: end = 1.5d3  ! in seconds

  ! RSGT & TSGT

  complex(kind(0e0)),allocatable :: stresssngl(:,:,:,:), displacementsngl(:,:,:,:),tmpsngl(:,:,:)
  complex(kind(0d0)), allocatable :: gt(:,:)
  real(kind(0e0)), allocatable :: ygt(:,:)
  real(kind(0d0)), allocatable :: tmpygt(:)
!-----------------------------------------------------------------------


  
  call pinputTra(outputDir,psvmodel,modelname,stationsinf,tlen,imin,imax,r0min,r0max,r0delta,r0lat,r0lon,itranslat,mt)
  np = imax
  
  parentdir = trim(outputDir)//"/ascii/"
  omegai = - dlog(1.d-2) / tlen



   
  open (1,file=stationsinf,status='old',action='read',position='rewind')
  !if(itranslat.eq.1) call translat (r0lat,r0lat)
  read(1,*)nsta
  r_n = nsta
  theta_n = nsta
  allocate(r_(1:r_n))
  allocate(theta(1:theta_n))
  allocate(phi(1:theta_n))
  allocate(stla(1:nsta))
  allocate(stlo(1:nsta))
  allocate(updown(1:nsta))

  do i = 1,nsta
     read(1,*) r_(i),stla(i),stlo(i),updown(i)
     r_(i) = 6371.d0 -r_(i)
     if(itranslat.eq.1) call translat(stla(i),stla(i))
     call calthetaphi(r0lat,r0lon,stla(i),stlo(i),theta(i),phi(i))
  enddo
  close(1)

  r0_n =  int((r0max-r0min)/r0delta)+1
  allocate(r0(1:r0_n))
  do i = 1, r0_n
     r0(i) = r0min + dble(i-1)*r0delta
  enddo
 


  ! lsmoothfinder for FFT

  np0 = np
  call lsmoothfinder(tlen, np0, samplingHz, lsmooth)

  i=1
  do while (i<lsmooth)
     i = i*2
  enddo
  lsmooth = i
  i = 0

  np1 = 1
  do while (np1<np0)
     np1 = np1*2
  enddo

  np1 = np1*lsmooth

  samplingHz = dble(2*np1)/tlen  
  iWindowStart = int(start*samplingHz)
  iWindowEnd   = int(end*samplingHz)



  
  allocate(gt(1:9,0:2*np1-1))
  allocate(ygt(1:9,0:2*np1-1))
  allocate(tmpygt(0:2*np1-1))
  allocate(stresssngl(1:6,1:6,1:nsta,imin:imax))
  allocate(displacementsngl(1:3,1:6,1:nsta,imin:imax))
  allocate(tmpsngl(1:6,1:6,1:nsta))

  stresssngl=cmplx(0.e0)
  displacementsngl=cmplx(0.e0)
  ir0 = r0_n

  do i =imin,imax
     if(iPS.ne.1) then ! PSV calculation
        write(coutfile, '(I5,".Stress_PSV")'), int(r0(ir0)*10.d0)
        do j = 1,7
           if (coutfile(j:j).eq.' ')coutfile(j:j) = '0'
        enddo
        
        coutfile = trim(modelname)//"."//coutfile
        coutfile = trim(outputDir)//"/Stress/"//coutfile
        
        open(1,file=coutfile,status='unknown',form='unformatted', &
             access = 'direct', recl=2*6*6*kind(0e0)*nsta)
        read(1,rec=i+1) tmpsngl(1:6,1:6,1:nsta)
        stresssngl(1:6,1:6,1:nsta,i)=stresssngl(1:6,1:6,1:nsta,i)+tmpsngl(1:6,1:6,1:nsta)
        !print *, i,stresssngl(1:6,1:1,i)
        close(1)         
        
        write(coutfile, '(I5,".Displacement_PSV")'), int(r0(ir0)*10.d0)
        do j = 1,7
           if (coutfile(j:j).eq.' ')coutfile(j:j) = '0'
        enddo
        
        coutfile = trim(modelname)//"."//coutfile
        coutfile = trim(outputDir)//"/Displacement/"//coutfile
        
        open(1,file=coutfile,status='unknown',form='unformatted', &
             access = 'direct', recl=2*3*6*kind(0e0)*nsta)
        read (1,rec=i+1) tmpsngl(1:3,1:6,1:nsta)
        displacementsngl(1:3,1:6,1:nsta,i)=displacementsngl(1:3,1:6,1:nsta,i)+tmpsngl(1:3,1:6,1:nsta)

        !print *, i,displacementsngl(1:3,1:6,1:1,i)
        close(1)     
     endif
       if(iPS.ne.2) then ! SH calculation
          write(coutfile, '(I5,".Stress_SH")'), int(r0(ir0)*10.d0)
        do j = 1,7
           if (coutfile(j:j).eq.' ')coutfile(j:j) = '0'
        enddo
        
        coutfile = trim(modelname)//"."//coutfile
        coutfile = trim(outputDir)//"/Stress/"//coutfile
        
        open(1,file=coutfile,status='unknown',form='unformatted', &
             access = 'direct', recl=2*6*6*kind(0e0)*nsta)
        read(1,rec=i+1) tmpsngl(1:6,1:6,1:nsta)
        stresssngl(1:6,1:6,1:nsta,i)=stresssngl(1:6,1:6,1:nsta,i)+tmpsngl(1:6,1:6,1:nsta)
        !print *, i,stresssngl(1:6,1:1,i)
        close(1)         
        
        write(coutfile, '(I5,".Displacement_SH")'), int(r0(ir0)*10.d0)
        do j = 1,7
           if (coutfile(j:j).eq.' ')coutfile(j:j) = '0'
        enddo
        
        coutfile = trim(modelname)//"."//coutfile
        coutfile = trim(outputDir)//"/Displacement/"//coutfile
        
        open(1,file=coutfile,status='unknown',form='unformatted', &
             access = 'direct', recl=2*3*6*kind(0e0)*nsta)
        read (1,rec=i+1) tmpsngl(1:3,1:6,1:nsta)
        displacementsngl(1:3,1:6,1:nsta,i)=displacementsngl(1:3,1:6,1:nsta,i)+tmpsngl(1:3,1:6,1:nsta)

        !print *, i,displacementsngl(1:3,1:6,1:1,i)
        close(1)  
     endif
  enddo
  
  
  do i = 1,nsta
     gt = cmplx(0.d0)
    
        do imt=1,6
           gt(1:6,imin:imax)=gt(1:6,imin:imax)+stresssngl(1:6,imt,i,imin:imax)*mt(imt)
           gt(7:9,imin:imax)=gt(7:9,imin:imax)+displacementsngl(1:3,imt,i,imin:imax)*mt(imt)
           
        enddo
        
        
        
        ygt = 0.e0
        call tensorFFT_real(9,imin,np1,gt,ygt,omegai,tlen)
        !do jj = 1,9
        !      call bwfilt(dble(ygt(jj,iWindowStart:iWindowEnd)),tmpygt(iWindowStart:iWindowEnd),1.d0/samplingHz,iWindowEnd-iWindowStart+1,1,4,1.d-2,5.d-1)
        !      ygt(3,:) = real(tmpygt(:))
        ! FOR TEST
        do jj = 1,9
           write(coutfile, '("green",I5.5,I3.3)') ,i,jj
           
           do j = 1,9
              if (coutfile(j:j).eq.' ')coutfile(j:j) = '0'
           enddo
           coutfile = trim(parentDir)//"/"//trim(coutfile)
           !print *, coutfile
           open(1,file=coutfile,status='unknown', form='formatted')
           do j = iWindowStart,iWindowEnd
              write(1,*) dble(j)/samplingHz,ygt(jj,j)
           enddo
           close(1)
        enddo
        ! FOR TEST
   
  enddo
 


end program TraFFT
