!
! This routine was developed to be self-contained with pgplot.
! It is based on $cmfdisp/spec_plt/subs/smear_v2.f and
!                $cmfdisp/spec_plt/subs/uvabs_v2.fi .
! that were originally developed by Jim Herald. Unnecessary 
! routines were deleted, and FFT option was removed.
! Calls to tune were also deleted.
!
! Created 22-MAr-2023
!
	MODULE MOD_SMEAR_PG
	  PUBLIC SMEAR_V3_PG
	  PRIVATE nonfftconvolve,fillresponse, fill_rot_response, linearize_v2,
	1         nu_to_lambda, lambda_to_lnlambda, map, map2, gauss,
	1         log2, larger_index , smaller_index 
        CONTAINS
!----------------------------------------------------------------------------
!
! Subroutine to convolve data down to a given resolution.
! Three main options are available?
!
! (1) Smooth to fixed resolution dLAM (=inst_res) over the indicated passband. 
!        Smoothing is done in wavelength plane.
! (2) Smooth so that Lambda/dLam (=Resolution) is fixed over the indicated
!       pasband. Smoothing is done in the log-wavelength plane, since smoothing
!       function is then independent of wavelngth.
! (3) Perform a crude rotational broadening with 2 free parameters --- VSINI
!       and a limb darkening parameter EPSILON. This is similar to (2) in that
!       Lambda/dLam is fixed over the indicated pasband. Smoothing is done in 
!       the log-wavelngth plane, since smoothing function is then independent 
!       of wavelngth.
!
! NB: dLAM is interpreted as the FWHM of the instrumental gaussian profile in
!       options (1) and (2).
!
	subroutine smear_v3_pg(origfreq,origflux,Norig,
	1                wave_max,wave_min,
	1                inst_res,resolution,vsini,epsilon,
	1                min_res_kms,num_res)
!
	implicit none
	include 'constants.inc'
!
! Altered 21-Mar-2023 : Extended range now activated when "resolution" or "vsini" is set.
! Altered 18-Apr-2008 : tempwave & tempflux now allocated with max(Nmod,Ntmp)
! 27-Jan-2000 : Cleaned; Converted to version 2.
!               Rotation option installed.
!
	integer Norig
	real*8 origwave(Norig),origflux(Norig)
	real*8 wave_max,wave_min
	real*8 inst_res     	! the instrumental resolution (dlambda)
	real*8 Resolution   	! resolving power (R = lambda/dlambda)
	real*8 min_res_kms  	! minimum resolution to work with
	real*8 num_res      	! num resolution elements to consider beyond
                            	! wavelength range.
	real*8 vsini		! Projected rotational velocity in km/s
	real*8 epsilon		! Limb darkenin parameter: 
	                        !         I/I(o)=1 - epsilon + epsilon mu
!
	real*8 origfreq(Norig)
	real*8, allocatable :: modwave(:)
	real*8, allocatable :: modflux(:)
	real*8, allocatable :: tempwave(:)
	real*8, allocatable :: tempflux(:)
!
	integer Nmod
	integer Ntmp
	real*8 cwave_max,cwave_min
	real*8 model_res    ! the model resolution
	real*8 min_dlam     ! the model resolution
	real*8 kernal_sig   ! the corresponding sigma (fwhm = 2.354sig)
	real*8 extend       ! wavelength extension to consider
                            ! in convolution
	real*8 t1
!
	integer alter_min ! min and max indices of original array to be
	integer alter_max ! changed.
	integer con_min   ! the min and max indices of original
        integer con_max   ! to be *considered*
	integer i_min,i_max
	integer ios
!
        call nu_to_lambda(origfreq,origwave,Norig,forward)
!
	if(vsini .ne. 0.0d0)then
	  extend=2*(wave_max+wave_min)*vsini/3.0e+05
	else if(resolution .ne. 0.0d0)then
	  extend=num_res * (wave_max+wave_min)/resolution
	else
	  extend = num_res * inst_res
	end if
	cwave_max = wave_max + extend
	cwave_min = wave_min - extend
	
	con_min = larger_index(origwave,Norig,cwave_min)
	con_max = smaller_index(origwave,Norig,cwave_max)
	alter_min = larger_index(origwave,Norig,wave_min)
	alter_max = smaller_index(origwave,Norig,wave_max)
	Nmod = con_max - con_min + 1
!
	write(*,*)' '
	write(*,*)'Considering wavelength range of ',cwave_min,' to ',cwave_max
	write(*,*)'For convolution, only ',wave_min,' to ',wave_max,' of'
	write(*,*)'       original data will be altered'
!
	if(allocated(modwave))deallocate(modwave)
	if(allocated(modflux))deallocate(modflux)
	allocate (modwave(Nmod),stat=ios)
	if(ios .eq. 0)allocate (modflux(Nmod),stat=ios)
	if(ios .ne. 0)then
	  write(6,*)'Error allocating memory for modwave in smear'
	  return
	end if
!
	modwave(1:Nmod) = origwave(con_min:con_max) !preserve original wave and
	modflux(1:Nmod) = origflux(con_min:con_max) !flux arrays for later
!
! We interpret RESOLUTION to be LAMBDA/dLAM where dLAM in the Full Width at
! Half Maximum (FWHM) of the instrumental profile.
!
	if (resolution .ne. zero .or. vsini .ne. 0) then
         call lambda_to_lnlambda(modwave,Nmod,forward)
	end if
	if(vsini .ne. 0)then
	 min_dlam=min_res_kms/c_kms
	else if(resolution .ne. 0)then
	 inst_res = dlog((two*resolution+one)/(two*resolution-one))
	 min_dlam=min_res_kms/c_kms
	else
	 min_dlam=wave_min*min_res_kms/c_kms
	endif
!
	model_res=minval(modwave(2:Nmod)-modwave(1:Nmod-1))
	if (model_res .lt. min_dlam)model_res = min_dlam
	write(*,*)'Working model_res = ',model_res
!
	if(vsini .ne. 0)then
	  model_res=min(0.2d0*vsini/c_kms,model_res)
	else if(resolution .ne. 0)then
	  model_res=min(0.2d0*inst_res,model_res)
	else
	  model_res=min(0.2d0*inst_res,model_res)
	end if
!
	ntmp=(modwave(nmod)-modwave(1))/model_res+1
	if(allocated(tempwave))deallocate(tempwave)
	if(allocated(tempflux))deallocate(tempflux)
	allocate (tempwave(Max(Ntmp,Nmod)),stat=ios)
	if(ios .eq. 0)allocate (tempflux(Max(Ntmp,Nmod)),stat=ios)
	if(ios .ne. 0)then
	  write(6,*)'Error allocating memory for tempwave in smear'
	  return
	end if
	tempwave(1:Nmod)=modwave(1:Nmod)
	tempflux(1:Nmod)=modflux(1:Nmod)
!
	call linearize_v2(tempwave,tempflux,Nmod,Ntmp,model_res)
!
! inst_res refers to the FWHM of the instrumental profile. This is related
! to sigma of the Gaussian profile ( EXP[-0.5 {(x-mu)/sig}^2]) as indicated
! below.
!
	kernal_sig = inst_res/2.35482
	if(resolution .eq. 0)then
	  t1=0.5D0*(wave_min+wave_max)
	  write(6,'(A,F7.2,A)')' Sigma of smoothing Gausian at middle of lambda range is ',kernal_sig*c_kms/t1,' km/s'
	else
	  write(6,'(A,F7.2,A)')' Sigma of smoothing Gausian is ',kernal_sig*c_kms,' km/s'
	end if
	write(6,*)' '
	call nonfftconvolve(tempwave,tempflux,Ntmp,kernal_sig,vsini,epsilon)
!                               
!	if (resolution .ne. zero) then
!	 call lambda_to_lnlambda(tempwave,ntmp,backward)
!	endif
!                     
	call map(tempwave,tempflux,Ntmp,modwave,modflux,Nmod)
!
	i_min = alter_min-con_min+1
	i_max = alter_max-con_min+1
	origflux(alter_min:alter_max) = modflux(i_min:i_max)
!
	return 
	end 

C-----------------------------------------------------------------------------
C subroutine to convolve data with gaussian of the given sigma (non-FFT).
C It assumes that data is evenly spaced in wavelength.
C

	subroutine nonfftconvolve(wave,flux,Nmod,sigma,vsini,epsilon)
	implicit none
	include 'constants.inc'

	integer Nmod
	real*8 wave(Nmod)
	real*8 flux(Nmod)
	real*8 sigma
	real*8 vsini
	real*8 epsilon
!
	real*8 answer(Nmod)
	real*8 response(Nmod)
	integer Nresponse,nron2
	integer i       ! index for convolution array
	integer j       ! index for response array
	integer k       ! index for model array
	integer l       ! temporay index for model array in concolution
!
	if(vsini .eq. 0)then
	  call fillresponse(wave,response,Nmod,Nresponse,sigma,izero)
	else
	  call fill_rot_response(wave,response,Nmod,Nresponse,vsini,epsilon,izero)
	end if
	nron2=Nresponse/2
!
	answer(1:Nmod)=0.0d0
!
! Three cases: At the boudaries we assume symmetry about the boundary.
!
! 1  +-----+
!   /\     |
!
	do j=-nron2,nron2
	 do i=1,nron2
	    k=j+1+nron2
	    l=i+j      
	    if(l .le. 0)l=2-l
	    answer(i) = answer(i) + flux(l)*response(k)
	  end do
	end do
!
! 2  +-----+
!    |  /\ |
!
	do j=-nron2,nron2
	 do i=1+nron2,nmod-nron2
	    k=j+1+nron2
	    answer(i) = answer(i) + flux(i+j)*response(k)
	  end do
	end do
!
! 2  +-----+
!    |     /\ |
!
	do j=-nron2,nron2
	  do i=nmod-nron2+1,nmod
	    k=j+1+nron2
	    l=i+j      
	    if(l .gt. nmod)l=2*Nmod-l
	    answer(i) = answer(i) + flux(l)*response(k)
	  end do
	end do
!
	flux(1:Nmod)=Answer(1:Nmod)
!	  
	return 
	end

C 
C
C-----------------------------------------------------------------------------
C function fills response array with a gaussian of the given sigma
C to NUM_SIGMA sigmas.  If wrap = 1, then store array in wrap-around order.
C
	subroutine fillresponse(wave,response,Nmod,Nresponse,sigma,wrap)

	implicit none
	include 'constants.inc'

	integer Nmod           ! size of data array
	real*8 wave(Nmod)        ! wavelength array
	real*8 response(Nmod)    ! response array
	real*8 sigma             ! sigma of gaussian 
	integer Nresponse      ! size of response array
	integer wrap           ! fill in wrap-around order?
!
	integer, parameter :: num_sigmas=5
	integer i
	real*8 lambda,dlam,cutoff,mu
	real*8 response_area

	dlam = (wave(Nmod)-wave(1))/(Nmod-1)
	response_area = zero

	cutoff = wave(1) + NUM_SIGMAS*sigma   ! fill response array out to
	                                      ! this wavelength
	i = 1
	do while (wave(i) .le. cutoff)        ! find out what index the 
	   i = i+1                            ! cutoff wavelength 
	end do                                ! corresponds to

	Nresponse = 2*(i-1)+1                 ! the size of the response array

	if (Nresponse .gt. Nmod) then
	   write(*,*)'error: response array larger than data array'
	   return
	end if

	if (wrap .eq. 1) then
	mu = zero
	lambda=zero
	response_area = gauss(lambda,mu,sigma)
	response(1) = gauss(lambda,mu,sigma)   ! fill the response array with
	lambda = dlam                          ! a gaussian of the given sigma
	do i=2, (Nresponse-1)/2+1              ! in wrap-around order, ie. if
	   response(i)=gauss(lambda,mu,sigma)  ! resp = ___/\___ it would
	   response(Nresponse+2-i)=response(i) ! be stored as: \______/
	   lambda = lambda + dlam              ! as per. numerical recipes.
	   response_area = response_area + two * response(i)
	end do
	endif

	if (wrap .eq. 0) then
	   mu = (Nresponse-1)/2 * dlam         ! fill in response array with
	   lambda = zero                       ! a gaussian centered about
	   do i=1,Nresponse                    ! mu characterized by sigma,
             response(i)=gauss(lambda,mu,sigma)! resp = __/\__
	     lambda = lambda + dlam
	     response_area = response_area + response(i)
	   end do
	endif
  
C       now we must multiply the discrete elements of the response function
C       so that their sum adds to one (Normalize the response function)

	do i=1,Nresponse
	   response(i) = response(i)/response_area
	end do

	return
	end
C
C-----------------------------------------------------------------------------
C function fills response array with a gaussian of the given sigma
C to NUM_SIGMA sigmas.  If wrap = 1, then store array in wrap-around order.
C
	subroutine fill_rot_response(wave,response,Nmod,Nresponse,
	1                              vsini,epsilon,wrap)

	implicit none
	include 'constants.inc'
!
	integer Nmod           ! size of data array
	real*8 wave(Nmod)     	! wavelength array
	real*8 response(Nmod) 	! response array
	real*8 vsini
	real*8 epsilon
	integer Nresponse      ! size of response array
	integer wrap           ! fill in wrap-around order?

	integer i
	real*8 lambda,dlam,cutoff,mu
	real*8 response_area
	real*8 t1,dlam_rot,a1,a2
!
	dlam = (wave(Nmod)-wave(1))/(Nmod-1)
	response_area = zero
!
	cutoff = wave(1) + vsini/c_kms        ! fill response array out to
	                                      ! this wavelength
	i = 1       
	do while (wave(i) .lt. cutoff)        ! find out what index the 
	   i = i+1                            ! cutoff wavelength 
	end do                                ! corresponds to
!                 
	Nresponse = 2*(i-1)+1                 ! the size of the response array
	if (Nresponse .gt. Nmod) then
	   write(*,*)'error: response array larger than data array'
	end if
!
	if(epsilon .eq. 1)then
	  a1=0
	  a2=1.0d0
	else
	  a1=2.0D0/jimpi
	  a2=0.5D0*epsilon/(1.0D0-epsilon)
	end if
	dlam_rot=vsini/c_kms
!
	if (wrap .eq. 1) then
	  mu = zero
	  lambda=zero
	  t1=1-(lambda/dlam_rot)**2
	  response_area = a1*sqrt(t1)+a2*t1
	  response(1) = response_area             ! fill the response array with
	  lambda = dlam                           ! a gaussian of the given sigma
	  do i=2, (Nresponse-1)/2+1               ! in wrap-around order, ie. if
	    t1=1-(lambda/dlam_rot)**2
	    if(t1 .le. 0)then
	      response(i)=0.0d0
	    else
	      response(i)=a1*sqrt(t1)+a2*t1         ! resp = ___/\___ it would
	    end if
	    response(Nresponse+2-i)=response(i)   ! be stored as: \______/
	    lambda = lambda + dlam                ! as per. numerical recipes.
	    response_area = response_area + two * response(i)
	  end do
	endif
!
	if (wrap .eq. 0) then
	  mu = (Nresponse-1)/2 * dlam         ! fill in response array with
	  lambda = 0.0D0                      ! a gaussian centered about
	  do i=1,Nresponse                    ! mu characterized by sigma,
	    t1=1.0d0-((lambda-mu)/dlam_rot)**2
	    IF(t1 .le. 0.0D0)then
	      response(i)=0.0d0
	    else
              response(i)=a1*sqrt(t1)+a2*t1       ! resp = __/\__
	    end if
	    lambda = lambda + dlam
	    response_area = response_area + response(i)
	  end do
	endif
!
! Now we must multiply the discrete elements of the response function
! so that their sum adds to one (Normalize the response function)
!
	do i=1,Nresponse
	   response(i) = response(i)/response_area
	end do

	return
	end

C-----------------------------------------------------------------------------
C
C Function returns a gauss shaped profile centered on zero and described
C by sigma
C

	function gauss(x,mu,sigma)

	implicit none
	include 'constants.inc'

	real*8 sigma,x,mu
	real*8 gauss

	gauss = exp(-onehalf*(((x-mu)/sigma)**two))/(sigma*sqrt(two*jimPI))
	return
	end
C---------------------------------------------------------------------------
C
C function returns the log-base-2 of number, rounded down.
C
	function log2(n)

	implicit none
	integer n,log2,junk

	log2 = 0
	junk = n/2
	do while (junk .ne. 0)
	   junk = junk/2
	   log2 = log2+1
	end do
	
	return
	end

C-----------------------------------------------------------------------------
C
C subroutine to make points linearly spaced in wavelength to a given
C delta lambda via linear interpolation
C

	subroutine linearize_v2(wave,flux,n,nmax,dlam)
	implicit none
!
	integer n
	integer nmax
	real*8 wave(NMAX),flux(NMAX)
	real*8 dlam 				! even spacing in wavelength
!
	real*8 tempwave(N),tempflux(N)
	real*8 dlamover2
	integer nnew
	integer i
!
	dlamover2 = 0.5d0*dlam
	nnew = (wave(n)-wave(1))/dlam  ! this will shorten wavelength 
                                          ! array by a fraction of dlam
!
	if (nnew .gt. NMAX) then
	   write(*,*)'Error, NMAX = ',NMAX
	   write(*,*)'too small for this resolution'
	   write(*,*)'Try setting NMAX = ',nnew
	   stop 
	endif
!
	tempwave(1:n)=wave(1:n)
	tempflux(1:n)=flux(1:n)
!	
	wave(1) = tempwave(1) + dlamover2
	do i=2,nnew
	  wave(i) = wave(i-1) + dlam
	end do
!
	call map2(tempwave,tempflux,n,wave,flux,nnew,dlam)
	nmax = nnew
!
	return
	end

C-----------------------------------------------------------------------------
C subroutine to build wavelength array from frequency array (if forward=true)
C or build frequency array from wavelength array (if forward=false)
C 
C
	subroutine nu_to_lambda(freq,wave,n,nutolambda)

	implicit none
	include 'constants.inc'

	integer*4 n
	real*8 freq(n)
	real*8 wave(n)
	logical nutolambda ! true=nu->lambda, false=lambda->nu

	integer*4 i

	if (nutolambda) then
	 do i=1,n
	   if (freq(i) .eq. zero) then
	      write(*,*)'zero frequency found at i=',i
	   endif
	   wave(i) = (c_kms*10D-3)/freq(i) ! wavelength in Angstroms
	 end do
	else
	 do i=1,n
          if (wave(i) .eq. zero) then
	    write(*,*)'zero wavelength found at i=',i
	  endif
	  freq(i) = (c_kms*10D-3)/wave(i) ! wavelength in Angstroms
         end do
        end if

	return
	end

C-----------------------------------------------------------------------------
C subroutine to convert lambda array to ln(lambda) array if forward=true
C or convert ln(lambda) array to lambda array if forward=false
C
	subroutine lambda_to_lnlambda(wave,n,ltolnl)

	implicit none
	include 'constants.inc'

	integer*4 n
	real*8 wave(n)
	logical ltolnl  ! if T, lambda->lnlambda, if F, nlambda->lambda.

	integer*4 i

	if (ltolnl) then
	 do i=1,n
	  wave(i) = dlog(wave(i))
	 end do
	else 
	 do i=1,n
	  wave(i) = dexp(wave(i))
	 end do
	end if 

	return
	end

C-----------------------------------------------------------------------------
C subroutine to map array A onto array B using linear interpolation

	subroutine map(Awave,Aflux,AN,Bwave,Bflux,BN)

	implicit none
	include 'constants.inc'

	integer*4 AN,BN
	real*8 Awave(AN),Aflux(AN)
	real*8 Bwave(BN),Bflux(BN)

	integer*4 i,j

	i = 1
	j = 1

	do while (Bwave(j) .le. Awave(1))
	 j = j+1
	end do
	do while ((j .lt. BN) .and. (Bwave(j) .lt. Awave(AN)))
	 do while (.not.((Bwave(j).ge.Awave(i)).and.(Bwave(j).le.Awave(i+1))))
	  i = i+1
	 end do
	 Bflux(j) = Aflux(i) + (Aflux(i+1)-Aflux(i))*
     1   	  (Bwave(j) - Awave(i))/(Awave(i+1) - Awave(i))
	 j = j+1
	end do

	return 
	end


C-----------------------------------------------------------------------------
C subroutine to map array A onto array B using weighted averaging of fluxes.
C To determine the flux of a point in B, all points withing +/- dlam/2 are
C considered, as well as the the points longward and shortward of this
C range.  The fluxes are weighted by their wavelength distance (?) from
C the considered point.

	subroutine map2(Awave,Aflux,AN,Bwave,Bflux,BN,dlam)

	implicit none
	include 'constants.inc'

	integer*4 AN,BN
	real*8 Awave(AN),Aflux(AN)
	real*8 Bwave(BN),Bflux(BN)
	real*8 dlam,dlamover2

	integer*4 i,j,i_min,i_max

	real*8 fluxsum
	real*8 dlamsum

	dlamover2=dlam/two
	i = 1
	j = 1

	do while (Bwave(j) .le. Awave(1))
	  j = j+1
	end do
!
	i_min=1
	i_max=1
	do while ((j .le. BN) .and. (Bwave(j) .lt. Awave(AN)))
C
C        when determining the flux of a point on the new grid (B), we
C        wish to consider all points within that points resolution range
C        (+/- dlam/2).
C        However, in the case where no old grid-points lie within this
C        range on one side or the other, we need to find the next nearest
C        point outside that range.  The new flux is the average
C        over these points, weighted by distance from the new point.
C        Warning: possible danger of coefficient blowing up if it is
C        very close to wavelength of point of interest (or same as.)
C
	 do while(awave(i_min) .lt. bwave(j)-dlamover2)
	   i_min=i_min+1
	 end do
	 if(awave(i_min) .gt. bwave(j))i_min=i_min-1
	 do while(awave(i_max) .lt. bwave(j)+dlamover2)
	   i_max=i_max+1
	 end do
	 if(awave(i_max) .lt. bwave(j))i_max=i_max+1
!
	 fluxsum = zero
	 dlamsum = zero    ! this is actulally the sum of dlam^-1
	 i = i_min
	 do while (i .le. i_max)
	   dlamsum = dlamsum + one/(abs(Bwave(j) - Awave(i)))
	   fluxsum = fluxsum + Aflux(i)/abs(Bwave(j) - Awave(i))
	   i = i+1
	 end do
	 Bflux(j) = fluxsum/dlamsum
	 j = j+1
	end do

	return 
	end

C-----------------------------------------------------------------------------
C function which returns the index of the first array element *larger* than
C or equal to the supplied number

	function larger_index(list,n,value)

	implicit none

	integer*4 larger_index,n
	real*8 list(n),value

	integer*4 i

	if (list(n) .lt. value) then
	   write(*,*)'Error - no element larger than or equal to ',value
	   ! some error handling function here
	endif
	i = 1
	do while (value .ge. list(i) .and. i .lt. n)
	   i = i+1
	end do
	larger_index = i

	return
	end

C-----------------------------------------------------------------------------
C function which returns the index of the last array element *smaller* or
C equal to the supplied number

	function smaller_index(list,n,value)

	implicit none

	integer*4 smaller_index,n
	real*8 list(n),value

	integer*4 i

	if (list(1) .gt. value) then
	   write(*,*)'Error - no element smaller than or equal to ',value
	   ! some error handling function here
	endif
	i = 1
	do while (list(i) .le. value .and. i .le. n)
	   i = i+1
	end do
	smaller_index = i-1

	return
	end function
!
	END MODULE MOD_SMEAR_PG
