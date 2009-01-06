!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      SUBROUTINE SOLVE_CMF_FORMAL_V2(CHI,ETA,IP,NU_DNU,
     *                     BOUNDARY,B_NUE,dBDTAU,ND,NP,NC)
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Determine the transfer matrix along the given ip ray.  Formulation
! is the Short Characteristic solution from Olson & Kunsasz (87) JQSRT
! and Hauschildt (92) JQSRT
!
! Written   11-94 DLM
! Altered 2/21/96 DLM Updated to F90 standard
! Altered 5/9/97  DLM installed MODULE for rpz and allocation
!                     of many variables
! Altered 5/16/97 DLM installed three possible inner boundary conditions
!                     'hollow'    - hollow core with I+=I-
!                     'gray'      - diffusion approximation for gray atmosphere
!                                   Mihalas (1980) 237 p 574
!                     'diffusion' - diffusion approximation with I=B(nu,T)+mu*dBdtau(nu,T)
!                                   Mihalas (1978) p 51
! Altered 6/12/97 DLM Name changed from solve_rel_formal.f
! Altered 31/12/04 DJH : Cleaned
!--------------------------------------------------------------------
!
      USE MOD_SPACE_GRID_V2
      IMPLICIT NONE
!
! Grid size variables
!
      INTEGER :: NP,ND,NC
      INTEGER :: I,IP
!
! Opacity and emissivity variables
!
      real*8, dimension(nd) :: chi,eta
!
! Frequency variable
!
      real*8 :: nu_dnu
!
! Boundary conditions
!
      real*8 :: B_nue,dBdtau
!
      character*(*) :: boundary
!
! Local variables
!
      LOGICAL, PARAMETER :: L_TRUE=.TRUE.
      LOGICAL, PARAMETER :: L_FALSE=.FALSE.
!
      real*8 ee,e0,e1,e2,alpha,beta,gamma,t1
      real*8, dimension(nd) :: chi_tau
      real*8, dimension(nd) :: source_prime
!
      integer :: iz
      integer :: nzz
      integer :: luer,error_lu
      external error_lu
!
!--------------------------------------------------------------------
!
! Determine transfer variables.
!
!   Returns:
!     source_prime=(eta+alpha*Iprev)/(alpha+chi_prime)
!     chi_tau=alpha+chi_prime
!
!     where:
!       chi_prime=chi+3*b_*
!       alpha=freq(k)*b_*/(freq(k-1)-freq(k))
!       b_*=gamma*((1-mu_*^2)*beta/r+gamma^2*mu_**(mu_*+beta)dbetadr)
!
!         where: _* is _p or _m for ray in plus or minus direction
!
!	WRITE(6,*)ND,NU_DNU
!	WRITE(6,*)IP,RAY(IP)%NZ
!	WRITE(6,*)ETA(1),CHI(1)
!	WRITE(6,*)RAY(IP)%B_M(1)
!	WRITE(6,*)RAY(IP)%I_M_PREV(1)
!	WRITE(6,*)CHI_TAU(1)
!	WRITE(6,*)SOURCE_PRIME(1)
!
      CALL REL_VARIABLES(ND,CHI,ETA,NU_DNU,RAY(IP)%B_M,RAY(IP)%I_M_PREV,CHI_TAU,SOURCE_PRIME)
!
! Determine "optical depth" from outside to center (s_m,mu_m)
!
      CALL OPTDEPTH_V2(CHI_TAU,RAY(IP)%NZ,IP,L_FALSE)
!
! Calculate transfer in inward direction, I- (mu=-1)
!
      RAY(IP)%I_m(1)=0.0d0
!
!---------------------------------------------------------------
!
! Formal integral from OUTSIDE to INSIDE
!
!---------------------------------------------------------------
!
      nzz=ray(ip)%nz
      do iz=2,ray(ip)%nz-1
!
        t1=dtau(iz-1)
        if(t1 .gt. 0.01D0)then
          ee=exp(-dtau(iz-1))
          e0=1.0d0-ee
          e1=dtau(iz-1)-e0
          e2=dtau(iz-1)*dtau(iz-1)-2.0d0*e1
        else
          e2=t1*t1*t1*(1.0D0-0.25D0*t1*(1-0.2D0*t1*(1.0D0-t1/6.0D0*
     1              (1.0D0-t1/7.0D0))))/3.0D0
          e1=0.5d0*(t1*t1-e2)
          e0=t1-e1
          ee=1.0d0-e0
        end if
!
        alpha=e0+(e2-(dtau(iz)+2.0d0*dtau(iz-1))*e1)/(dtau(iz-1)*(dtau(iz)+dtau(iz-1)))
        beta=((dtau(iz)+dtau(iz-1))*e1-e2)/(dtau(iz-1)*dtau(iz))
        gamma=(e2-dtau(iz-1)*e1)/(dtau(iz)*(dtau(iz)+dtau(iz-1)))
     	t1=alpha*source_prime(iz-1)+beta*source_prime(iz)+gamma*source_prime(iz+1)
        if(t1 .lt. 0)then
	  t1=source_prime(iz)+dtau(iz)*(source_prime(iz)-source_prime(iz-1))/dtau(iz-1)
     	  t1=alpha*source_prime(iz-1)+beta*source_prime(iz)+gamma*t1
	end if
	ray(ip)%I_m(iz)=ray(ip)%I_m(iz-1)*ee+t1
!
      end do
!
! If only 2 points along the ray or at inner boundary, must use linear
! source function
!
      t1=dtau(nzz-1)
      if(t1 .gt. 0.01d0)then
        ee=exp(-dtau(nzz-1))
        e0=1.0d0-ee
        e1=dtau(nzz-1)-e0
        e2=dtau(nzz-1)*dtau(nzz-1)-2.0d0*e1
      else
        e2=t1*t1*t1*(1.0D0-0.25D0*t1*(1-0.2D0*t1*(1.0D0-t1/6.0D0*
     1              (1.0D0-t1/7.0D0))))/3.0D0
        e1=0.5d0*(t1*t1-e2)
        e0=t1-e1
        ee=1.0d0-e0
      end if
!
      alpha=e0-e1/dtau(nzz-1)
      beta=e1/dtau(nzz-1)
      ray(ip)%I_m(nzz)=ray(ip)%I_m(nzz-1)*ee+
     *     alpha*source_prime(nzz-1)+beta*source_prime(nzz)
!
!---------------------------------------------------------------
!
! Calculate transfer in outward direction, I+ (mu=1)
!
! At inner boundary use I+ = I- for core and non-core rays
!
      if((ip.gt.nc+1).or.(boundary.eq.'hollow'))then
        ray(ip)%I_p(nzz)=ray(ip)%I_m(nzz)
      elseif((boundary.eq.'gray').or.(boundary.eq.'diffusion'))then
        ray(ip)%I_p(nzz)=B_nue+ray(ip)%mu_p(nzz)*dBdtau
      else
	luer=error_lu()
        write(luer,*)' Unknown inner boundary condition in SOLVE_CMF_FORMAL'
        write(luer,*)'Passed boundary condition is ',trim(boundary)
        stop
      end if
!
! Determine transfer variables.
!
      call rel_variables(nd,chi,eta,nu_dnu,ray(ip)%b_p,
     *       ray(ip)%I_p_prev,chi_tau,source_prime)
!
! Determine "optical depth" from outside to center (s_m,mu_m)
!
      CALL OPTDEPTH_V2(CHI_TAU,RAY(IP)%NZ,IP,L_TRUE)
!
!---------------------------------------------------------------
!
! Formal integral from INSIDE to OUTSIDE
!
!---------------------------------------------------------------
!
      do iz=nzz-1,2,-1
!
        t1=dtau(iz)
        if(t1 .gt. 0.01d0)then
          ee=exp(-dtau(iz))
          e0=1.0d0-ee
          e1=dtau(iz)-e0
          e2=dtau(iz)*dtau(iz)-2.0d0*e1
        else
           e2=t1*t1*t1*(1.0D0-0.25D0*t1*(1-0.2D0*t1*(1.0D0-t1/6.0D0*
     1              (1.0D0-t1/7.0D0))))/3.0D0
           e1=0.5d0*(t1*t1-e2)
           e0=t1-e1
           ee=1.0d0-e0
         end if
!
        alpha=(e2-dtau(iz)*e1)/
     *       (dtau(iz-1)*(dtau(iz)+dtau(iz-1)))
        beta=((dtau(iz)+dtau(iz-1))*e1-e2)/(dtau(iz-1)*dtau(iz))
        gamma=e0+(e2-(dtau(iz-1)+2.0d0*dtau(iz))*e1)/(dtau(iz)*(dtau(iz)+dtau(iz-1)))
        t1=alpha*source_prime(iz-1)+beta*source_prime(iz)+gamma*source_prime(iz+1)
        if(t1 .lt. 0)then
	  t1=source_prime(iz)+dtau(iz-1)*(source_prime(iz)-source_prime(iz+1))/dtau(iz)
     	  t1=alpha*t1+beta*source_prime(iz)+gamma*source_prime(iz+1)
	end if
        ray(ip)%I_p(iz)=ray(ip)%I_p(iz+1)*ee+t1
!
      end do
!
! If only 2 points along the ray or at outer boundary, must use linear
! source function
!
      t1=dtau(1)
      if(t1 .gt. 0.01d0)then
        ee=exp(-dtau(1))
        e0=1.0d0-ee
        e1=dtau(1)-e0
        e2=dtau(1)*dtau(1)-2.0d0*e1
      else
        e2=t1*t1*t1*(1.0D0-0.25D0*t1*(1-0.2D0*t1*(1.0D0-t1/6.0D0*
     1              (1.0D0-t1/7.0D0))))/3.0D0
        e1=0.5d0*(t1*t1-e2)
        e0=t1-e1
        ee=1.0d0-e0
      end if
!
      beta=e1/dtau(1)
      gamma=e0-e1/dtau(1)
      ray(ip)%I_p(1)=ray(ip)%I_p(2)*ee+beta*source_prime(1)+gamma*source_prime(2)
!
      return
      end
