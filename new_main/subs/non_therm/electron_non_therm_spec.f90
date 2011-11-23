	SUBROUTINE ELECTRON_NON_THERM_SPEC(ND)
	USE MOD_NON_THERM
	USE MOD_CMFGEN
	IMPLICIT NONE
!
!	INTEGER NUM_IONS
	INTEGER ND
!
	REAL*8, SAVE, ALLOCATABLE :: LELEC(:)
	REAL*8, SAVE, ALLOCATABLE :: SOURCE(:)
	REAL*8, SAVE, ALLOCATABLE :: RHS(:)
	REAL*8, SAVE, ALLOCATABLE :: DETAI_DE(:)
	REAL*8, SAVE, ALLOCATABLE :: MAT(:,:)
	REAL*8, SAVE, ALLOCATABLE :: Qnn(:)
	REAL*8, SAVE, ALLOCATABLE :: POP_VEC(:)
!
	INTEGER, SAVE, ALLOCATABLE :: INDX(:)
	LOGICAL, SAVE, ALLOCATABLE :: DO_THIS_ION(:)
	logical, save, allocatable :: do_and_1st_time(:)
	integer, save, allocatable :: nlow_maxs(:)
!
	REAL*8 EKT
	REAL*8 EKTP
	REAL*8 EMIN
	REAL*8 EMAX
	REAL*8 E_INIT
	REAL*8 dE
!
	REAL*8 T1,T2
	REAL*8 SIG1,SIG2
	REAL*8 XCROSS
	REAL*8 XION_POT
	REAL*8 NATOM
!
	INTEGER GET_INDX_DP
	REAL*8 INTSIGC
	EXTERNAL INTSIGC,GET_INDX_DP
!
	REAL*8, PARAMETER :: XKT_MIN=1.0D0
	REAL*8, PARAMETER :: XKT_MAX=1.0D+03
	REAL*8, PARAMETER :: DELTA_ENR_SOURCE=30.0D0
	REAL*8, PARAMETER :: Hz_to_eV=13.60D0/3.2897D0
	LOGICAL, PARAMETER :: VERBOSE_OUTPUT=.FALSE.
!
	INTEGER ISPEC
	INTEGER ID
	INTEGER IT
	INTEGER IKT
	INTEGER IKT0
	INTEGER IKTP
	INTEGER IKTPN
	INTEGER DPTH_INDX
	INTEGER IOS
	INTEGER I, J, K
	INTEGER IST
	INTEGER NL,NUP
	INTEGER L
	INTEGER MAX_LOW_LEV
!
	INTEGER, PARAMETER :: LU_ER=6
	INTEGER, PARAMETER :: LU_TH=8
	integer, parameter :: lu_bethe = 100
	LOGICAL, SAVE :: FIRST_TIME=.TRUE.
	LOGICAL, PARAMETER :: L_TRUE=.TRUE.
	LOGICAL, PARAMETER :: L_FALSE=.FALSE.
!
	LOGICAL INJECT_DIRAC_AT_EMAX
	LOGICAL INCLUDE_EXCITATION
	LOGICAL INCLUDE_IONIZATION
	CHARACTER(LEN=10) SOURCE_TYPE
	CHARACTER(LEN=40) TMP_NAME
!
	SOURCE_TYPE='CONSTANT'       !or 'BELL SHAPE'
	INJECT_DIRAC_AT_EMAX=.TRUE.
	INCLUDE_IONIZATION=.TRUE.    !.TRUE.
	INCLUDE_EXCITATION=.TRUE.    !.TRUE.
!
	WRITE(LU_ER,*)' '
	WRITE(LU_ER,*)'Entering ELECTRON_NON_THERM_SPEC'
	OPEN(UNIT=LU_TH,FILE='NON_THERM_SPEC_INFO',STATUS='UNKNOWN')
!	CALL SET_LINE_BUFFERING(LU_TH)
!
! Allocate arrays that will be used for each depth:.
!
! The following vectors/arrays are stored in MOD_NON_THERM
!
	IF(.NOT. ALLOCATED(XKT))THEN
	  WRITE(LU_TH,*)'Allocating memory to describe non-thermal electron distribution'
	  ALLOCATE (XKT(NKT),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (dXKT(NKT),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (YE(NKT,ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (FRAC_ELEC_HEATING(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (FRAC_ION_HEATING(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (FRAC_EXCITE_HEATING(ND),STAT=IOS)
	  IF (IOS.NE.0) THEN
            WRITE(LU_ER,*) 'Error allocating XKT',IOS
	    STOP
	  END IF
	  XKT(1:NKT)=0.0D0; dXKT(1:NKT)=0.0D0
	END IF
!
! The following vectors/arrrays are local:
!
! LELEC : Energy loss through Coulomb interaction
! SOURCE : Electron source
!
	WRITE(LU_TH,*)'Checking MAT memory allocation'
	IF(.NOT. ALLOCATED(MAT))THEN
	  WRITE(LU_TH,*)'Allocating memory in ELECTRON_NON_THERMAL_SPEC'
	  ALLOCATE (LELEC(NKT),STAT=IOS) 
	  IF(IOS .EQ. 0)ALLOCATE (SOURCE(NKT),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (RHS(NKT),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (Qnn(NKT),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (dETAI_dE(NKT),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (MAT(NKT,NKT),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (POP_VEC(NUM_THD))
	  IF(IOS .EQ. 0)ALLOCATE (INDX(NUM_THD))
	  IF(IOS .EQ. 0)ALLOCATE (DO_THIS_ION(NUM_IONS))
	  if(ios .eq. 0)allocate (do_and_1st_time(num_ions))
	  if(ios .eq. 0)allocate (nlow_maxs(num_ions))
	  IF(IOS .NE. 0)THEN
	    WRITE(LU_ER,*) 'Error allocating LELEC, SOURCE, etc in ELECTRON_NON_TERMAL_SPEC',IOS
	    STOP
	  END IF
	END IF
!
! Set the energy array. Use a linear scaling since we need to resolve the high-energy source
! but also the transition near the ionization/excitation potentials at about 10eV.
!
	IF(FIRST_TIME)THEN
!
	  WRITE(LU_TH,*)'Setting XKT array'
	  CALL SET_XKT_ARRAY(XKT_MIN,XKT_MAX,NKT,XKT,dXKT,'lin')
!
! If needed, allocate memory for collisional ioinzation cross-sections.
!
	  IOS=0
	  DO IT=1,NUM_THD
	    IF(THD(IT)%PRES)ALLOCATE (THD(IT)%CROSS_SEC(NKT),STAT=IOS) 
	    IF(IOS.NE.0) THEN
              WRITE(LU_ER,*) 'Error allocating THD%CROSS_SEC',IOS
	      STOP
	    END IF
	  END DO
	  WRITE(LU_TH,*)'Computing cross-sections'
	  CALL ARNAUD_CROSS_V3()
!
	  I=0
	  DO IT=1,NUM_THD
	    IF(THD(IT)%PRES)I=I+1
	  END DO
	  WRITE(LU_TH,*)'Number of active species in ELECTRON_NON_TERMAL_SPEC is',I
	END IF
!
	FIRST_TIME=.FALSE.
	FRAC_ELEC_HEATING=0.0D0
	FRAC_ION_HEATING=0.0D0
	FRAC_EXCITE_HEATING=0.0D0
	YE=0.0D0
!
! Set how the lectrons are ejected at high enegry. This will be the smae for all
! depths. We can either iject at EMAX, or use a Bell-shaped curve near EMAX.
! E_INIT will be a normalisation constant to yield fractions for ionization etc.
!
	WRITE(LU_TH,*)'Constructing SOURCE'
	IF (INJECT_DIRAC_AT_EMAX) THEN ! zero otherwise
	  E_INIT = XKT_MAX
	ELSE
	  T1 = 0.0D0
	  T2 = DELTA_ENR_SOURCE
       	  IF (SOURCE_TYPE .EQ. 'CONSTANT')THEN
	    DO IKT=1,NKT
	      IF (XKT(IKT) .GE. (XKT_MAX-T2)) THEN
	        SOURCE(IKT) = 1.0D0
	        T1 = T1 + SOURCE(IKT)*dXKT(IKT)
              END IF
	    END DO
	  ELSE IF(SOURCE_TYPE .EQ. 'BELL_SHAPE')THEN
	    DO IKT=1,NKT
	      IF (XKT(IKT) .GE. (XKT_MAX-T2)) THEN
                SOURCE(IKT) = 6.0D0 / T2**3 * (T2**2/4.D0 - (XKT(IKT)-XKT_MAX+T2/2.0D0)**2)
	        T1 = T1 + SOURCE(IKT)*dXKT(IKT)
              END IF
	    END DO
          END IF
	  SOURCE(:) = SOURCE(:) / T1
!
	  T1 = 0.0D0
	  T2 = 0.0D0
	  DO IKT=1,NKT
	    T1 = T1 + SOURCE(IKT)*dXKT(IKT)
	    T2 = T2 + XKT(IKT)*SOURCE(IKT)*dXKT(IKT)
	  END DO
	  E_INIT = T2
	END IF
!
!	open(unit=lu_bethe,file='bethe_cross_chk',status='unknown')
	do_and_1st_time(:)=l_false
	DO DPTH_INDX=1,ND
!
	  WRITE(LU_TH,'(/,X,A,I3)')'Starting depth index: ',DPTH_INDX
!
! Decide on which species we will handle when computing the electron distribution.
! For each route, we sum up over all levels of that term. This approach will
! be independent of the SL assignments.
!
	  DO IT=1,NUM_THD
	    THD(IT)%N_ATOM=0.0D0
	    IF(THD(IT)%PRES)THEN
	      ID=THD(IT)%LNK_TO_ION
	      DO J=1,THD(IT)%N_STATES
	        K=THD(IT)%ATOM_STATES(J)
	        THD(IT)%N_ATOM=THD(IT)%N_ATOM+ATM(ID)%XzV_F(K,DPTH_INDX)
	      END DO
	    END IF
	    POP_VEC(IT)=THD(IT)%N_ATOM
	  END DO
!
! DO_THIS_ION and THD()%DO_THIS_ION_ROUTE contain similar information.
!   DO_THIS_ION is used for deciding by ION
!   DO_THIS_ION_ROUTE is used for deciding by the ionization route.
!
	  WRITE(LU_TH,'(X,A)')'Determining the top 10 ionization routes at this depth'
	  WRITE(LU_TH,'(2X,A,4X,A,7X,A,9X,A)')'IT','ID','Ion','Pop'
	  CALL INDEXX(NUM_THD,POP_VEC,INDX,L_FALSE)
	  DO_THIS_ION(:)=L_FALSE
	  THD(:)%DO_THIS_ION_ROUTE=L_FALSE
	  DO IT=1,10
	    ID=THD(INDX(IT))%LNK_TO_ION
	    THD(INDX(IT))%DO_THIS_ION_ROUTE=.TRUE.
	    DO_THIS_ION(ID)=.TRUE.
	    do_and_1st_time(id)=.true.
	    WRITE(LU_TH,'(I4,2X,I4,4X,A6,ES12.3)')IT,ID,TRIM(ION_ID(ID)),POP_VEC(INDX(IT))
	  END DO
!
! Building the upper-diagonal matrix mat(1:nkt,1:nkt)
!
	  WRITE(LU_TH,*)'Constructing MAT'
!
! Get the energy losses through Coulomb interaction
!
	  WRITE(LU_TH,*)'Electron density=',ED(DPTH_INDX)  !ED(DPTH_INDX)
	  CALL GET_LELEC(LELEC,XKT,NKT,ED(DPTH_INDX))  !,ED(DPTH_INDX))
	  IF(VERBOSE_OUTPUT)THEN
	    DO IKT=1,NKT
	      WRITE(100,*)IKT,XKT(IKT),LELEC(IKT)
	    END DO
	  END IF
!
! Initialize MAT, and the add Diagonal term for Coulomb interaction with 
! thermal electrons.
!
	  MAT=0.0D0
	  DO IKT=1,NKT
	    MAT(IKT,IKT) = LELEC(IKT)
	  END DO
!
	  CALL TUNE(1,'MAT_ION')
	  IF(INCLUDE_IONIZATION)THEN
!
! Loop over all species and ionization stages
!
	    DO IT=1,NUM_THD
	      IF(THD(IT)%PRES .AND. THD(IT)%DO_THIS_ION_ROUTE)THEN
	        NATOM=THD(IT)%N_ATOM
                XION_POT = THD(IT)%ION_POT
	        WRITE(LU_TH,'(A,I3,A,F7.2,A)')'Ionization potential for IT=',IT,' is ',XION_POT,' eV'
!
	        IF(XION_POT .LT. XKT(1))THEN
	           IST=1
	        ELSE IF(XION_POT .GT. XKT(NKT))THEN
	           IST=NKT+1
	        ELSE
	           IST=GET_INDX_DP(XION_POT,XKT,NKT)
	        END IF
! 
!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(IKT,IKTP,EKT,EKTP,EMIN,EMAX,XCROSS)
	        DO IKTP=IST,NKT
	          EkTP = XKT(IKTP)
                  XCROSS = THD(IT)%CROSS_SEC(IKTP)
!
	          DO IKT=1,IKTP
	            IF(XCROSS .EQ. 0.0D0)EXIT
	            EKT = XKT(IKT)
	            EMIN = MAX(EKTP-EKT,XION_POT)
	            EMAX = MIN(0.5D0*(EKTP+XION_POT),XKT(NKT))
	            IF(EMAX .GT. EMIN)THEN
	              MAT(IKT,IKTP) = MAT(IKT,IKTP) + &
                            NATOM*dXKT(IKTP)*XCROSS*INTSIGC(EKTP,XION_POT,EMIN,EMAX)
	            END IF
!
                    IF (EKTP .GT. (2.0D0*EKT+XION_POT)) THEN
                      EMIN = EKT+XION_POT
                      EMAX = MIN(0.5d0*(EKTP+XION_POT),XKT(NKT))
                      MAT(IKT,IKTP) = MAT(IKT,IKTP) &
                       - NATOM*dXKT(IKTP)*XCROSS*INTSIGC(EKTP,XION_POT,EMIN,EMAX)
                    END IF
	          END DO
!
	        END DO
!$OMP END PARALLEL DO
	      END IF    !Species present
	    END DO	!Looping over ionization stage / shell
	  END IF
	  CALL TUNE(2,'MAT_ION')
!
! Collisional excitations.
!
	CALL TUNE(1,'MAT_EXCITE')
	IF (INCLUDE_EXCITATION)THEN
	  WRITE(LU_TH,*)'Beginning excitation section'
!	  open(unit=lu_bethe,file='bethe_cross_chk',status='unknown',action='write')
	  DO ID=1,NUM_IONS
	    IF(ATM(ID)%XzV_PRES .AND. DO_THIS_ION(ID))THEN
!
! Determine levels from which excitation will occur.
!
	      T1=SUM(ATM(ID)%XzV_F(:,DPTH_INDX))
	      MAX_LOW_LEV=ATM(ID)%NXzV_F
	      DO I=ATM(ID)%NXzV_F,1,-1
	        IF(ATM(ID)%XzV_F(I,DPTH_INDX)/T1 .GT. 1.0D-04)EXIT
	        MAX_LOW_LEV=I
	      END DO
	      nlow_maxs(id)=max(nlow_maxs(id),max_low_lev)
!
	      DO J=2,ATM(ID)%NXzV_F
	        DO I=1,MIN(J-1,MAX_LOW_LEV)
	          NL=I; NUP=J
!	          write(lu_bethe,'(x,I5,2X,I4,4X,A6,ES16.8)')dpth_indx,ID,TRIM(ION_ID(ID)), &
!	                                                atm(id)%xzv_f(nl,dpth_indx)
!	          CALL BETHE_APPROX(Qnn,NL,NUP,XKT,dXKT,NKT,ID,DPTH_INDX,FAST_BETHE_METHOD)
	          CALL BETHE_APPROX_V3(Qnn,NL,NUP,XKT,dXKT,NKT,ID,DPTH_INDX)
!
!
!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(IKT,IKT0,IKTP,IKTPN,EKT,EMAX)
	          DO IKT=1,NKT
	            EKT = XKT(IKT)
	            IKT0 = IKT
!
! Need to find the upper bound iktpn for which ektp = ekt + Excitation_Energy
!
	            EMAX = EKT+Hz_TO_eV*(ATM(ID)%EDGEXZV_F(I)-ATM(ID)%EDGEXZV_F(J))
	            IKTPN=IKT0 
	            IF (EMAX .GT. XKT(NKT-1)) then
	              IKTPN = NKT
	            ELSE
	              DO WHILE (XKT(IKTPN) .LT. EMAX)
	                IKTPN = IKTPN + 1
	              END DO
	            END IF
!
	            DO IKTP=IKT0,IKTPN
	               MAT(IKT0,IKTP) = MAT(IKT0,IKTP) + Qnn(IKTP)
	            END DO
!
	          END DO	!Loop over energy
!OMP END PARALLEL DO
!
	        END DO		!Loop over lower level
	      END DO		!Loop over upper level
	    END IF		!Ionization stage present
	  END DO		!Loop over ionization stage
!	  close(lu_bethe)
	END IF			!Include excitations?
	CALL TUNE(2,'MAT_EXCITE')
!
! Check no diagonal term of mat is zero
!
	 DO IKT=1,NKT
	   IF (MAT(IKT,IKT) .EQ. 0.0D0) THEN
	     WRITE(LU_ER,*) 'Zero diagonal element in mat ',IKT,SOURCE(IKT)
	     STOP
	   END IF
	 END DO
!
! We now build the right-hand side and the matrix corresponding to
! Eq.7 of KF92.
!
	 IF (INJECT_DIRAC_AT_EMAX) THEN
	   RHS(1:NKT) = 1.0D0
	 ELSE
	   RHS(NKT) = 0.0D0
	   IF (SOURCE_TYPE .EQ. 'CONSTANT') RHS(NKT) = SOURCE(NKT) * dXKT(NKT)
	   DO IKT=NKT-1,1,-1
	     RHS(IKT) = RHS(IKT+1) + SOURCE(IKT) * dXKT(IKT)
	   END DO
	 END IF
!
! Solve for the degradation function YE
!
	 YE(NKT,DPTH_INDX)=0.0D0
	 IF(INJECT_DIRAC_AT_EMAX)YE(NKT,DPTH_INDX)=1.0D0/LELEC(NKT)
	 DO IKT=NKT-1,1,-1
	   T1 = 0.0D0
	   DO IKTP=IKT+1,NKT
	      T1 = T1 + MAT(IKT,IKTP)*YE(IKTP,DPTH_INDX)
	   END DO
	   YE(IKT,DPTH_INDX) = (RHS(IKT) - T1) / MAT(IKT,IKT) ! diagonal term is never zero (lelec contrib.)
	 END DO
!
! Compute electron heating
!
	  DO IKT=1,NKT
	    FRAC_ELEC_HEATING(DPTH_INDX)=FRAC_ELEC_HEATING(DPTH_INDX) + YE(IKT,DPTH_INDX)*LELEC(IKT)*dXKT(IKT)
	  END DO
	  FRAC_ELEC_HEATING(DPTH_INDX)=FRAC_ELEC_HEATING(DPTH_INDX)/E_INIT
	  WRITE(LU_TH,*)'Fraction into electron heating is',FRAC_ELEC_HEATING(DPTH_INDX)
	  WRITE(LU_TH,*)'Fraction into electron heating is',XKT(1)*LELEC(1)*YE(1,DPTH_INDX)/E_INIT
	  FRAC_ELEC_HEATING(DPTH_INDX)=FRAC_ELEC_HEATING(DPTH_INDX) + XKT(1)*LELEC(1)*YE(1,DPTH_INDX)/E_INIT
!
! Compute the ionization and heating contributions
!
	  CALL TUNE(1,'CHK_ION')
	  WRITE(LU_TH,*)'Checking ionization energy'
	  WRITE(LU_TH,'(3X,A,2X,A,4X,A,2X,A,3X,A)')'IT','IST','SIG(IST)','SIG(IST+1)','dE(ION)/E'
	  IF (INCLUDE_IONIZATION) THEN
	    DO IT=1,NUM_THD
	      IF(THD(IT)%PRES .AND. THD(IT)%DO_THIS_ION_ROUTE)THEN
	        NATOM=THD(IT)%N_ATOM
                XION_POT = THD(IT)%ION_POT
!
	        IF(XION_POT .LT. XKT(1))THEN
	           IST=1
	           SIG1=THD(IT)%CROSS_SEC(1); SIG2=0.0D0
	        ELSE IF(XION_POT .GT. XKT(NKT))THEN
	           IST=NKT+1
	           SIG1=0.0D0; SIG2=0.0D0
	        ELSE
	           IST=GET_INDX_DP(XION_POT,XKT,NKT)
	           SIG1=THD(IT)%CROSS_SEC(IST); SIG2=THD(IT)%CROSS_SEC(IST+1)
	        END IF
	        T1 = 0.0D0
	        DO IKTP=IST,NKT
	          XCROSS = THD(IT)%CROSS_SEC(IKTP)
	          DETAI_DE(IKTP) = DETAI_DE(IKTP) + &
	                NATOM * XION_POT / E_INIT * YE(IKTP,DPTH_INDX)*XCROSS
	          T1 = T1 + YE(IKTP,DPTH_INDX)*dXKT(IKTP)*XCROSS
	        END DO
	        WRITE(LU_TH,'(2I5,3E12.2)')IT,IST,SIG1,SIG2,NATOM * XION_POT * T1 / E_INIT
	        FRAC_ION_HEATING(DPTH_INDX) = FRAC_ION_HEATING(DPTH_INDX) + NATOM * XION_POT * T1 / E_INIT
	      END IF
	    END DO
	    WRITE(LU_TH,*)'Fraction into ionization heating is',FRAC_ION_HEATING(DPTH_INDX)
	  END IF
	  CALL TUNE(2,'CHK_ION')
!
	CALL TUNE(1,'CHK_EXCITE')
	IF(INCLUDE_EXCITATION)THEN
	  DO ID=1,NUM_IONS
	    IF(ATM(ID)%XzV_PRES .AND. DO_THIS_ION(ID))THEN
	      DO J=2,ATM(ID)%NXzV_F
	        DO I=1,MIN(J-1,MAX_LOW_LEV)
	          NL=I; NUP=J
!	          CALL BETHE_APPROX(Qnn,NL,NUP,XKT,dXKT,NKT,ID,DPTH_INDX,L_TRUE)
	          CALL BETHE_APPROX_V3(Qnn,NL,NUP,XKT,dXKT,NKT,ID,DPTH_INDX)
	          dE=Hz_TO_eV*(ATM(ID)%EDGEXZV_F(I)-ATM(ID)%EDGEXZV_F(J))
	          DO IKT=1,NKT
	            FRAC_EXCITE_HEATING(DPTH_INDX)=FRAC_EXCITE_HEATING(DPTH_INDX)+dE*Qnn(IKT)*YE(IKT,DPTH_INDX)
	          END DO
!	          write(lu_bethe,*)''
!	          write(lu_bethe,'(x,i3,x,i4,x,i5,x,i5)')id,dpth_indx,nl,nup
!	          write(lu_bethe,'(1X,1P8E16.7)')(Qnn(ikt)/atm(id)%xzv_f(nl,dpth_indx),ikt=1,nkt)
	        END DO
	      END DO
	    END IF
	  END DO		!Loop over ionization stage
	END IF			!Include excitations?
	CALL TUNE(2,'CHK_EXCITE')
	FRAC_EXCITE_HEATING(DPTH_INDX)=FRAC_EXCITE_HEATING(DPTH_INDX)/E_INIT
	WRITE(LU_TH,*)'Fraction into excitation heating is',FRAC_EXCITE_HEATING(DPTH_INDX)
!
	WRITE(LU_TH,*)'Total heating is ', FRAC_ELEC_HEATING(DPTH_INDX) + &
	                               FRAC_ION_HEATING(DPTH_INDX)+ &
	                               FRAC_EXCITE_HEATING(DPTH_INDX)
!
!	  CALL DP_CURVE(NKT,XKT,YE(1,DPTH_INDX))
!	  CALL GRAMON_PGPLOT('E(keV)','Y(E)',' ',' ')
!
	END DO 		!loop over depth
!	open(unit=lu_bethe,file='bethe_cross_chk',status='unknown')
!	  if(include_excitation)then
!	    do id=1,num_ions
!	      if(atm(id)%xzv_pres .and. do_and_1st_time(id))then
!	        do_and_1st_time(id) = .false.
!	        do j=2,atm(id)%nxzv_f
!	          do i=1,min(j-1,nlow_maxs(id))
!	            nl=i; nup=j
!	            CALL bethe_cross(Qnn,NL,NUP,XKT,dXKT,NKT,ID)
!	            write(lu_bethe,*)''
!	            write(lu_bethe,'(x,a6,x,i3,x,i5,x,i5)')trim(ion_id(id)),id,nl,nup
!	            write(lu_bethe,'(1X,1P8E16.7)')(Qnn(ikt),ikt=1,nkt)
!	          end do
!	        end do
!	      end if
!	    end do
!	  end if
!	close(unit=lu_bethe)
!
!
	WRITE(LU_TH,'(/,X,A,/)')'Summary of energy channels'
	WRITE(LU_TH,'(4(5X,A))'),' dE(ELEC)/E','  dE(ION)/E','  dE(EXC)/E','dE(TOTAL)/E'
	DO I=1,ND
	  WRITE(LU_TH,'(4ES16.6)')FRAC_ELEC_HEATING(I),FRAC_ION_HEATING(I),FRAC_EXCITE_HEATING(I), &
	  FRAC_ELEC_HEATING(I)+FRAC_ION_HEATING(I)+FRAC_EXCITE_HEATING(I)
	END DO
!
! Below we compute the number of non thermal electrons arising from the input of non-thermal energy.
! This is the number on electrons per 1eV of energy input.
!
	WRITE(LU_TH,*)' '
	WRITE(LU_TH,*)'Integral of normalized YE as a function of depth'
	WRITE(LU_TH,'(1X,A,5X,A,10X,A)')'Depth','Int(YE)','Ne'
	DO DPTH_INDX=1,ND
	  T1=0.0D0
	  DO IKT=1,NKT
	   T1=T1+YE(IKT,DPTH_INDX)/SQRT(XKT(IKT))*dXKT(IKT)
	  END DO
	  T1=T1*SQRT(9.109389D-28/2.0D0/1.602177D-12)/E_INIT
	  WRITE(LU_TH,'(I6,2ES12.3)')DPTH_INDX,T1,ED(DPTH_INDX)
	END DO
!
! Close diagnostic file, forcing all output.
!
	CLOSE(UNIT=LU_TH)
	WRITE(LU_ER,*)'Exiting ELECTRON_NON_THERM_SPEC'
!
! Before leaving the routine we normalize YE by the E_INIT. We can later scale by the actual energy input.
!
	YE=YE/E_INIT
!
! Output the degradation spectrum
!
	open(unit=lu_bethe,file='NON_THERM_DEGRADATION_SPEC',status='unknown')
	write(lu_bethe,'(1x,A,I5)')'Number of energy grid: ',NKT
	write(lu_bethe,'(1x,A,I5)')'Number of depth: ',ND
	write(lu_bethe,*)''
	write(lu_bethe,'(A)')'Energy grids (eV)'
	write(lu_bethe,'(1x,1p8e16.7)')(XKT(IKT),IKT=1,NKT)
	write(lu_bethe,*)''
	write(lu_bethe,'(A)')'Source'
	write(lu_bethe,'(1x,1p8e16.7)')(SOURCE(IKT),IKT=1,NKT)
	write(lu_bethe,*)''
	write(lu_bethe,'(A)')'Lelec'
	write(lu_bethe,'(1x,1p8e16.7)')(LELEC(IKT),IKT=1,NKT)
	write(lu_bethe,*)''
	do dpth_indx=1,nd
	  write(lu_bethe,'(A,I5)')'Depth: ',dpth_indx
	  write(lu_bethe,'(1x,1p8e16.7)')(YE(IKT,DPTH_INDX),IKT=1,NKT)
	  write(lu_bethe,*)''
	end do
	close(unit=lu_bethe)
!
	RETURN
	END
