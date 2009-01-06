	SUBROUTINE GET_MOMS_NON_REL(R, V, FREQ, N_TYPE, ND)
	USE MOD_RAY_MOM_STORE
	IMPLICIT NONE
!
	INTEGER ND
	REAL*8 R(ND)
	REAL*8 V(ND)
	REAL*8 FREQ
	CHARACTER(LEN=*) N_TYPE
!
! Local vectors
!
	REAL*8 TMP_VEC(1:ND)
	REAL*8 RSQ_JNU(1:ND)
	REAL*8 RSQ_HNU_MID(1:ND)
	REAL*8 RSQ_NNU_MID(1:ND)
!
	REAL*8 T1,T2
	INTEGER, PARAMETER :: IONE=1
	INTEGER I
	INTEGER NS
	INTEGER NDM1
	INTEGER LUER,ERROR_LU
	EXTERNAL ERROR_LU
!
	NS=ND_STORE
	IF(.NOT. ALLOCATED(K_ON_J))THEN
	   ALLOCATE (RMID(ND))
!
	   ALLOCATE (K_ON_J(ND))
	   ALLOCATE (NMID_ON_HMID(ND))
	   ALLOCATE (NMID_ON_J(ND))
!
	   ALLOCATE (K_ON_J_PREV(ND))
	   ALLOCATE (NMID_ON_HMID_PREV(ND))
	   ALLOCATE (NMID_ON_J_PREV(ND))
!
	   ALLOCATE (K_ON_J_SAVE(ND))
	   ALLOCATE (NMID_ON_HMID_SAVE(ND))
	   ALLOCATE (NMID_ON_J_SAVE(ND))
	 END IF
	 DO I=1,ND-1
	   RMID(I)=0.5D0*(R(I)+R(I+1))
	 END DO
!
! Compute FEDD=K/J. Since J & K are always computed at the nodes,
! this calculaton is independent of HN_DEF_AT_NODES.
!
! We only need RSQ_JNU if we are computing N(mid)/J(node)
!
	IF(ND .NE. NS)THEN
	  IF(N_TYPE .NE. 'G_ONLY')THEN
	    TMP_VEC(1:NS)=R_STORE(1:NS)*R_STORE(1:NS)*JNU_STORE(1:NS)
	    CALL MON_INTERP(RSQ_JNU,ND,IONE,R,ND,TMP_VEC,NS,R_STORE,NS)
	  END IF
	  TMP_VEC(1:NS)=KNU_STORE(1:NS)/JNU_STORE(1:NS)
	  CALL MON_INTERP(K_ON_J,ND,IONE,R,ND,TMP_VEC,NS,R_STORE,NS)
	ELSE
	  IF(N_TYPE .NE. 'G_ONLY')THEN
	    RSQ_JNU(1:NS)=JNU_STORE(1:NS)*R_STORE(1:NS)*R_STORE(1:NS)
	  END IF
	  K_ON_J(1:NS)=KNU_STORE(1:NS)/JNU_STORE(1:NS)
	END IF
!
	IF(HN_DEF_ON_NODES)THEN
!
! Compute r^2.{H,N} at the mid points. These are then used to compute
! Eddington-like factor at the grid mid points.
!
	  NDM1=ND-1
	  TMP_VEC(1:NS)=R_STORE(1:NS)*R_STORE(1:NS)*HNU_STORE(1:NS)
	  CALL MON_INTERP(RSQ_HNU_MID,NDM1,IONE,RMID,NDM1,TMP_VEC,NS,R_STORE,NS)
	  TMP_VEC(1:NS)=R_STORE(1:NS)*R_STORE(1:NS)*NNU_STORE(1:NS)
	  CALL MON_INTERP(RSQ_NNU_MID,NDM1,IONE,RMID,NDM1,TMP_VEC,NS,R_STORE,NS)
!
	ELSE IF(ND .EQ. NS)THEN

	  DO I=1,NS-1
	    RSQ_HNU_MID(I)=HNU_STORE(I)*RMID_STORE(I)*RMID_STORE(I)
	    RSQ_NNU_MID(I)=NNU_STORE(I)*RMID_STORE(I)*RMID_STORE(I)
	  END DO

	ELSE 			!ND .NE. NS 
!
! Compute H/J at the nodes. We interpolate in r^2 H. Note that H_IN and
! H_OUT are multiplied by r^2 in REGRID_H.
!
	  NDM1=ND-1
	  TMP_VEC(1)=HNU_AT_OB; TMP_VEC(NS+1)=HNU_AT_IB
	  DO I=1,NS-1
	    TMP_VEC(I+1)=HNU_STORE(I)*RMID_STORE(I)*RMID_STORE(I)
	  END DO
	  I=NS+1
	  CALL MON_INTERP(RSQ_HNU_MID,NDM1,IONE,RMID,NDM1,TMP_VEC,I,EXT_RMID_STORE,I)
!
	  TMP_VEC(1)=NNU_AT_OB; TMP_VEC(NS+1)=NNU_AT_IB
	  DO I=1,NS-1
	    TMP_VEC(I+1)=NNU_STORE(I)*RMID_STORE(I)*RMID_STORE(I)
	  END DO
	  I=NS+1
	  CALL MON_INTERP(RSQ_NNU_MID,NDM1,IONE,RMID,NDM1,TMP_VEC,I,EXT_RMID_STORE,I)
!
	END IF
!
! Now compute the "G" Eddington factors. Because H may be zero we have to
! be careful. Depending on N_TYPE we can specify N in terms of J, N in
! terms of H, or N in terms of J and H.
!
	NDM1=ND-1
	IF(N_TYPE .EQ. 'N_ON_J')THEN
	  NMID_ON_HMID(1:NDM1)=0.0D0
	  DO I=1,NDM1
	    NMID_ON_J(I)=RSQ_NNU_MID(I)/(RSQ_JNU(I)+RSQ_JNU(I+1))
	  END DO
	ELSE IF(N_TYPE .EQ. 'MIXED')THEN
	  NMID_ON_J(1:NDM1)=0.0D0
	  DO I=1,NDM1
	    IF(RSQ_HNU_MID(I) .NE. 0)THEN
	      T1=RSQ_NNU_MID(I)/RSQ_HNU_MID(I)
	    ELSE
	      T1=100.0D0
	    END IF
	    IF(T1 .GT. 1.1 .OR. T1 .LT. 0.05)THEN
	      NMID_ON_J(I)=RSQ_NNU_MID(I)/(RSQ_JNU(I)+RSQ_JNU(I+1))
	      NMID_ON_HMID(I)=0.0D0
	    ELSE
	      NMID_ON_HMID(I)=T1
	    END IF
	  END DO
	ELSE IF(N_TYPE .EQ. 'G_ONLY')THEN
!
! Compute G Eddington factor storing in N. We also check the validity of
! the Eddington factor in case strange values are occurring because H
! is near zero (switching sign?).
!
	  NMID_ON_J(1:NDM1)=0.0D0
	  DO I=1,NDM1
	    IF(RSQ_HNU_MID(I) .NE. 0)THEN
	       NMID_ON_HMID(I)=RSQ_NNU_MID(I)/RSQ_HNU_MID(I)
	       IF(NMID_ON_HMID(I) .GT. 1.0D0)THEN
	         NMID_ON_HMID(I)=1.0D0
	       ELSE IF( RSQ_NNU_MID(I) .LT. 0.01) THEN
	         NMID_ON_HMID(I)=0.01D0
	       END IF
	    ELSE
	      NMID_ON_HMID(I)=0.0D0              !Later replaced by average.
	      LUER=ERROR_LU()
	      WRITE(LUER,'(1X,A,1PE16.8)')'HNU zero for frequency:',FREQ
	      WRITE(LUER,'(1X,A,I4)')'Error occurred at depth:',I
	     END IF
	  END DO
	END IF
!
	RETURN
	END
