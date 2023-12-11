	SUBROUTINE GET_RSQH_REL(RSQH,R, V, FREQ, ND)
	USE SET_KIND_MODULE
	USE MOD_RAY_MOM_STORE
	IMPLICIT NONE
!
	INTEGER ND
	REAL(KIND=LDP) RSQH(ND)
	REAL(KIND=LDP) R(ND)
	REAL(KIND=LDP) V(ND)
	REAL(KIND=LDP) FREQ
!
! Local vectors
!
	REAL(KIND=LDP) TMP_VEC(1:ND+1)
	REAL(KIND=LDP) RSQ_JNU(1:ND)
	REAL(KIND=LDP) RSQ_HNU_MID(1:ND)
	REAL(KIND=LDP) RSQ_NNU_MID(1:ND)
!
	REAL(KIND=LDP) T1,T2
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
	   ALLOCATE (H_ON_J(ND))
	   ALLOCATE (N_ON_J(ND))
	   ALLOCATE (K_ON_J(ND))
	   ALLOCATE (KMID_ON_J(ND))
	   ALLOCATE (NMID_ON_J(ND))
	   ALLOCATE (NMID_ON_HMID(ND))
	   ALLOCATE (dlnGRSQJdlnR(ND))
!
	   ALLOCATE (H_ON_J_PREV(ND))
	   ALLOCATE (N_ON_J_PREV(ND))
	   ALLOCATE (K_ON_J_PREV(ND))
	   ALLOCATE (KMID_ON_J_PREV(ND))
	   ALLOCATE (NMID_ON_J_PREV(ND))
	   ALLOCATE (NMID_ON_HMID_PREV(ND))
!
	   ALLOCATE (H_ON_J_SAVE(ND))
	   ALLOCATE (N_ON_J_SAVE(ND))
	   ALLOCATE (K_ON_J_SAVE(ND))
	   ALLOCATE (KMID_ON_J_SAVE(ND))
	   ALLOCATE (NMID_ON_J_SAVE(ND))
	   ALLOCATE (NMID_ON_HMID_SAVE(ND))
!
	 END IF
!
	 DO I=1,ND-1
	   RMID(I)=0.5_LDP*(R(I)+R(I+1))
	 END DO
!
! The following behaviour depends on whether H & N are defined at the nodes,
! or at the mid-points.
!
	IF(HN_DEF_ON_NODES)THEN
!
! Compute r^2.{H,N} at the mid points. These are then used to compute
! Eddington-like factor at the grid mid points.
!
	  NDM1=ND-1
	  TMP_VEC(1:NS)=R(1:NS)*R(1:NS)*HNU_STORE(1:NS)
	  CALL MON_INTERP(RSQH,NDM1,IONE,RMID,NDM1,TMP_VEC,NS,R_STORE,NS)
!	
	ELSE
!
	  IF(ND .EQ. NS)THEN
	    DO I=1,NS-1
	      RSQH(I)=RMID_STORE(I)*RMID_STORE(I)*HNU_STORE(I)
	    END DO
	  ELSE
!
! Compute H/J at the mid points in the `new grid.
!
	    I=NS+1
	    TMP_VEC(1)=HNU_AT_OB; TMP_VEC(NS+1)=HNU_AT_IB
	    DO I=1,NS-1
	      TMP_VEC(I+1)=RMID_STORE(I)*RMID_STORE(I)*HNU_STORE(I)
	    END DO
	    CALL MON_INTERP(RSQH,ND,IONE,RMID,ND,TMP_VEC,I,EXT_RMID_STORE,I)
!
	  END IF
	END IF
!
	RETURN
	END
