	SUBROUTINE DO_FULL_R_GRID_V1(R_OLD,DONE_R_REV,LU,ND)
	USE SET_KIND_MODULE
	USE MOD_CMFGEN
	IMPLICIT NONE
!
	INTEGER ND
	INTEGER LU
	REAL(KIND=LDP) R_OLD(ND)
	LOGICAL DONE_R_REV
!
	REAL(KIND=LDP) DTAU(ND)
	REAL(KIND=LDP) OLD_TAU(ND)
	REAL(KIND=LDP) TCHI(ND)
	REAL(KIND=LDP) dCHIdR(ND)
!
	REAL(KIND=LDP) T1
	REAL(KIND=LDP) R_SCALE_FAC
	REAL(KIND=LDP) TAU_SCALE_FAC
	REAL(KIND=LDP) IB_RAT
	REAL(KIND=LDP) OB_RAT
	REAL(KIND=LDP) DTAU2_ON_DTAU1
	REAL(KIND=LDP) dLOGT_MAX
!
	INTEGER I
	INTEGER NIB
	INTEGER NOB
	INTEGER NEW_ND
!
	CHARACTER(LEN=10) METHOD
	CHARACTER(LEN=10) OPAC_MEAN
!
	LOGICAL, PARAMETER :: L_FALSE=.FALSE.
!
! Set reasonable defaults.
!
	NIB=2; NOB=1
	NEW_ND=ND
	R_SCALE_FAC=1.4
	IB_RAT=2.0D0; OB_RAT=1.5D0
	TAU_SCALE_FAC=0.0D0
	DTAU2_ON_DTAU1=100.0D0
	dLOGT_MAX=0.04D0
	OPAC_MEAN='ROSS'
	METHOD='LOGLOG'
	DONE_R_REV=.FALSE.
!
	WRITE(6,'(A)')' '
	CALL RD_STORE_DBLE(R_SCALE_FAC,'R_SCL_FAC',L_FALSE,'Factor (>1) to enhance maximum dLog(R) spacing')
	CALL RD_STORE_DBLE(dLOGT_MAX,'dLOGT_MAX',L_FALSE,'Maximum fractional change in the temperature')
!
	CALL RD_STORE_INT(NIB,'INS_N_IB',L_FALSE,'Number of depth points to insert at INNER boundary')
	CALL RD_STORE_DBLE(IB_RAT,'IB_RAT',L_FALSE,'Ratio in optical depth increments at INNER boundry (>1)')
!
	CALL RD_STORE_INT(NOB,'INS_N_OB',L_FALSE,'Number of depth points to insert at OUTER boundary')
	CALL RD_STORE_DBLE(OB_RAT,'OB_RAT',L_FALSE,'Ratio in optical depth increments at OUTER boundry (>1)')
!
	CALL RD_STORE_CHAR(OPAC_MEAN,'OPAC_MEAN',L_FALSE,'Which opacity to use for optical depth sclae')
	CALL RD_STORE_DBLE(TAU_SCALE_FAC,'TAU_SCL_FAC',L_FALSE,'Factor to reduce optical depth: TAU=TAU-SF*TAU(1)')
	CALL RD_STORE_DBLE(DTAU2_ON_DTAU1,'D2OND1',L_FALSE,'~DTAU(2)/DTAU(1) at outer boudary')
!
	IF(OPAC_MEAN(1:4) .EQ. 'ROSS')THEN
	  TCHI(1:ND)=ROSS_MEAN(1:ND)*CLUMP_FAC(1:ND)
	  CALL DERIVCHI(dCHIdR,TCHI,R,ND,METHOD)
	  CALL NORDTAU(DTAU,TCHI,R,R,dCHIdR,ND)
	ELSE IF(OPAC_MEAN(1:4) .EQ. 'FLUX')THEN
	  DO I=1,ND
	    T1=6.65D-15*ED(I)
	    TCHI(I)=MAX(FLUX_MEAN(I),T1)*CLUMP_FAC(I)
	  END DO
	  CALL DERIVCHI(dCHIdR,TCHI,R,ND,METHOD)
	  CALL NORDTAU(DTAU,TCHI,R,R,dCHIdR,ND)
	ELSE IF(OPAC_MEAN(1:4) .EQ. 'SQRT')THEN
	  DO I=1,ND
	    TCHI(I)=SQRT(MAX(FLUX_MEAN(I),ROSS_MEAN(I))*ROSS_MEAN(I))*CLUMP_FAC(I)
	  END DO
	  CALL DERIVCHI(dCHIdR,TCHI,R,ND,METHOD)
	  CALL NORDTAU(DTAU,TCHI,R,R,dCHIdR,ND)
	END IF
!
	T1=LOG(DENSITY(5)/DENSITY(1))/LOG(R(1)/R(5))
	T1=TCHI(1)*R(1)/MAX(2.0D0,T1-1.0D0)*(1.0D0-TAU_SCALE_FAC)
	OLD_TAU(1)=T1
	DO I=2,ND
	  OLD_TAU(I)=OLD_TAU(I-1)+DTAU(I-1)
	END DO
	NEW_ND=ND
	CALL ADJUST_SN_R_GRID(R,R_OLD,T,OLD_TAU,R_SCALE_FAC,dLOGT_MAX,
	1         IB_RAT,OB_RAT,DTAU2_ON_DTAU1,NIB,NOB,NEW_ND,ND)
!
	DONE_R_REV=.TRUE.
!
	RETURN
	END
