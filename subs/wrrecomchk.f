C
C Subroutine to write out the recombination and collision rates for
C any ion. The valuse are output to an existing file (overwriting), or
C to a new file.
C
C
	SUBROUTINE WRRECOMCHK(PR,RR,CPR,CRR,DIERECOM,
	1                 ADDRECOM,X_RECOM_1,X_RECOM_2,
	1                 R,T,ED,DHYD,NETRR,TOTRR,N,ND,LU,
	1                 FILNAM,STRDESC)
	USE SET_KIND_MODULE
	IMPLICIT NONE
C
C Altered 26-Jun-1996 - Call to GEN_ASCI_OPEN installed.
C Altered 28-May-1996 : Calls to DP_ZERO removed.
C Altered 22-Jul-1994 - X_RECOM_1 and X_RECOM_2 installed
C Created 15-Feb-1988 - Based on WRPRRRDIE. ADDRECOM now included which
C                         allows for recombinations to levels not explicitly
C                         included.
C Created  5-Oct-1987 (Based on WRPRRRGEN)
C
	INTEGER N,ND,LU,ML,MF,I,MS,J,IOS
	REAL(KIND=LDP) PR(N,ND),RR(N,ND),R(ND),T(ND),ED(ND),DHYD(ND)
	REAL(KIND=LDP) TOTRR(ND),NETRR(ND),CPR(ND),CRR(ND),DIERECOM(ND)
	REAL(KIND=LDP) ADDRECOM(ND),X_RECOM_1(ND),X_RECOM_2(ND)
	CHARACTER*(*) FILNAM,STRDESC
	REAL(KIND=LDP) ABS_SUM
C
	INTEGER ERROR_LU,LUER
	INTEGER, PARAMETER :: IZERO=0
	EXTERNAL ERROR_LU
C
	NETRR(:)=0.0_LDP                 !ND
	TOTRR(:)=0.0_LDP                 !ND
C
	CALL GEN_ASCI_OPEN(LU,FILNAM,'UNKNOWN',' ',' ',IZERO,IOS)
	IF(IOS .NE. 0)THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'Error opening RECOM file',FILNAM
	  WRITE(LUER,*)'IOSTAT=',IOS
	  RETURN
	END IF
C
	MS=1
	DO 10 ML=0,ND-1,10
	  MF=ML+10
	  IF(MF .GT. ND)MF=ND
	  IF(ML .NE. 0)WRITE(LU,'(1H1)')
C
C
	  WRITE(LU,'(/,3X,''Radius [1.0E+10cm] '')')
	  WRITE(LU,999)(R(J),J=MS,MF)
	  WRITE(LU,'(/,3X,''Temperature [1.0E+4K] '')')
	  WRITE(LU,999)(T(J),J=MS,MF)
	  WRITE(LU,'(/,3X,''Electron Density'')')
	  WRITE(LU,999)(ED(J),J=MS,MF)
	  WRITE(LU,'(/,3X,''Ion Density '')')
	  WRITE(LU,999)(DHYD(J),J=MS,MF)
C
	  WRITE(LU,'(/,3X,(A),'' Photoionization Rates'')')STRDESC
	  DO I=1,N
	    WRITE(LU,999)(PR(I,J),J=MS,MF)
	  END DO
C
	  WRITE(LU,'(/3X,''Colisional Ionization Rate '') ')
	  WRITE(LU,999)(CPR(J),J=MS,MF)
C
	  WRITE(LU,'(/,3X,(A),'' Recombination Rates'')')STRDESC
	  DO  I=1,N
	    WRITE(LU,999)(RR(I,J),J=MS,MF)
	  END DO
C
	  WRITE(LU,'(/3X,''Colisional Recombination Rate '') ')
	  WRITE(LU,999)(CRR(J),J=MS,MF)
C
	  IF(DIERECOM(1) .NE. 0 .AND. DIERECOM(ND) .NE. 0)THEN
	    WRITE(LU,'(/3X,''Dielectronic Recombination Rate '') ')
	    WRITE(LU,999)(DIERECOM(J),J=MS,MF)
	  END IF
C
	  IF(ADDRECOM(1) .NE. 0 .AND. ADDRECOM(ND) .NE. 0)THEN
	    WRITE(LU,'(/3X,''Implicit Recombination Rate '') ')
	    WRITE(LU,999)(ADDRECOM(J),J=MS,MF)
	  END IF
C
C Write out net X-ray recombinations from i+1 to i-1.
C (eg If NIV, net recom's from NV to NIII).
C
	  IF(X_RECOM_1(1) .NE. 0 .AND. X_RECOM_1(ND) .NE. 0)THEN
	    WRITE(LU,'(/3X,
	1   ''Net X-ray recombination rate (to previous ionization state)'') ')
	    WRITE(LU,999)(X_RECOM_1(J),J=MS,MF)
	  END IF
C
C Write out net X-ray recombinations from i+2 to i.
C (eg If NIV, net recom's from NVI to NIV).
C
	  IF(X_RECOM_2(1) .NE. 0 .AND. X_RECOM_2(ND) .NE. 0)THEN
	    WRITE(LU,'(/3X,''Net X-ray recombination rate to '',A,
	1	'' g.s.)'') ')STRDESC
	    WRITE(LU,999)(X_RECOM_2(J),J=MS,MF)
	  END IF
C
	  DO J=MS,MF
	    ABS_SUM=0.0_LDP
	    DO I=1,N
	      NETRR(J)=NETRR(J)-PR(I,J)+RR(I,J)
	      TOTRR(J)=TOTRR(J)+RR(I,J)
	      ABS_SUM=ABS_SUM+RR(I,J)+PR(I,J)
	    END DO
	    ABS_SUM=ABS_SUM+CRR(J)+CPR(J)+ABS(DIERECOM(J))+ABS(ADDRECOM(J))+
	1             ABS(X_RECOM_1(J))+ABS(X_RECOM_2(J))
	    NETRR(J)=200.0_LDP*(NETRR(J)+(CRR(J)-CPR(J))+DIERECOM(J)+ADDRECOM(J)
	1                    +X_RECOM_1(J)+X_RECOM_2(J))/ABS_SUM
	    DIERECOM(J)=DIERECOM(J)/(TOTRR(J)+ADDRECOM(J))
	    ADDRECOM(J)=ADDRECOM(J)/ED(J)/DHYD(J)
	    TOTRR(J)=TOTRR(J)/ED(J)/DHYD(J)
	  END DO
C
	  WRITE(LU,'(/3X,''Net Recombination Rate (% of total) '') ')
	  WRITE(LU,999)(NETRR(J),J=MS,MF)
	  WRITE(LU,'(/3X,''Radiative Recombination Coefficient for '',
	1                       ''explicitly treated levels.'') ')
	  WRITE(LU,999)(TOTRR(J),J=MS,MF)
	  IF(ADDRECOM(1) .NE. 0 .AND. ADDRECOM(ND) .NE. 0)THEN
	    WRITE(LU,'(/3X,''Radiative Recombination for'',
	1                         ''implicit levels.'') ')
	    WRITE(LU,999)(ADDRECOM(J),J=MS,MF)
	  END IF
	  IF(DIERECOM(1) .NE. 0 .AND. DIERECOM(ND) .NE. 0)THEN
	    WRITE(LU,'(/3X,''Ratio of Total Dielectronic Recombination'',
	1                       '' to Total Radiatve Recombination'') ')
	    WRITE(LU,999)(DIERECOM(J),J=MS,MF)
	  END IF
C
	  MS=MS+10
10	CONTINUE
C
	CLOSE(UNIT=LU)
	RETURN
999	FORMAT(1X,1P,10E12.4)
	END
