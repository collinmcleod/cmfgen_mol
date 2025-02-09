C
	SUBROUTINE NEW_WRITEDC_V4(HYD,HYDLTE,WHYD,
	1               EDGEHYD,GHYD,NHYD,
	1               DHYD,GION,NION,R,T,ED,V,CLUMP_FAC,
	1               DO_DPTH,LUM,ND,FILENAME,OPTION,FORM)
	IMPLICIT NONE
C
C 19-Jun-2000 - DO_DPTH inserted (changed to V4)
C               TEXCITE calculation improved.
C 05-Mar-1999 - TEXCITE now implicitly dimensioned by NHYD, not 200
C               WHYD (Level dissolution coeffiecients) inclued in call for
C                 forrect evaluation of TX.
C 07-Jul-1997 - CLUMP_FAC installed in call, and called _V2.
C                 R now written with 7 digits of precision.
C 17-Nov-1989 -Based on WRITEDC : TX option installed.
C
	INTEGER*4 NHYD,NION,ND,FORM
	REAL*8 HYD(NHYD,ND),HYDLTE(NHYD,ND),WHYD(NHYD,ND)
	REAL*8 EDGEHYD(NHYD),GHYD(NHYD)
	REAL*8 DHYD(NION,ND),R(ND),T(ND),ED(ND),V(ND),CLUMP_FAC(ND)
	REAL*8 LUM
	REAL*8 GION
!
! DPTH_IND indicates which depths are to be output. Useful for
! debugging purposes when running a new model.
!
	LOGICAL DO_DPTH(ND)
	CHARACTER*(*)FILENAME,OPTION
	CHARACTER*90 FMT
!
	REAL*8 T1,T2,DELTA_T
	REAL*8 TEXCITE(NHYD)
	INTEGER*4 I,J,COUNT,ND_CNT
!
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/LINE/ OPLIN,EMLIN
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ,OPLIN,EMLIN
C
	  IF(DHYD(1,ND) .NE. 0)THEN
C
C 1 = H, HeII
C 2 = HeI Sing
C 3 = HeI triplets, CIV,NV
C 4 = CIII,NIV
C 5 = HeI (Singlets and Triplets)
C
	    IF(FORM .EQ. 1)FMT='(X,1P5E15.5)'
	    IF(FORM .EQ. 2)FMT='(X,1P1E15.5,:/X,2E15.5,:/X,3E15.5,
	1                        :/,(X4E15.5))'
	    IF(FORM .EQ. 3)FMT='(X,1P2E15.5,:/X,3E15.5,:/,(X4E15.5))'
	    IF(FORM .EQ. 4)FMT='(X,1P1E15.5,:/X,2E15.5,:/X,3E15.5,
	1                        :/,(X6E15.5))'
	    IF(FORM .EQ. 5)FMT='(X,1P1E15.5,:/X,2E15.5,:/X,2E15.5,:/X,
	1            3E15.5,:/X,3E15.5,/:,(X,6E15.5))'
C
	    OPEN(UNIT=9,STATUS='NEW',FILE=FILENAME)
!
! Determine number points that wll be output.
!
	    ND_CNT=0
	    DO I=1,ND
	      IF(DO_DPTH(I))ND_CNT=ND_CNT+1
	    END DO
!
	    WRITE(9,'(/,X,A,T40,A)')'07-Jul-1997','!Format date'
	    WRITE(9,2120)R(ND),LUM,NHYD,ND_CNT
	    IF(OPTION .EQ. 'DC')THEN
	      DO I=1,ND
	        IF(DO_DPTH(I))THEN
	          T1=0.0D0
	          T2=0.0D0
	          DO J=1,NHYD
	            T1=T1+HYD(J,I)
	          END DO
	          DO J=1,NION
	            T2=T2+DHYD(J,I)
	          END DO
	          T1=T1/T2
	          WRITE(9,2122)R(I),DHYD(1,I),ED(I),T(I),T1,V(I),CLUMP_FAC(I)
	          WRITE(9,FMT)((HYD(J,I)/HYDLTE(J,I)),J=1,NHYD)
	        END IF
	      END DO
	    ELSE IF(OPTION .EQ. 'TX')THEN
	      DO I=1,ND_CNT
	        IF(DO_DPTH(I))THEN
	          T1=0.0D0
	          T2=0.0D0
	          DO J=1,NHYD
	            T1=T1+HYD(J,I)
	          END DO
	          DO J=1,NION
	            T2=T2+DHYD(J,I)
	          END DO
	          T1=T1/T2
	          WRITE(9,2122)R(I),DHYD(1,I),ED(I),T(I),T1,V(I),CLUMP_FAC(I)
	          DO J=1,NHYD
	            T1=LOG( HYD(J,I)*GION/GHYD(J)/WHYD(J,I)/
	1                             2.07D-22/ED(I)/DHYD(1,I) )
	            T2=HDKT*EDGEHYD(J)
	            DELTA_T=10.0D0
	            COUNT=0
	            TEXCITE(J)=T(I)
	            DO WHILE( ABS(DELTA_T) .GT. 1.0E-06 .AND. COUNT .LT. 100 )
	              COUNT=COUNT+1
	              DELTA_T=( T1- T2/TEXCITE(J) + 1.5D0*LOG(TEXCITE(J)) )*
	1                     TEXCITE(J)/(T2/TEXCITE(J)+1.5D0)
	              IF(DELTA_T .GT. 0.8*TEXCITE(J))DELTA_T=0.8*TEXCITE(J)
	              IF(DELTA_T .LT. -0.8*TEXCITE(J))DELTA_T=-0.8*TEXCITE(J)
	              TEXCITE(J)=TEXCITE(J)-DELTA_T
	            END DO
	            IF(COUNT .EQ. 100)THEN
	              WRITE(6,*)'Error - TEXC didnt converge in 100 iterations'
	              WRITE(6,*)'I,J=',I,J
	              RETURN
	            END IF
	          END DO
	          WRITE(9,FMT)(TEXCITE(J),J=1,NHYD)
	        END IF
	      END DO
	    ELSE
	      DO I=1,ND_CNT
	        IF(DO_DPTH(I))THEN
	          T1=0.0D0
	          T2=0.0D0
	          DO J=1,NHYD
	            T1=T1+HYD(J,I)
	          END DO
	          DO J=1,NION
	            T2=T2+DHYD(J,I)
	          END DO
	          T1=T1/T2
	          WRITE(9,2122)R(I),DHYD(1,I),ED(I),T(I),T1,V(I),CLUMP_FAC(I)
	          WRITE(9,FMT)(HYD(J,I),J=1,NHYD)
	        END IF
	      END DO
	    END IF
	    CLOSE(UNIT=9)
	  END IF
C
2120	  FORMAT(/,X,F9.4,5X,1PE10.4,5X,0P,I4,5X,I4)
2122	  FORMAT(/,X,1P,E15.7,6E15.5)
C
	  RETURN
	  END
