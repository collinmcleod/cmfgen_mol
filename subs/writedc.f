	SUBROUTINE WRITEDC(HYD,HYDLTE,NHYD,DHYD,NION,R,T,ED,V,
	1               LUM,ND,FILENAME,OPTION,FORM)
	USE SET_KIND_MODULE
	IMPLICIT NONE
C
C Altered 26-Jun-1996 : CALL GEN_ASCI_OPEN installed.
C Altered 28-May-1996 : Removed for [jdh.disp]SETVEC routine
C                       DOUBLE PRECISION declaration removed.
C
C Altered  4-Aug-1988 : Write Departure coefficients out - not b-1.
C
	INTEGER NHYD,NION,ND,FORM,I,J,IOS
	INTEGER, PARAMETER :: IZERO=0
	REAL(KIND=LDP) HYD(NHYD,ND),HYDLTE(NHYD,ND),DHYD(NION,ND)
	REAL(KIND=LDP) R(ND),T(ND),ED(ND),V(ND),LUM,T1,T2
	CHARACTER*(*)FILENAME,OPTION
	CHARACTER*30 NEWNAME
	CHARACTER*90 FMT
C
	INTEGER ERROR_LU,LUER
	EXTERNAL ERROR_LU
C
C	NEWNAME=FILENAME//'.'//OPTION
	NEWNAME=FILENAME
C
	  IF(DHYD(1,ND) .NE. 0)THEN
C
C 1 = H, HeII
C 2 = HeI Sing
C 3 = HeI triplets, CIV,NV
C 4 = CIII,NIV
C 5 = HeI (Singlets and Triplets)
C
	    IF(FORM .EQ. 1)FMT='(1X,1P5E15.5)'
	    IF(FORM .EQ. 2)FMT='(1X,1P1E15.5,:/X,2E15.5,:/X,3E15.5,
	1                        :/,(X4E15.5_LDP))'
	    IF(FORM .EQ. 3)FMT='(1X,1P2E15.5,:/X,3E15.5,:/,(X4E15.5))'
	    IF(FORM .EQ. 4)FMT='(1X,1P1E15.5,:/X,2E15.5,:/X,3E15.5,
	1                        :/,(X6E15.5_LDP))'
	    IF(FORM .EQ. 5)FMT='(1X,1P1E15.5,:/X,2E15.5,:/X,2E15.5,:/X,
	1            3E15.5_LDP,:/X,3E15.5_LDP,/:,(X,6E15.5_LDP))'
C
	    I=9
	    CALL GEN_ASCI_OPEN(I,NEWNAME,'REPLACE',' ',' ',IZERO,IOS)
	    IF(IOS .NE. 0)THEN
	      LUER=ERROR_LU()
	      WRITE(LUER,*)'Error opening D.C, file',NEWNAME
	      WRITE(LUER,*)'IOSTAT=',IOS
	      RETURN
	    END IF
	    WRITE(9,2120)R(ND),LUM,NHYD,ND
	    IF(OPTION .EQ. 'DC')THEN
	      DO I=1,ND
	        T1=0.0_LDP
	        T2=0.0_LDP
	        DO J=1,NHYD
	          T1=T1+HYD(J,I)
	        END DO
	        DO J=1,NION
	          T2=T2+DHYD(J,I)
	        END DO
	        T1=T1/T2
	        WRITE(9,2122)R(I),DHYD(1,I),ED(I),T(I),T1,V(I)
	        WRITE(9,FMT)((HYD(J,I)/HYDLTE(J,I)),J=1,NHYD)
	      END DO
	    ELSE
	      DO I=1,ND
	        T1=0.0_LDP
	        T2=0.0_LDP
	        DO J=1,NHYD
	          T1=T1+HYD(J,I)
	        END DO
	        DO J=1,NION
	          T2=T2+DHYD(J,I)
	        END DO
	        T1=T1/T2
	        WRITE(9,2122)R(I),DHYD(1,I),ED(I),T(I),T1,V(I)
	        WRITE(9,FMT)(HYD(J,I),J=1,NHYD)
	      END DO
	    END IF
	    CLOSE(UNIT=9)
	  END IF
C
2120	  FORMAT(/,1X,F9.4,4X,1PE11.4,5X,0P,I4,5X,I4)
2122	  FORMAT(/,1X,1P6E15.5)
C
	  RETURN
	  END
