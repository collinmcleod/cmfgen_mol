	SUBROUTINE WRITE_LINE_OLD(TRAN_NAME,MOD_NAME,
	1            DIF,CONT_INT,FREQ,AMASS,R,V,SIGMA,T,
	1            MASS_DENSITY,CLUMP_FAC,
	1            ETA,CHI,ESEC,CHIL,ETAL,ND,LUOUT,CHI_TH_PASSED)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Altered: 21-Feb-2020   Added optional variable CHI_TH_PASSED
! Altered: 12-Dec-2018   CHIL, ETAL no longer altered in rountine if TRAN_NAME set to continuum.
!
	CHARACTER(LEN=*) TRAN_NAME
	CHARACTER(LEN=*) MOD_NAME
	LOGICAL NEW_FILE
!
! CHI_TH_PASSED indicates whether CHI (total op.)  or just the thermal opacity has
! been passed in the variable CHI.
!
	LOGICAL, OPTIONAL :: CHI_TH_PASSED
	LOGICAL DIF
	REAL(KIND=LDP) CONT_INT
	REAL(KIND=LDP) FREQ
	REAL(KIND=LDP) AMASS
!
	INTEGER ND
	INTEGER LUOUT
!
	REAL(KIND=LDP) R(ND)
	REAL(KIND=LDP) V(ND)
	REAL(KIND=LDP) T(ND)
	REAL(KIND=LDP) SIGMA(ND)
	REAL(KIND=LDP) ETA(ND)
	REAL(KIND=LDP) CHI(ND)
	REAL(KIND=LDP) ESEC(ND)
	REAL(KIND=LDP) ETAL(ND)
	REAL(KIND=LDP) CHIL(ND)
	REAL(KIND=LDP) MASS_DENSITY(ND)
	REAL(KIND=LDP) CLUMP_FAC(ND)
	REAL(KIND=LDP) TA(ND)			!Work array
!
	REAL(KIND=LDP) VAC_LAM
	REAL(KIND=LDP) SPEED_OF_LIGHT
	EXTERNAL SPEED_OF_LIGHT
!
	INTEGER I
	CHARACTER(LEN=20) TMP_STR
!
	I=LEN_TRIM(MOD_NAME)
	TMP_STR=MOD_NAME
!
	VAC_LAM=1.0E-07_LDP*SPEED_OF_LIGHT()/FREQ
	IF(I .GT. 20)TMP_STR=MOD_NAME(I-19:I)
   	WRITE(LUOUT,'(1X,A,T25,A,T40,A)')'16-Sep-2017','[Date]','Revised format date'
   	WRITE(LUOUT,'(1X,A,T25,A,T40,A)')TRIM(TMP_STR),'[MOD_ID]','Model'
   	WRITE(LUOUT,'(1X,A,T25,A,T40,A)')TRIM(TRAN_NAME),'[TR_ID]','Transition identification'
 	WRITE(LUOUT,'(1X,A,T25,A,T40,A)')'TRUE','[DIF]','Diffusion approximation'
	WRITE(LUOUT,'(1X,1PE15.8,T25,A,T40,A)')FREQ,'[FREQ]','Frequency (10^15 Hz)'
	WRITE(LUOUT,'(1X,1PE12.5,T25,A,T40,A)')VAC_LAM,'[LAM]','Vacuum wavelength (Ang)'
	WRITE(LUOUT,'(1X,1PE15.8,T25,A,T40,A)')AMASS,'[AMASS]','Atomic mass (AMUs)'
	WRITE(LUOUT,'(1X,1PE15.8,T25,A,T40,A)')CONT_INT,'[IC]','Schuster intensity'
	WRITE(TMP_STR,'(I5)')ND; TMP_STR=ADJUSTL(TMP_STR)
	WRITE(LUOUT,'(1X,A,T25,A,T40,A)')TRIM(TMP_STR),'[ND]','Number of depth points'
	WRITE(LUOUT,*)' '
	WRITE(LUOUT,
	1       '(1X,T9,A,T23,A,T35,A,T51,A,T64,A,T76,A,T91,A,T105,A,T119,A,T133,A,T147,A)')
	1            'R','T','SIGMA','V','ETA','CHI_TH','ESEC',
	1            'ETAL','CHIL','ROH','f'
!
	TA=CHI-ESEC
	IF(PRESENT(CHI_TH_PASSED))THEN
	  IF(CHI_TH_PASSED)TA=CHI
	END IF
!
       IF(TRAN_NAME.EQ. 'Continuum')THEN
         DO I=1,ND
           WRITE(LUOUT,'(1X,11ES14.6)')R(I),T(I),SIGMA(I),V(I),
	1                 ETA(I),TA(I),
	1                 ESEC(I),1.0D-10,1.0D-10,MASS_DENSITY(I),CLUMP_FAC(I)
	  END DO
	ELSE
          DO I=1,ND
            WRITE(LUOUT,'(1X,11ES14.6)')R(I),T(I),SIGMA(I),V(I),
	1                 ETA(I),TA(I),
	1                 ESEC(I),ETAL(I),CHIL(I),MASS_DENSITY(I),CLUMP_FAC(I)
	  END DO
	END IF
	CLOSE(UNIT=LUOUT)
!
	RETURN
	END
