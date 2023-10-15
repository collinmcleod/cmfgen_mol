	SUBROUTINE WRITE_LINE_OLD(TRAN_NAME,MOD_NAME,
	1            DIF,CONT_INT,FREQ,AMASS,R,V,SIGMA,T,
	1            MASS_DENSITY,CLUMP_FAC,
	1            ETA,CHI,ESEC,CHIL,ETAL,ND,LUOUT,CHI_TH_PASSED)
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
	REAL(10) CONT_INT
	REAL(10) FREQ
	REAL(10) AMASS
!
	INTEGER ND
	INTEGER LUOUT
!
	REAL(10) R(ND)
	REAL(10) V(ND)
	REAL(10) T(ND)
	REAL(10) SIGMA(ND)
	REAL(10) ETA(ND)
	REAL(10) CHI(ND)
	REAL(10) ESEC(ND)
	REAL(10) ETAL(ND)
	REAL(10) CHIL(ND)
	REAL(10) MASS_DENSITY(ND)
	REAL(10) CLUMP_FAC(ND)
	REAL(10) TA(ND)			!Work array
!
	REAL(10) VAC_LAM
	REAL(10) SPEED_OF_LIGHT
	EXTERNAL SPEED_OF_LIGHT
!
	INTEGER I
	CHARACTER(LEN=20) TMP_STR
!
	I=LEN_TRIM(MOD_NAME)
	TMP_STR=MOD_NAME
!
	VAC_LAM=1.0D-07*SPEED_OF_LIGHT()/FREQ
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
