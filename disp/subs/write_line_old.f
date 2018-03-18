	SUBROUTINE WRITE_LINE_OLD(TRAN_NAME,MOD_NAME,
	1            DIF,CONT_INT,FREQ,AMASS,R,V,SIGMA,T,
	1            MASS_DENSITY,CLUMP_FAC,
	1            ETA,CHI,ESEC,CHIL,ETAL,ND,LUOUT)
	IMPLICIT NONE
!
	CHARACTER(LEN=*) TRAN_NAME
	CHARACTER(LEN=*) MOD_NAME
	LOGICAL NEW_FILE
!
	LOGICAL DIF
	REAL*8 CONT_INT
	REAL*8 FREQ
	REAL*8 AMASS
!
	INTEGER ND
	INTEGER LUOUT
!
	REAL*8 R(ND)
	REAL*8 V(ND)
	REAL*8 T(ND)
	REAL*8 SIGMA(ND)
	REAL*8 ETA(ND)
	REAL*8 CHI(ND)
	REAL*8 ESEC(ND)
	REAL*8 ETAL(ND)
	REAL*8 CHIL(ND)
	REAL*8 MASS_DENSITY(ND)
	REAL*8 CLUMP_FAC(ND)
!
	REAL*8 VAC_LAM
	REAL*8 SPEED_OF_LIGHT
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
	WRITE(LUOUT,'(1X,1PE15.8,T25,A,T40,A)')AMASS,'[AMASS]','Atmoic mass (AMUs)'
	WRITE(LUOUT,'(1X,1PE15.8,T25,A,T40,A)')CONT_INT,'[IC]','Schuster intensity'
	WRITE(TMP_STR,'(I5)')ND; TMP_STR=ADJUSTL(TMP_STR)
	WRITE(LUOUT,'(1X,A,T25,A,T40,A)')TRIM(TMP_STR),'[ND]','Number of depth points'
	WRITE(LUOUT,*)' '
	WRITE(LUOUT,
	1       '(1X,T9,A,T23,A,T35,A,T51,A,T64,A,T76,A,T91,A,T105,A,T119,A,T133,A,T147,A)')
	1            'R','T','SIGMA','V','ETA','CHI_TH','ESEC',
	1            'ETAL','CHIL','ROH','f'
        IF(TRAN_NAME.EQ. 'Continuum')THEN
          DO I=1,ND
            CHIL(I)=1.0D-10
            ETAL(I)=1.0D-10
          END DO
       END IF
       DO I=1,ND
         WRITE(LUOUT,'(1X,11ES14.6)')R(I),T(I),SIGMA(I),V(I),
	1                 ETA(I),(CHI(I)-ESEC(I)),
	1                 ESEC(I),ETAL(I),CHIL(I),MASS_DENSITY(I),CLUMP_FAC(I)
	END DO
	CLOSE(UNIT=LUOUT)
!
	RETURN
	END
