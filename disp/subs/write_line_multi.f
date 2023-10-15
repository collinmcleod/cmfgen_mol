	SUBROUTINE WRITE_LINE_MULTI(TRAN_NAME,MOD_NAME,NEW_FILE,
	1            DIF,CONT_INT,FREQ,AMASS,R,V,SIGMA,T,
	1            MASS_DENSITY,CLUMP_FAC,
	1            ETA,CHI,ESEC,CHIL,ETAL,ND,LUOUT)
	IMPLICIT NONE
!
! Altered 27-May-2020 : Now write out N_CONT and N_LINE strings (values will need to be updated).
!
	CHARACTER(LEN=*) TRAN_NAME
	CHARACTER(LEN=*) MOD_NAME
	LOGICAL NEW_FILE
!
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
	IF(NEW_FILE)THEN
	  IF(I .GT. 20)TMP_STR=MOD_NAME(I-19:I)
   	  WRITE(LUOUT,'(1X,A,T25,A,T40,A)')'11-Oct-2017','[Date]','Revised format date'
   	  WRITE(LUOUT,'(1X,A,T25,A,T40,A)')TRIM(TMP_STR),'[MOD_ID]','Model'
 	  WRITE(LUOUT,'(1X,A,T25,A,T40,A)')'TRUE','[DIF]','Diffusion approximation'
 	  WRITE(LUOUT,'(1X,A,T25,A,T40,A)')'TRUE','[FADJ]','Opacities etc scaled by clumping factor'
	  WRITE(TMP_STR,'(I5)')ND; TMP_STR=ADJUSTL(TMP_STR)
	  WRITE(LUOUT,'(1X,A,T25,A,T40,A)')0,'[N_CONT]','Number of continuum frequencies'
	  WRITE(LUOUT,'(1X,A,T25,A,T40,A)')0,'[N_LINE]','Number of line frequencies'
	  WRITE(LUOUT,'(1X,A,T25,A,T40,A)')TRIM(TMP_STR),'[ND]','Number of depth points'
	  WRITE(LUOUT,*)' '
!
	  WRITE(LUOUT,'(1X,A,/,(1P,1X,9E14.6))')'R',(R(I),I=1,ND)
	  WRITE(LUOUT,'(1X,A,/,(1P,1X,9E14.6))')'V',(V(I),I=1,ND)
	  WRITE(LUOUT,'(1X,A,/,(1P,1X,9E14.6))')'SIGMA',(SIGMA(I),I=1,ND)
	  WRITE(LUOUT,'(1X,A,/,(1P,1X,9E14.6))')'T',(T(I),I=1,ND)
	  WRITE(LUOUT,'(1X,A,/,(1P,1X,9E14.6))')'Mass density',(MASS_DENSITY(I),I=1,ND)
	  WRITE(LUOUT,'(1X,A,/,(1P,1X,9E14.6))')'Clumping factor',(CLUMP_FAC(I),I=1,ND)
	  WRITE(LUOUT,'(1X,A,/,(1P,1X,9E14.6))')'ESEC',(ESEC(I),I=1,ND)
	END IF
!
	VAC_LAM=1.0D-07*SPEED_OF_LIGHT()/FREQ
	IF(TRAN_NAME .EQ. 'Continuum')THEN
	  WRITE(LUOUT,*)' '
	  WRITE(LUOUT,'(1X,A,T25,A,T40,A)')TRIM(TRAN_NAME),'[TR_ID]','Transition identification'
	  WRITE(LUOUT,'(1X,1PE15.8,T25,A,T40,A)')CONT_INT,'[IC]','Schuster intensity'
	  WRITE(LUOUT,'(1X,1PE15.8,T25,A,T40,A)')FREQ,'[FREQ]','Frequency (10^15 Hz)'
	  WRITE(LUOUT,'(1X,1PE12.5,T25,A,T40,A)')VAC_LAM,'[LAM]','Vacuum wavelength (Ang)'
	  WRITE(LUOUT,'(1X,A,/,(1P,1X,9E14.6))')'ETA',(ETA(I),I=1,ND)
	  WRITE(LUOUT,'(1X,A,/,(1P,1X,9E14.6))')'CHI_TH',( (CHI(I)-ESEC(I)) ,I=1,ND)
	ELSE
	  WRITE(LUOUT,*)' '
	  WRITE(LUOUT,'(1X,A,T25,A,T40,A)')TRIM(TRAN_NAME),'[TR_ID]','Transition identification'
	  WRITE(LUOUT,'(1X,1PE15.8,T25,A,T40,A)')CONT_INT,'[IC]','Schuster intensity'
	  WRITE(LUOUT,'(1X,1PE15.8,T25,A,T40,A)')AMASS,'[AMASS]','Atomic mass (AMUs)'
	  WRITE(LUOUT,'(1X,1PE15.8,T25,A,T40,A)')FREQ,'[FREQ]','Frequency (10^15 Hz)'
	  WRITE(LUOUT,'(1X,1PE12.5,T25,A,T40,A)')VAC_LAM,'[LAM]','Vacuum wavelength (Ang)'
	  WRITE(LUOUT,'(1X,A,/,(1P,1X,9E14.6))')'ETAL',(ETAL(I),I=1,ND)
	  WRITE(LUOUT,'(1X,A,/,(1P,1X,9E14.6))')'CHIL',(CHIL(I),I=1,ND)
	END IF
	CLOSE(UNIT=LUOUT)
!
	RETURN
	END
