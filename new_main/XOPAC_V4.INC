!
! Compute opacities and emissivities due to K (& L) shell ionization. In all cases
! it is assumed that 2 electrons are ejected. The K (& L) shell cross-sections are
! assumed to be independent of the level of the valence electron. In practice,
! ionizations will generally be determined by the population of the ground
! configuration.
!
! NB: This differs from XOPAC.RAW used by CMFGEN in that the full
! arrays, not the super level atoms, are accessed.
!
! Since the cross-sections are level independent, we can sum over the levels before
! we add the contribution to the opacity/emissivity.
!
! Altered 12-Oct-2003 : We now explicitly check if the XRAY cross-section is zero.
!                         If zero, we can skip calculations. Avoids floating overflow 
!                         problems at loww temperatures, and is computationally more efficient.
! Altered: 19/May/2002  LST_DEPTH_ONLY IF installed to save time in DTDR computation.
!                       Variable T4 introduced to save computational effort.
! Altered: Bux fix. NCI replaced by NXzV. Provided CI was included in
!                     Compile code, minimal effect since opacity dominated by
!                     lowest levels.
!
	IF(LST_DEPTH_ONLY)THEN
	  DO ID=1,NUM_IONS
	    IF(ATM(ID)%XzV_PRES .AND. ATM(ID+1)%XzV_PRES)THEN
	      T2=AT_NO(SPECIES_LNK(ID))+1-ATM(ID)%ZXzV
	      T1=XCROSS_V2(FL,AT_NO(SPECIES_LNK(ID)),T2,IZERO,IZERO,L_FALSE,L_FALSE)
	      IF(T1 .NE. 0.0D0)THEN
	        I=ND
	        T2=0.0D0			!Temporary CHI
	        T3=0.0D0			!Temporary ETA
	        T4=ATM(ID+1)%XzVLTE_F(1,I)/ATM(ID+1)%XzV_F(1,I)*EMHNUKT(I)
	        DO J=1,ATM(ID)%NXzV_F
		  T2=T2+ATM(ID)%XzV_F(J,I)
	          T3=T3+ATM(ID)%XzVLTE_F(J,I)
	        END DO
	        CHI(I)=CHI(I)+T1*(T2-T3*T4)
	        ETA(I)=ETA(I)+T1*T3*T4*TWOHCSQ*(FL**3)
	      END IF
	    END IF
	  END DO
	ELSE
	  DO ID=1,NUM_IONS
	    IF(ATM(ID)%XzV_PRES .AND. ATM(ID+1)%XzV_PRES)THEN
	      T2=AT_NO(SPECIES_LNK(ID))+1-ATM(ID)%ZXzV
	      T1=XCROSS_V2(FL,AT_NO(SPECIES_LNK(ID)),T2,IZERO,IZERO,L_FALSE,L_FALSE)
	      IF(T1 .NE. 0.0D0)THEN
	        DO I=1,ND
	          T2=0.0D0			!Temporary CHI
	          T3=0.0D0			!Temporary ETA
	          T4=ATM(ID+1)%XzVLTE_F(1,I)/ATM(ID+1)%XzV_F(1,I)*EMHNUKT(I)
	          DO J=1,ATM(ID)%NXzV_F
		    T2=T2+ATM(ID)%XzV_F(J,I)
	            T3=T3+ATM(ID)%XzVLTE_F(J,I)
	          END DO
	          CHI(I)=CHI(I)+T1*(T2-T3*T4)
	          ETA(I)=ETA(I)+T1*T3*T4*TWOHCSQ*(FL**3)
	        END DO
	      END IF
	    END IF
	  END DO
	END IF 
C
