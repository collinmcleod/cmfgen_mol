C
C***************************************************************
C
	FUNCTION VOIGT(A,V)
C
C Subroutine to evaluate a VOIGT profile, normalized so that the integral
C over V is unity.
C
C Routine is designed for speed.
C
C Adapted from a Midas which was based on a routined obtained from
C  S L Wright /UCL/ and a routine from Peter Hofflich which inturn
C was based on HUI,A.K.,ARMSTRONG,E.H.,WRAY,A.A. 1978,JQSRT 19,P.509
C
C Accuracy:
C         < 0.1% for 0.1 < a 
C         < 0.5% for a < 0.1
C         < 0.3% for a < 0.05
C
	IMPLICIT NONE
C                                                                   
	REAL*8  VOIGT,V,A
C
C Variables for interpolation (a < 0.1).
C
	INTEGER N,N1
	REAL*8  V0,V2
	REAL*8  H1(81)
C
C variables for functional form (a > 0.1)
C
	COMPLEX ZH,F
	REAL*8 C(7),D(7)
C
C Interpolation constants.
C
	DATA H1/
	1  -1.1283800,-1.1059600,-1.0404800,-0.9370300,-0.8034600,-0.6494500,
	1  -0.4855200,-0.3219200,-0.1677200,-0.0301200, 0.0859400, 0.1778900,
	1   0.2453700, 0.2898100, 0.3139400, 0.3213000, 0.3157300, 0.3009400,
	1   0.2802700, 0.2564800, 0.2317260, 0.2075280, 0.1848820, 0.1643410,
	1   0.1461280, 0.1302360, 0.1165150, 0.1047390, 0.0946530, 0.0860050,
	1   0.0785650, 0.0721290, 0.0665260, 0.0616150, 0.0572810, 0.0534300,
	1   0.0499880, 0.0468940, 0.0440980, 0.0415610, 0.0392500, 0.0351950,
	1   0.0317620, 0.0288240, 0.0262880, 0.0240810, 0.0221460, 0.0204410,
	1   0.0189290, 0.0175820, 0.0163750, 0.0152910, 0.0143120, 0.0134260,
	1   0.0126200, 0.0118860, 0.0112145, 0.0105990, 0.0100332, 0.0095119,
	1   0.0090306, 0.0085852, 0.0081722, 0.0077885, 0.0074314, 0.0070985,
	1   0.0067875, 0.0064967, 0.0062243, 0.0059688, 0.0057287, 0.0055030,
	1   0.0052903, 0.0050898, 0.0049006, 0.0047217, 0.0045526, 0.0043924,
	1   0.0042405, 0.0040964, 0.0039595/
C
C Functional constants.
C
	DATA C/122.6079,214.3823,181.9285,93.15558,30.18014,5.912626,
	1     5.641896E-01/
	DATA D/122.6079,352.7306,457.3344,348.7039,170.3540,53.99291,
	1     10.47986/
C
	IF(A .GT. 0.1)THEN
          ZH = CMPLX(A,(-ABS(V)))
          F = ((((((C(7)*ZH+C(6))*ZH+C(5))*ZH+C(4))*ZH+C(3))*ZH+C(2))
	1      *ZH+C(1))
	1      /(((((((ZH+D(7))*ZH+D(6))*ZH+D(5))*ZH+D(4))*ZH+D(3))*ZH+D(2)
	1      )*ZH+D(1))
	  VOIGT = REAL(F)
	ELSE IF(A .EQ. 0)THEN
	  VOIGT=EXP(-V*V)
	ELSE
C
C In original UCL routine, H0 and H2 were interpolated. Better to compute ---
C not much slower and improves accuracy from 3% to beter than 0.5% for
C a=0.001.
C
	  V0=ABS(V)*10.D0
	  N=V0
          IF(N.LT.40)THEN
	    V2=V0-N                     
	    N=N+1
	    N1=N+1
	    VOIGT=(1.0D0+A*A*(1-2*V*V))*EXP(-V*V) + 
	1            A*(H1(N)+V2*(H1(N1)-H1(N)) )
	  ELSE IF(N.LT.120)THEN
	    N=N/2+20
     	    V2=20.0D0+0.5D0*V0-N
	    N=N+1
	    N1=N+1
	    VOIGT=A*((H1(N1)-H1(N))*V2+H1(N))
	    IF(V0 .LT. 50)VOIGT=VOIGT+EXP(-V*V)
	  ELSE
            VOIGT=(0.56419D0+0.846D0/(V*V))/(V*V)*A
	  END IF
	END IF
	VOIGT=VOIGT/1.772453851	  	!SQRT(PI)
C
	RETURN
	END
