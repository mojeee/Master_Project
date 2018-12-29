c %%%%%%%%%%%%%%%%%%%%%%
c
c GETS GRADIENTS OF C
c
c %%%%%%%%%%%%%%%%%%%%%%
c
      PROGRAM GRADC
      IMPLICIT NONE

      INTEGER N, NBIN, NSIM, NC, I, II, III, SDRCOUNT(100)
      PARAMETER(N = 2122, NBIN = 50, NSIM = 3500)
      DOUBLE PRECISION DC, C_VEC(N,NSIM), GRADC_VEC((N-2),2),
     1                 DX,TDX,X_VEC(N),
     2                 BINS(NBIN),BINSB(NBIN+1),
     3                 TEMP,
     4                 SDRVECTOR(NBIN)
      LOGICAL SORTED(N),VECTOR_ENTRY(NBIN)
CCC SET BINS AND BINSB(ORDER) VECTORS
      DC = 1.D0/( DFLOAT(NBIN) )
      DO I = 1, NBIN
         BINSB(I) = (I-1)*DC
      ENDDO
      BINSB(NBIN+1) = 1.2D0!purposely make last bin bigger 
      DO I = 1, (NBIN-1)
         BINS(I) = (BINSB(I+1)+BINSB(I))/2D0
      ENDDO
      BINS(NBIN) = (BINSB(NBIN)+1D0)/2D0
CCC OPEN INPUT FILES 
      OPEN(1, file="x.csv" , status="unknown")!distance data 
      OPEN(2, file="interpedc.csv", status="unknown")!lem data
CCC CREATE OUTPUT FILES
      OPEN(3, FILE="gradc.csv", status="unknown")
CCC READ LEM DATA INTO ARRAYS
      DO II = 1, N
         READ(1,*) X_VEC(II)
      ENDDO
      DO I = 1, NSIM
         DO II = 1, N
            READ(2,*) C_VEC(II,I)
         ENDDO
      ENDDO  
CCC CHECK THAT DX IS CONSTANT AND GET DX
      DX = X_VEC(10) - X_VEC(9)
      TEMP = X_VEC(100) - X_VEC(99)
      IF( (DX - TEMP) .GT. 1D-7 ) THEN
         PRINT *, 'CHECK DX: DX IS NOT CONSTANT'
         STOP
      ELSE
         DX = X_VEC(2) - X_VEC(1)
         TDX = 2D0*DX
      ENDIF    
CCC INITIALIZE SDRCOUNT VECTOR AS ZEROS
ccc (counts number of SDRs included for each bin)
      DO I = 1, NBIN
         SDRCOUNT(I) = 0
      ENDDO
CCC SET THE DEFAULT SDRVEDCTOR VALUES TO 0.D0
      DO I = 1, NBIN
         SDRVECTOR(I) = 0.D0
      ENDDO
C
C
CCC START CALCULATING GRADIENTS 
C
C
      DO I = 1, NSIM
CCC CALC GRADIENTS
         CALL NORMDATA(N,I,C_VEC)
         CALL GET_GRADIENTS(N,I,TDX,C_VEC,GRADC_VEC)
         CALL SORT(N,NBIN,BINSB,GRADC_VEC,SDRVECTOR,SDRCOUNT)
      ENDDO
CCC AVERAGE PDF ACROSS REALIZATIONS AND WRITE TO FILE
      DO I = 1, NBIN
         IF( SDRCOUNT(I) .GT. 0 ) THEN
            TEMP = SDRVECTOR(I)/SDRCOUNT(I)
            WRITE(3,*) BINS(I), TEMP
         ELSE
            PRINT *, 'EMPTY BIN @ I = ', I
         ENDIF
      ENDDO
C
      END PROGRAM
c
c -----------------------SUBROUTINES--------------------------
c
c
c
c ------------------------------------------------------------
c
C NORMALIZE DATA
C
      SUBROUTINE NORMDATA(N,I,C_VEC)
C
      INTEGER N,I,II
      DOUBLE PRECISION T_MAX,T_MIN,C_VEC(N,*),DIFF
C
      T_MAX = C_VEC(N,I)
      T_MIN = C_VEC(1,I)
      DO II = 2,N
         IF( C_VEC(II,I) > T_MAX ) T_MAX = C_VEC(II,I)
         IF( C_VEC(II,I) < T_MIN ) T_MIN = C_VEC(II,I)            
      ENDDO
      DIFF = T_MAX - T_MIN
      DO II = 1, N
         C_VEC(II,I) = (C_VEC(II,I)-T_MIN)/DIFF
      ENDDO
C
      RETURN
C     
      END
C
c ------------------------------------------------------------
c
CCC CALCULATE GRADIANTS OF C AND STORE VALUES
C
      SUBROUTINE GET_GRADIENTS(N,I,TDX,C_VEC,GRADC_VEC)
C
      INTEGER N,I,II
      DOUBLE PRECISION GRADC_VEC(N-2,*),C_VEC(N,*),TDX
C
      DO II = 1, (N-2)
         GRADC_VEC(II,1) = C_VEC(II+1,I)
         GRADC_VEC(II,2) = ABS(C_VEC(II+2,I) - C_VEC(II,I))/TDX
      ENDDO
C   
      RETURN
C
      END
C
C -----------------------------------------------------------
C
CCC SORT THE GRADIENTS
C
      SUBROUTINE SORT(N,NBIN,BINSB,GRADC_VEC,SDRVECTOR,SDRCOUNT)
C
      INTEGER N,NC,II,III,SDRCOUNT(*)
      DOUBLE PRECISION GRADC_VEC(N-2,2),BINSB(*),
     1                 TEMP,SDRVECTOR(*)
      LOGICAL SORTED(N-2),VECTORENTRY(NBIN)
C SET LOGICAL VALUES TO FALSE 
      DO II = 1, (N-2)
         SORTED(II) = .FALSE.
      ENDDO
      DO II = 1, NBIN 
         VECTORENTRY = .FALSE.
      ENDDO
c BEGIN SEARCH
      DO III = 1, NBIN!mean
         NC = 0! counter for number of pdfs put into specific bin
         TEMP = 0D0
         DO II = 1, N-2
            IF (SORTED(II) .EQ. .FALSE.) THEN
               IF(GRADC_VEC(II,1) .GE. BINSB(III)   .AND.
     1            GRADC_VEC(II,1) .LT. BINSB(III+1) ) THEN
                  TEMP=TEMP+GRADC_VEC(II,2)
                  NC = NC + 1!add one SDR to counter
                  SORTED(II) = .TRUE.
                  VECTORENTRY(III) = .TRUE.
               ENDIF
            ENDIF
         ENDDO
C AVERAGE CONTRIBUTIONS AND STORE IN SDRVECTOR
         IF( NC .GT. 0) THEN
            SDRVECTOR(III) = SDRVECTOR(III)+TEMP/DFLOAT(NC)
         ENDIF           
      ENDDO
C USE LOGICAL TO COUNT UP CONTRIBUTIONS
      DO III = 1, NBIN
         IF( VECTORENTRY(III) .EQ. .TRUE. ) THEN
            SDRCOUNT(III) = SDRCOUNT(III) + 1
         ENDIF
      ENDDO
C
      RETURN
C
      END
         

















      
