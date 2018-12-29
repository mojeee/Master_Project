c %%%%%%%%%%%%%%%%%%%%%%
c
c BUILD PDFS WITH LEM SOLUTION
c %%%%%%%%%%%%%%%%%%%%%%
c
c bmp10.f
c
c Modified bmp9.f such that it works on clusters.
c New parameter is NTHREADS, where it must match number of
c cores designated @ bqsubmit.dat for clusters.
c
c
c bmp9.f
c
c Uses a parameter "SPACE" to space out the pdf construction.
c This is needed for large data files, or else not enough
c memory. Essentially bmp8 or older uses a "SPACE" of 1 by 
c default. Now it's adjustable.
c
c bmp8.f 
c
c changed normalized data from O2 to temperature
c
c
c bmp7.f
c
c parallel code, designed for 8 threads
c if segmentation fault happens, increase stack size using
c the following command in terminal:
c ulimit -s 2048000
c
c
c bmp6.f
c
c slight modifications for LEM
c -normalizes each realization individually before proceeding
c
c
c
c
c bmp5.f
c
c Makes sure all PDFs do not have gaps of 0 due
c to insufficient points in the laminar flame.
c Interpolates empty bins within the distribution
c linearly.
c
c bmp4.f:
c
c almost the same as bmp3 (below)
c changes:
c reverting back to 51 bins first bin @ c = 0, last bin @ c = 1
c   interior bins centered at 0.02, 0.04, etc, 49 of them
c   total 51 bins
c progress variable is now O2 instead of CO2
c   normalization:
c            ( O_ini - O_@x ) / ( O_ini - O_final )
c
c
c bmp3.f:
c
c PDF bins are constant at dc = 0.02 width, 
c centered @ c = 0.01 + dc for all 50 bins: 0.01, 0.03, 0.05, ... 0.99
C
c PDF means are tabulated from 0.01 to 0.99 using bin widths of 0.02
c PDF variances (normalized using segragation) are tabulated from
c    0.01 to 0.99 using bin widths of 0.02
c end result is a 50 x 50 table of means and variances
c
C LINES WITH **** SHOULD BE LOOKED AT / USED WITH CAUTION
C
c
c THEORY:
c
c first and last bins are special (done according to Jin Bei's 2008 PDF CSE/PDF paper)
c first bin includes: delta function at c = 0 and interior pdf for 0.01 < c < 0.02
c interior bin boundaries: 0.02 < c < 0.04, 0.04 < c < 0.06, ... 0.96 < c < 0.98
c last bin includes: 0.98 < 0.99 and delta function at c = 1
c
c IMPLEMENTATION:
c
c 1) reads simulations from interpedc.csv (CO2 data) and x2.csv (x coordinate data)
c 2) normalizes interpedc and calcs gradients of c
c       - store in matrix: CGC
c 3) build all cases (I-IV) for each realization individually
c 4) loop through the realizations one at a time to build/average 
c    PDFs from each case
c
c
c
c %%%%%%%%%%%%%%%%%%%%%%%

      PROGRAM BUILDPDF
      IMPLICIT NONE

      INTEGER N, NBIN, DELTA_CUTOFF, ROW,
     1        XMAX, XMIN, X1, X2,
     2        I, II, III, J, JJ, NC, BIG, K, KK,
     3        NSIM,PDFCOUNT(50,50), S1, SE, TID,
     4        OMP_GET_THREAD_NUM,NTHREADS,
     5        SPACE
      PARAMETER(N = 2122, NBIN = 51, NSIM = 3500, 
     1    DELTA_CUTOFF = 100, BIG = 400000, SPACE = 2,
     2    NTHREADS = 24)!messing with BIG values
      DOUBLE PRECISION UL,T_MIN,T_MAX, 
     1                 DC,CMIN,CMAX,C_VEC(N,NSIM),
     2                 DX,TDX,X_VEC(N),
     3                 CGC(N-2,2),BINS(NBIN+1),SOL(NBIN,2),
     4                 MB(50),VB(50),MBB(51),VBB(51),
     5                 TEMP,TEMP_VEC(N),
     6                 A,B,C,I0,I1,I2,
     7                 PDFTABLE(50,50,NBIN),
     8                 ALLCASE(NTHREADS,BIG,NBIN+2),
     9                 SEG
      PARAMETER(UL = 1D6)
      LOGICAL SORTED(BIG),TABLEENTRY(50,50)
CCC SET BINS VECTOR AND SOL MATRIX
      DC = 1.D0/( DFLOAT(NBIN) )
      CMIN = DC/2D0
      CMAX = 1-DC/2D0
      BINS(1) = DC/2d0!border of bins
c USING JIN'S DEFINITION, THERE WILL ONLY BE 49 interior bins.
c This means we'll only need 50 interior borders.
c left border of bin one @ DC/2 ****
c there is method to this madness
c other interior bins have borders @ 0.03, 0.05, +DC ...
      DO J = 2, NBIN-1
         BINS(J) = BINS(1) + DC*(J-1)
      ENDDO
c last bin has right border @ 1-DC/2 ****
c      BINS(NBIN) = 1D0-DC/2d0
      SOL(1,1) = 0D0! FOR CONVENIENCE ****
      DO J = 2, (NBIN-1)
         SOL(J,1) = (J-1)*DC
      ENDDO
      SOL(NBIN,1) = 1D0! AGAIN, FOR CONVENIENCE ****
CCC OPEN SOME FILES 
      OPEN(1, file="x.csv" , status="unknown")!distance data 
      OPEN(2, file="interpedc.csv", status="unknown")!flamelet data
      OPEN(3, FILE="PDF.csv", status="unknown")
      OPEN(8, FILE="sum.csv", status="unknown")
      OPEN(9, FILE="test.csv", status="unknown")
CCC READ LEM DATA INTO ARRAYS
      DO J = 1,N
         READ(1,*) X_VEC(J)
      ENDDO
      DO I = 1, NSIM
         DO J = 1,N
            READ(2,*) C_VEC(J,I)
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
CCC INITIALIZE PDFCOUNT MATRIX AS ZEROS
ccc (counts number of PDFs inclded for each mean and variance bin)
      DO I = 1, 50
         DO J = 1, 50
            PDFCOUNT(I,J) = 0
         ENDDO
      ENDDO
C   SET MB, MBB, BV, VBB (mean bin, mean bin border, var bin, etc)
      DO I = 1, 51!BORDERS
         MBB(I) = (I-1)*0.02d0
         VBB(I) = (I-1)*0.02d0
      ENDDO
c purposely make the last variance bin border larger so
c possible ALL solutions will be averaged (in case some
c pdf variances exceed a normalized value of 1)
      VBB(51) = 1.2D0! ****
      DO I = 1, 50!BINS
         MB(I) = 0.01D0 + (I-1)*0.02D0
         VB(I) = 0.01D0 + (I-1)*0.02D0
      ENDDO
C SET THE DEFAULT PDFTABLE VALUES TO 0.D0
      DO I = 1,50
         DO J = 1,50
            DO K = 1,NBIN
               PDFTABLE(I,J,K) = 0.D0
            ENDDO
         ENDDO
      ENDDO
C
C
CCC LOOP THROUGH ALL NSIMS WITH THE SAME PDF BUILDING PROCESS
C
      CALL OMP_SET_NUM_THREADS(NTHREADS)
!$OMP PARALLEL FIRSTPRIVATE(T_MAX,T_MIN,I,II,III,J,JJ,KK,
!$OMP& TEMP_VEC,TEMP,TID,XMIN,XMAX,X1,X2,CGC,ROW,SOL,NC,SEG,A,B,C,I0,
!$OMP& TABLEENTRY,SORTED)
!$OMP DO SCHEDULE(DYNAMIC,1)
      DO JJ = 1, NSIM
         TID = OMP_GET_THREAD_NUM()
         DO K = 0, (NTHREADS-1)
            IF ( TID == K ) THEN
               CALL PDF(N,JJ,XMIN,XMAX,ROW,NBIN,BIG,SPACE,NTHREADS,
     1              PDFCOUNT,TID,MB,MBB,VBB,BINS,C_VEC,ALLCASE,
     2              PDFTABLE,CMIN,CMAX,UL,TDX,DC,CGC,SOL)
            ENDIF
         ENDDO
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
      PRINT *, 'FINISHED BUILDING PDFS, NOW AVERAGING SOLUTIONS'
C   WRITE SOME HEADINGS TO OUTPUT FILE
      WRITE(3,200) ' Dimensions'
 200  FORMAT('#',A)
      WRITE(3,*) '51 PDF bins ', '50 mean bins ', '50 variance bins'

      WRITE(3,200) ' ctilde'
      WRITE(3,300) (MB(I),I=1,50)
 300  FORMAT(6(f5.2))

      WRITE(3,200) ' cvar'
      WRITE(3,300) (VB(I),I=1,50)
CCC AVERAGE PDF ACROSS REALIZATIONS
      NC = 0!count blank bins
      DO I = 1, 50
         DO J = 1, 50  
            TEMP = 0d0
            WRITE(3,400) MB(I), VB(J)
 400        FORMAT('# mean =', f5.2, ' variance =', f5.2)
            DO K = 1, NBIN
               TEMP_VEC(K) = 0d0
            ENDDO 
            IF( PDFCOUNT(I,J) .GT. 0 ) THEN
               DO K = 1, NBIN
                  TEMP_VEC(K) = PDFTABLE(I,J,K)/PDFCOUNT(I,J)
                  TEMP = TEMP + TEMP_VEC(K)
               ENDDO 
               WRITE(8,*) I, J, TEMP
            ELSE
               WRITE(8,*) I, J, 0D0
               NC = NC + 1
CCC WRITE TO FILE
            ENDIF
            WRITE(3,500) (TEMP_VEC(K),K=1,NBIN)
 500        FORMAT(6(f20.16))
            WRITE(9,600) (TEMP_VEC(K),K=1,NBIN)
 600        FORMAT(51(f20.16))
CCC CHECK FOR WEIRD EMPTY BINS IN THE MIDDLE OF THE PDF
ccc CHECK FOR FIRST BIN WITH DATA
            S1 = 2
            DO WHILE ( SOL(S1,2) == 0D0 )
               S1 = S1 + 1
            ENDDO
CCC CHECK FOR LAST BIN WITH DATA
            SE = NBIN-1
            DO WHILE ( SOL(SE,2) == 0D0 )
               SE = SE - 1
            ENDDO
CCC CHECK ALL BINS IN THE MIDDLE TO MAKE SURE THERE AIN'T EMPTIES
            DO K = S1,SE
               IF( SOL(K,2) .EQ. 0D0 ) THEN
                  WRITE(8,*) 'ERROR IN PDF, EMPTY BIN'
               ENDIF
            ENDDO
            
         ENDDO
      ENDDO
      PRINT *, 'Number of empty bins: ', NC
CCC DONE!
      END PROGRAM
c
c
c --------------------SUBROUTINES--------------------------
c
c
c ---------------------------------------------------------
c
C MAKES SURE SOL (PDF) DOESNT HAVE GAPS WITHIN BINS
C CALCULATES I0
C
C RETURNS UPDATED SOL AND IO VALUES
C
      SUBROUTINE CALC_I0(I0,NBIN,DC,SOL)
c
      INTEGER NBIN, I, S1, S2, SE
      DOUBLE PRECISION I0, DC, SOL(NBIN,2)
C 
ccc CHECK FOR FIRST BIN WITH DATA
      S1 = 2
      DO WHILE ( SOL(S1,2) == 0D0 )
         S1 = S1 + 1
      ENDDO
CCC CHECK FOR LAST BIN WITH DATA
      SE = NBIN-1
      DO WHILE ( SOL(SE,2) == 0D0 )
         SE = SE - 1
      ENDDO
CCC CHECK IN BETWEEN THE LIMITS TO MAKE SURE NOTHING IS 0D0
CCC IF ANY BIN IS EMPTY, INTERPOLATE WITH NEIGHBOURS!!
      DO I = (S1+1),(SE-1)
         IF( SOL(I,2) == 0D0 ) THEN
            S2 = I+1
            DO WHILE( SOL(S2,2) .EQ. 0D0 ) 
               S2 = S2 + 1
            ENDDO
            SOL(I,2) = SOL(I-1,2)+(SOL(S2,2) - SOL(I-1,2))/(S2-I+1)
         ENDIF
      ENDDO
CCC CALCULATE I0 FOR NORMALIZATION!
      I0 = 0d0
      DO I = 2,(NBIN-1)
         I0 = I0 + SOL(i,2)
      ENDDO
      I0 = I0*DC   
C
      RETURN
C
      END         
c
c ------------------------------------------------------------
c
C CALCULATES ALL INSTANCES OF CASE ONE PDFS - NO DELTAS
C STORES THEM IN ALLCASE MATRIX
C
      SUBROUTINE CASE1(TID,BIG,N,XMIN,XMAX,NBIN,SPACE,NTHREADS,
     &           BINS,DC,CGC,SOL,ALLCASE,ROW)
C
      INTEGER X1,X2,XMIN,XMAX,I,J,NC,KK,ROW,N,BIG,NBIN,TID,
     1        SPACE,NTHREADS
      DOUBLE PRECISION TEMP,SOL(NBIN,2),CGC(N-2,2),DC,I0,B,
     1                 ALLCASE(NTHREADS,BIG,NBIN+2),BINS(*)
c
      ROW = 1
C      DO X2 = (XMIN+3), (XMAX-1)
      X2 = (XMIN+3)
      DO WHILE (X2 .LE. (XMAX-1))
         X1 = (XMIN+1)
C         DO X1 = (XMIN+1), (X2-1)
         DO WHILE(X1 .LE. (X2-1))
            IF( X2 - X1 .GE. 10 ) THEN
               DO I = 1, NBIN-2
                  SOL(I+1,2) = 0D0
                  TEMP = 0D0   
                  NC = 0
                  DO J = X1, X2
                     IF( CGC(J,1) .GE. BINS(I) .AND.
     1                   CGC(J,1) .LT. BINS(I+1) ) THEN
                        TEMP = TEMP + CGC(J,2)
                        NC = NC + 1
                     ENDIF   
                  ENDDO
                  IF (NC .GT. 0) THEN
                     SOL(I+1,2) = TEMP/NC
                  ENDIF
               ENDDO
CCC CALCULATE I0
               CALL CALC_I0(I0,NBIN,DC,SOL)
               IF( I0 .gt. 0d0 ) then
CCC CALCULATE B
                  B = 1D0/I0;
CCC CALCULATE MEAN
                  ALLCASE((TID+1),ROW,NBIN+1) = 0D0
                  DO I = 2,(NBIN-1)
                     ALLCASE((TID+1),ROW,NBIN+1) = 
     1                      ALLCASE((TID+1),ROW,NBIN+1) 
     2                      + B*SOL(I,1)*SOL(I,2)
                  ENDDO
                  ALLCASE((TID+1),ROW,NBIN+1) = 
     1               ALLCASE((TID+1),ROW,NBIN+1)*DC
CCC CALCULATE VARIANCE
                  ALLCASE((TID+1),ROW,NBIN+2) = 0D0
                  DO I = 2,(NBIN-1)        
                     ALLCASE((TID+1),ROW,NBIN+2) = 
     1                  ALLCASE((TID+1),ROW,NBIN+2)
     2                  + B*SOL(I,2)*( SOL(I,1) 
     3                  - ALLCASE((TID+1),ROW,NBIN+1) )**2
                  ENDDO
                  ALLCASE((TID+1),ROW,NBIN+2) = 
     1               ALLCASE((TID+1),ROW,NBIN+2)*DC
CCC STORE PDF IN CASE VECTOR
                  ALLCASE((TID+1),ROW,1) = 0d0
                  DO KK = 2, (NBIN-1) 
                     ALLCASE((TID+1),ROW,KK) = B*SOL(KK,2)*DC
                  ENDDO
                  ALLCASE((TID+1),ROW,NBIN) = 0d0
CCC WRITE DATA TO FILE
C                  WRITE(3,*) ALLCASE((TID+1),ROW,NBIN+1),ALLCASE((TID+1),ROW,NBIN+2)
                  ROW = ROW + 1
               ENDIF
            ENDIF
            X1 = X1 + SPACE
         ENDDO
         X2 = X2 + SPACE
      ENDDO
c
      RETURN
C
      END
c
c ------------------------------------------------------------
c
C CALCULATES ALL INSTANCES OF CASE TWO PDFS - LEFT DELTA
C STORES THEM IN ALLCASE MATRIX
C
      SUBROUTINE CASE2(TID,BIG,N,XMIN,XMAX,NBIN,SPACE,NTHREADS,
     &           BINS,DC,CGC,SOL,ALLCASE,ROW)
C
      INTEGER X1,X2,XMIN,XMAX,I,J,NC,KK,ROW,N,BIG,NBIN,TID,
     1        SPACE,NTHREADS
      DOUBLE PRECISION TEMP,SOL(NBIN,2),CGC(N-2,2),DC,I0,A,B,
     1                 ALLCASE(NTHREADS,BIG,NBIN+2),BINS(*)
C
      X1 = XMIN
      X2 = (XMIN+3)
C      DO x2 = (xmin+3),(xmax-1)
      DO WHILE( X2 .LE. (XMAX-1))
c sort gc into bins and average
         DO I = 1, NBIN-2
            SOL(I+1,2) = 0D0
            TEMP = 0D0   
            NC = 0
            DO J = X1, X2
               IF( CGC(J,1) .GE. BINS(I) .AND.
     1             CGC(J,1) .LT. BINS(I+1) ) THEN
                   TEMP = TEMP + CGC(J,2)
                   NC = NC + 1
               ENDIF   
            ENDDO
            IF (NC .GT. 0) THEN
               SOL(I+1,2) = TEMP/NC
            ENDIF
         ENDDO
c calculate I's:
         CALL CALC_I0(I0,NBIN,DC,SOL)
c calculate A and B
         DO ii = 0,50
            A = 0.02d0*ii
            B = (1-A)/I0
CCC CALCULATE MEAN
            ALLCASE((TID+1),ROW,NBIN+1) = 0D0
            DO I = 2,(NBIN-1)
               ALLCASE((TID+1),ROW,NBIN+1) = 
     1            ALLCASE((TID+1),ROW,NBIN+1) 
     2            + B*SOL(I,1)*SOL(I,2)
            ENDDO
            ALLCASE((TID+1),ROW,NBIN+1) = 
     1         ALLCASE((TID+1),ROW,NBIN+1)*DC
CCC CALCULATE VARIANCE
            ALLCASE((TID+1),ROW,NBIN+2) = 0D0
            DO I = 2,(NBIN-1)        
               ALLCASE((TID+1),ROW,NBIN+2) = 
     1            ALLCASE((TID+1),ROW,NBIN+2)
     2          + B*SOL(I,2)*( SOL(I,1) 
     3          - ALLCASE((TID+1),ROW,NBIN+1) )**2
            ENDDO
            ALLCASE((TID+1),ROW,NBIN+2) = 
     1         ALLCASE((TID+1),ROW,NBIN+2)*DC
     2         + A*( ALLCASE((TID+1),row,NBIN+1) )**2!add delta function @ c = 0
CCC STORE PDF IN CASE VECTOR
            ALLCASE((TID+1),ROW,1) = A
!FIRST BIN IS COMBINATION OF DELTA FUNCTION AND INTERIOR PDF
            DO KK = 2, (NBIN-1) 
               ALLCASE((TID+1),ROW,KK) = B*SOL(KK,2)*DC
            ENDDO
            ALLCASE((TID+1),ROW,NBIN) = 0d0
CCC WRITE DATA TO FILE
c            WRITE(4,*) ALLCASE((TID+1),ROW,NBIN+1), ALLCASE((TID+1),ROW,NBIN+2)
            ROW = ROW + 1
         ENDDO
         X2 = X2 + SPACE
      ENDDO
C      
      RETURN
C
      END
c
c ------------------------------------------------------------
c
C CALCULATES ALL INSTANCES OF CASE THREE PDFS - RIGHT DELTA
C STORES THEM IN ALLCASE MATRIX
C
      SUBROUTINE CASE3(TID,BIG,N,XMIN,XMAX,NBIN,SPACE,NTHREADS,
     &           BINS,DC,CGC,SOL,ALLCASE,ROW)
C
      INTEGER X1,X2,XMIN,XMAX,I,J,NC,KK,ROW,N,BIG,NBIN,TID,
     1        SPACE,NTHREADS
      DOUBLE PRECISION TEMP,SOL(NBIN,2),CGC(N-2,2),DC,I0,B,C,
     1                 ALLCASE(NTHREADS,BIG,NBIN+2),BINS(*)
C
      X2 = XMAX
      X1 = (XMIN+1)
C      DO x1 = (xmin+1),(xmax-1)
      DO WHILE( X1 .LE. (XMAX-1))
C Sort gc into bins and average
         DO I = 1, NBIN-2
            SOL(I+1,2) = 0D0
            TEMP = 0D0   
            NC = 0
            DO J = X1, X2
               IF( CGC(J,1) .GE. BINS(I) .AND.
     1             CGC(J,1) .LT. BINS(I+1) ) THEN
                   TEMP = TEMP + CGC(J,2)
                   NC = NC + 1
               ENDIF   
            ENDDO
            IF (NC .GT. 0) THEN
               SOL(I+1,2) = TEMP/NC
            ENDIF
         ENDDO
c calculate I's:
         CALL CALC_I0(I0,NBIN,DC,SOL)
         IF( I0 .gt. 0d0 ) then
c calculate B and C
            DO II = 0,50
               C = 0.02*II
               B = (1d0-C)/I0
CCC CALCULATE MEAN
               ALLCASE((TID+1),ROW,NBIN+1) = 0d0
               DO I = 2,(NBIN-1)
                  ALLCASE((TID+1),ROW,NBIN+1) = 
     1               ALLCASE((TID+1),ROW,NBIN+1) 
     2               + B*SOL(I,1)*SOL(I,2)
               ENDDO
               ALLCASE((TID+1),ROW,NBIN+1) = 
     1            ALLCASE((TID+1),ROW,NBIN+1)*DC + C
CCC CALCULATE VARIANCE
               ALLCASE((TID+1),ROW,NBIN+2) = 0D0
               DO I = 2,(NBIN-1)        
                  ALLCASE((TID+1),ROW,NBIN+2) = 
     1               ALLCASE((TID+1),ROW,NBIN+2)
     2             + B*SOL(I,2)*( SOL(I,1) 
     3             - ALLCASE((TID+1),ROW,NBIN+1) )**2
               ENDDO
               ALLCASE((TID+1),ROW,NBIN+2) = 
     1            ALLCASE((TID+1),ROW,NBIN+2)*DC
     2          + C*( 1d0 - ALLCASE((TID+1),row,NBIN+1) )**2
                 !add delta function @ c = 0
CCC STORE PDF IN CASE VECTOR
               ALLCASE((TID+1),ROW,1) = 0d0
               DO KK = 2, (NBIN-1) 
                  ALLCASE((TID+1),ROW,KK) = B*SOL(KK,2)*DC
               ENDDO
               ALLCASE((TID+1),ROW,NBIN) = C
CCC WRITE DATA TO FILE
c               WRITE(5,*) ALLCASE((TID+1),ROW,NBIN+1), ALLCASE((TID+1),ROW,NBIN+2)
               ROW = ROW + 1
            ENDDO
         ENDIF
         X1 = X1 + SPACE
      ENDDO
C
      RETURN
C
      END
c
c ------------------------------------------------------------
c
C CALCULATES ALL INSTANCES OF CASE FOUR PDFS - DOUBLE DELTAS
C STORES THEM IN ALLCASE MATRIX
C
      SUBROUTINE CASE4(TID,BIG,N,XMIN,XMAX,NBIN,SPACE,NTHREADS,
     &           BINS,DC,CGC,SOL,ALLCASE,ROW)
C
      INTEGER X1,X2,XMIN,XMAX,I,J,NC,KK,ROW,N,BIG,NBIN,TID,
     1        SPACE,NTHREADS
      DOUBLE PRECISION TEMP,SOL(NBIN,2),CGC(N-2,2),DC,I0,A,B,C,
     1                 ALLCASE(NTHREADS,BIG,NBIN+2),BINS(*)
C 
      X1 = XMIN
      X2 = XMAX
c sort gc into bins and average
      DO I = 1, NBIN-2
c         TID = OMP_GET_THREAD_NUM()
c         WRITE(99,*) TID, ROW
         SOL(I+1,2) = 0D0
         TEMP = 0D0   
         NC = 0
         DO J = X1, X2
            IF( CGC(J,1) .GE. BINS(I) .AND.
     1          CGC(J,1) .LT. BINS(I+1) ) THEN
                TEMP = TEMP + CGC(J,2)
                NC = NC + 1
            ENDIF   
         ENDDO
         IF (NC .GT. 0) THEN
            SOL(I+1,2) = TEMP/NC
         ENDIF
      ENDDO
c calculate I's:
      CALL CALC_I0(I0,NBIN,DC,SOL)
c calculate B and C
      DO II = 0,50
         DO III = 0,50
            A = 0.02d0*DFLOAT(II)
            C = 0.02d0*DFLOAT(III)
            IF( (A + C) .LE. 1D0 ) THEN
               B = (1d0-A-C)/I0
CCC CALCULATE MEAN
               ALLCASE((TID+1),ROW,NBIN+1) = 0d0
               DO I = 2,(NBIN-1)
                  ALLCASE((TID+1),ROW,NBIN+1) = 
     1               ALLCASE((TID+1),ROW,NBIN+1) 
     2               + B*SOL(I,1)*SOL(I,2)
               ENDDO
               ALLCASE((TID+1),ROW,NBIN+1) = 
     1            ALLCASE((TID+1),ROW,NBIN+1)*DC + C
CCC CALCULATE VARIANCE
               ALLCASE((TID+1),ROW,NBIN+2) = 0D0
               DO I = 2,(NBIN-1)        
                  ALLCASE((TID+1),ROW,NBIN+2) = 
     1               ALLCASE((TID+1),ROW,NBIN+2)
     2             + B*SOL(I,2)*( SOL(I,1) 
     3             - ALLCASE((TID+1),ROW,NBIN+1) )**2
               ENDDO
               ALLCASE((TID+1),ROW,NBIN+2) = 
     1            ALLCASE((TID+1),ROW,NBIN+2)*DC
     2          + A*( ALLCASE((TID+1),row,NBIN+1) )**2!add delta function @ c = 0
     3          + C*( 1d0 - ALLCASE((TID+1),row,NBIN+1) )**2!add delta function @ c = 1
CCC STORE PDF IN CASE VECTOR
!FIRST BIN IS COMBINATION OF DELTA FUNCTION AND INTERIOR PDF
               ALLCASE((TID+1),ROW,1) = A
               DO KK = 2, (NBIN-1) 
                  ALLCASE((TID+1),ROW,KK) = abs(B)*SOL(KK,2)*DC
               ENDDO
               ALLCASE((TID+1),ROW,NBIN) = C
CCC WRITE DATA TO FILE
c               WRITE(6,*) ALLCASE((TID+1),ROW,NBIN+1), ALLCASE((TID+1),ROW,NBIN+2)
               ROW = ROW + 1
            ENDIF
         ENDDO
      ENDDO
C
      RETURN
C 
      END
c
c ------------------------------------------------------------
c
C NORMALIZE DATA
C
      SUBROUTINE NORMDATA(N,JJ,C_VEC)
C
      INTEGER N,J,JJ
      DOUBLE PRECISION T_MAX,T_MIN,C_VEC(N,*),DIFF
C
      T_MAX = C_VEC(N,JJ)
      T_MIN = C_VEC(1,JJ)
      DO J = 2,N
         IF( C_VEC(J,JJ) > T_MAX ) T_MAX = C_VEC(J,JJ)
         IF( C_VEC(J,JJ) < T_MIN ) T_MIN = C_VEC(J,JJ)            
      ENDDO
      DIFF = T_MAX - T_MIN
      DO J = 1,N
         C_VEC(J,JJ) = (C_VEC(J,JJ)-T_MIN)/DIFF
      ENDDO
C
      RETURN
C     
      END
c
c ------------------------------------------------------------
c
CCC FIND WHERE IN X DOES C MATCH WITH CMIN AND CMAX
C
      SUBROUTINE GET_XMINMAX(N,JJ,CMIN,CMAX,C_VEC,XMIN,XMAX)
C
      INTEGER J,JJ,N,XMIN,XMAX
      DOUBLE PRECISION TEMP,TEMP_VEC(N),C_VEC(N,*),CMIN,CMAX
C
      DO J = 1, N
         TEMP_VEC(J) = ABS( C_VEC(J,JJ) - CMIN )
      ENDDO
      TEMP = TEMP_VEC(1)
      DO J = 2, N
         IF( TEMP_VEC(J) < TEMP) THEN
            TEMP = TEMP_VEC(J)
            XMIN = J
         ENDIF
      ENDDO
      DO J = 1, N
         TEMP_VEC(J) = ABS( C_VEC(J,JJ) - CMAX )
      ENDDO
      TEMP = TEMP_VEC(1)
      DO J = 2, N
         IF( TEMP_VEC(J) < TEMP) THEN
            TEMP = TEMP_VEC(J)
            XMAX = J
         ENDIF
      ENDDO
C
      RETURN
C
      END
c
c ------------------------------------------------------------
c
CCC CALCULATE GRADIANTS OF C AND STORE VALUES
C
      SUBROUTINE GET_GRADIENTS(N,JJ,UL,TDX,C_VEC,CGC)
C
      INTEGER N,NM2,X1,X2,J,JJ
      DOUBLE PRECISION CGC(N-2,2),C_VEC(N,*),TDX,UL
C
      NM2 = N-2
      X1 = 1
      X2 = NM2
      DO J = 1, NM2
         CGC(J,1) = C_VEC(J+1,JJ)
         CGC(J,2) = TDX/ABS(C_VEC(J+2,JJ) - C_VEC(J,JJ))
         IF( CGC(J,2) > UL ) CGC(J,2) = UL 
      ENDDO
C   
      RETURN
C
      END
c
c ------------------------------------------------------------
c
C SORT AND AVERAGE PDFS BY MEAN AND VARIANCE
C 50 BINS IN MEAN, 50 BINS IN VARIANCE
C 2500 BINS TOTAL
C CHECK ALL PDFS WITH EACH BIN, THEN AVG IF NECESSARY
C
      SUBROUTINE SORT_PDFS(TID,BIG,NBIN,ROW,MB,MBB,
     &           VBB,NTHREADS,ALLCASE,PDFTABLE,PDFCOUNT)
C
      INTEGER I,J,K,KK,ROW,NC,NBIN,BIG,PDFCOUNT(50,50),TID,NTHREADS
      DOUBLE PRECISION SEG,MB(*),TEMP_VEC(NBIN),
     1                 ALLCASE(NTHREADS,BIG,NBIN+2),
     2                 MBB(*),VBB(*),PDFTABLE(50,50,NBIN)
      LOGICAL TABLEENTRY(50,50),SORTED(BIG)
C SET logical table to false for all bins BEFORE SORTING
      DO I = 1, 50
         DO J = 1, 50
            TABLEENTRY(I,J) = .FALSE.
         ENDDO
      ENDDO
C SET SORTED TO FALSE BEFORE
      DO I = 1, ROW
         SORTED(I) = .FALSE.
      ENDDO
c begin search
      DO I = 1, 50!mean
         SEG = MB(I)*(1D0 - MB(I))
         DO J = 1, 50!variance
            NC = 0! counter for number of pdfs put into specific bin
            DO KK = 1, NBIN
               TEMP_VEC(KK) = 0d0!temperary vector
            ENDDO
            DO K = 1, ROW
               IF (SORTED(K) .EQ. .FALSE.) THEN
                  IF(ALLCASE((TID+1),K,NBIN+1) .ge. MBB(I)   .AND.
     1               ALLCASE((TID+1),K,NBIN+1) .lt. MBB(I+1) .AND.
     2               ALLCASE((TID+1),K,NBIN+2) .ge. VBB(J)*SEG   .AND.
     3               ALLCASE((TID+1),K,NBIN+2) .lt. VBB(J+1)*SEG) THEN
                     DO II = 1, NBIN
                        TEMP_VEC(II)=TEMP_VEC(II)+ALLCASE((TID+1),K,II)
                     ENDDO
                     NC = NC + 1!add one pdf to counter
                     SORTED(K) = .TRUE.
                     TABLEENTRY(I,J) = .TRUE.
                  ENDIF
               ENDIF
            ENDDO
C   AVERAGE CONTRIBUTIONS AND STORE IN PDFTABLE MATRIX
            IF( NC .GT. 0) THEN
               DO II = 1, NBIN
                  PDFTABLE(I,J,II) = PDFTABLE(I,J,II) 
     1                               + TEMP_VEC(II)/DFLOAT(NC)
               ENDDO
            ENDIF
         ENDDO
      ENDDO
C use the logical check to count up how many realizations
c contributed to each mean/variance bin
      DO i = 1, 50
         DO j = 1, 50
            IF( tableentry(i,j) .eq. .true. ) then
               pdfcount(i,j) = pdfcount(i,j) + 1
            ENDIF
         ENDDO
      ENDDO
CCC RESET SORTED TO FALSE FOR NEXT REALIZATION
      DO I = 1, ROW
         SORTED(I) = .FALSE.
      ENDDO
C
      RETURN
C
      END
c
c ------------------------------------------------------------
c
C BUILD PDF
C
      SUBROUTINE PDF(N,JJ,XMIN,XMAX,ROW,NBIN,BIG,SPACE,NTHREADS,
     1  PDFCOUNT,TID,MB,MBB,VBB,BINS,C_VEC,ALLCASE,
     2  PDFTABLE,CMIN,CMAX,UL,TDX,DC,CGC,SOL)
C
      INTEGER N,JJ,XMIN,XMAX,ROW,NBIN,BIG,PDFCOUNT(50,50),TID,
     1        SPACE,NTHREADS
      DOUBLE PRECISION MB(*),MBB(*),VBB(*),BINS(*),C_VEC(N,*),
     1       ALLCASE(NTHREADS,BIG,NBIN+2),PDFTABLE(50,50,NBIN),
     2       CMIN,CMAX,UL,TDX,DC,CGC(N-2,2),SOL(NBIN,2)  
C
      PRINT *, TID, JJ
CCC NORMALIZE DATA
      CALL NORMDATA(N,JJ,C_VEC)
CCC FIND WHERE IN X DOES C MATCH WITH CMIN AND CMAX
      CALL GET_XMINMAX(N,JJ,CMIN,CMAX,C_VEC,XMIN,XMAX)
CCC CALCULATE GRADIANTS OF C AND STORE VALUES
      CALL GET_GRADIENTS(N,JJ,UL,TDX,C_VEC,CGC)
CCC BUILD PDFS OF VARIOUS MEANS AND VARIANCES ACCORDING TO 4 CASES
C   CASE 1: XMIN < X1 < X2 < XMAX
      CALL CASE1(TID,BIG,N,XMIN,XMAX,NBIN,SPACE,NTHREADS,
     1           BINS,DC,CGC,SOL,ALLCASE,ROW)
c   CASE 2: x1 < xmin < x2 < xmax
      CALL CASE2(TID,BIG,N,XMIN,XMAX,NBIN,SPACE,NTHREADS,
     1           BINS,DC,CGC,SOL,ALLCASE,ROW)
c   CASE 3: xminn < x1 < xmaxx < x2
      CALL CASE3(TID,BIG,N,XMIN,XMAX,NBIN,SPACE,NTHREADS,
     1           BINS,DC,CGC,SOL,ALLCASE,ROW)
c   CASE 4: x1 < xminn < xmaxx < x2
      CALL CASE4(TID,BIG,N,XMIN,XMAX,NBIN,SPACE,NTHREADS,
     1           BINS,DC,CGC,SOL,ALLCASE,ROW)
CCC SORT AND AVERAGE PDFS BY MEAN AND VARIANCE
      CALL SORT_PDFS(TID,BIG,NBIN,ROW,MB,MBB,VBB,NTHREADS,
     1              ALLCASE,PDFTABLE,PDFCOUNT)
C
      RETURN
C  
      END

