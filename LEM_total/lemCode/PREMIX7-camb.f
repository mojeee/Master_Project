ccc
c
c PREMIX 
c v.7 - added features to deal with more variable number of cells/cm
C       flame stabilizes at around 1/3 of the domain from inlet
c     - NSIM and number of timesteps per sim is now read from input file
c     - input.dat needs to be modded as well - changed name to input7.dat
c v.6 - added integral length to params array
c     - can do multiple mappings per timestep
c v.5 - same as premix4, except optimized for collecting data
c       instead of testing
c v.4 - uses premixed equations (By Kee from Sandia) entirely
c     - cleaned the rest of the code (primarily aesthetics)
c     - fixed a bunch of subroutines to make them faster than before
c
c v.3 - in cgs units 
c     - uses chemkin for thermal conductivity and 
c       species diffusivities
c
c v.2 - solutions are now combined into 
c       2D arrays instead of separate vectors
c             - makes it easier to increase number of species
c     - when number of species changes, must edit the following places:
c             - parameter line for chemkin input
c             - readinput subroutine - change total parameters to be read
c             - writeic/readic subroutines - change species matrices
c             - writefiles/writepv subroutines
c             - update_ities 3 
c             - wdot
c             - enthalpy
c             - advy - need to tweak fuel reaction procedure
c             - advt
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      PROGRAM PREMIX7
      implicit none
CCC CHEMKIN COMPONENT
      integer nspc, nspcm1,
     1 lrckwk, lickwk, lcckwk, lenimc, lenrmc
      parameter(nspc=10,nspcm1=nspc-1,
     1          lrckwk=227,lickwk=71,lcckwk=14,
     2          lenimc=42 ,lenrmc=2300)
      character*16 cckwk(lcckwk)
      character*10 filenm
      double precision rckwk(lrckwk), rmcwrk(lenrmc), XMWT(nspc)
      integer ickwk(lickwk), imcwrk(lenimc)
CCC PREMIX COMPONENTS
      integer nc,ncm1,ncp1,i,j,k,kk,NTS,NTS_COUNT,NFL,
     1 NSIM_COUNT, NSIM, NTSPSIM
      parameter(nc=2121,ncp1=nc+1,ncm1=nc-1) !nc = number of cells
      double precision 
     1 DX, DOM, TEND, DT, XMDT, WT, XMDOT, P, params(33), GFAC,
     2 X(ncp1),RHO(ncp1),U(ncp1),T(ncp1),YZ(nspc,ncp1),
     3 COND(ncp1),D(nspc,ncp1),xxwt
CCC LOGICAL FLAGS
      LOGICAL LRIC,LTM,LWI,LOWIC,LCF,LCT
ccc LEM component
      integer seed, NTS_PE, MTS
      double precision PDFA, PDFB
ccc INITIALIZE CHEMKIN AND OTHER PARAMETERS
c      OPEN(2, file="d.csv", status="unknown")!all data (for testing)
      OPEN(3, file="pv.csv",status="unknown")!progress variable data
c      OPEN(8, file='tester.csv', status="unknown")
      OPEN(9, file="z.dat",  status="unknown")
      OPEN(910,file='cklink',form='unformatted',status='unknown')
      OPEN(920,file='tplink',form='unformatted',status='unknown')
      CALL CKINIT(lickwk,lrckwk,lcckwk,910,911,ickwk,rckwk,cckwk)
      CALL MCINIT(920,921,lenimc,lenrmc,imcwrk,rmcwrk)
      CALL INIT_ALL(NSPC,NCP1,NCM1,NFL,PARAMS,GFAC,
     1        LRIC,LTM,LWI,LOWIC,LCF,LCT,
     2        XMWT,RHO,U,T,YZ,XMDOT,P,
     3        lickwk,lrckwk,lcckwk,ickwk,rckwk,cckwk,
     4        lenimc,lenrmc,imcwrk,rmcwrk,
     5        TEND,DT,NTS,XMDT,DOM,DX,X,
     6        NSIM,NTSPSIM)
ccc   INITIALIZE LEM PARAMETERS
      CALL INIT_LEM(NC,PARAMS,DOM,TEND,DT,U,
     1           NTS_COUNT,NTS_PE,PDFA,PDFB,MTS)
C
      CALL INIT_RANDOM_SEED()       
ccc write down conditions before first timestep
      IF( LWI ) CALL WRITEFILES(NSPC,NCP1,RHO,U,T,YZ)
      print *, 'And so we start...'
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
c
c                   BEGIN TIMESTEPPING              
c
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      NSIM_COUNT = 0
C IF DEAD, READ IC AND RUN AGAIN
 100  CALL READIC(NSPC,NCP1,DX,XMDOT,X,RHO,U,T,YZ)
      DO WHILE ( NSIM_COUNT .LT. NSIM )
c 
         DO i = 1, NTSPSIM
C PREMIX ADVANCEMENT IN TIME
            CALL PREMIXADV(i,nspc,nc,dx,dt,P,XMDOT,XMWT,GFAC,
     1           ICKWK, RCKWK, IMCWRK, RMCWRK, T, YZ, COND, D)
C TRIPLET MAP
            IF( LTM .AND. mod(NTS_COUNT,NTS_PE) .eq. 0  ) THEN
               CALL TM(NSPC,NC,NCP1,DX,PDFA,PDFB,T,YZ,MTS)
            ENDIF
C CHECK/ADJUST FUEL FLOW
            IF( LCF .AND. mod(i,100) .eq. 0 ) THEN
               CALL CFUEL(NSPC,NCP1,NFL,DX,XMDOT,RHO(1),U(1),YZ)
            ENDIF
C increase timestep counter by 1 (used by triplet map exclusively)
            NTS_COUNT = NTS_COUNT + 1
C FINISHED ONE REALIZATION
         ENDDO
C CHECK FOR ERROR IN SOLUTION VECTORS
         CALL CKMMWY(YZ(1,nc),ickwk,rckwk,XXWT)!weight in g/mole
         U(nc)   = XMDOT*(8.314462d7*T(nc))/(P*XXWT) 
         IF( .not. (U(nc) .le. 2.d4) .OR. (U(nc) .lt. 0.d0) ) THEN
            WRITE(9,*) 'DIED AT:', NSIM_COUNT
            PRINT *, 'DIED AT:', NSIM_COUNT
            GO TO 100
         ENDIF
C WRITE TO PV FILE
         CALL WRITEPV(NCP1,T)
C WRITE IC FILES IN CASE OF RESTART
         CALL GET_RHO_U(NC,NCP1,NSPC,XMDOT,T,YZ,P,ICKWK,RCKWK,RHO,U)
         CALL WRITEICFILES(NSPC,NCP1,X,RHO,U,T,YZ)
C CHECK FUEL CONSUMPTION AND SET INLET VELOCITY
C         IF( LCF ) THEN
C            CALL SET_VELOCITY(NSPC,NCP1,NFL,DX,XMDOT,GFAC,P,
C     1           ICKWK,RCKWK,RHO(1),U(1),XMWT,T,YZ)
C         ENDIF
C WRITE DATA TO D.CSV
c         IF( LWI ) THEN
c            CALL GET_RHO_U(NC,NCP1,NSPC,XMDOT,T,YZ,P,ICKWK,RCKWK,RHO,U)
c            CALL WRITEFILES(NSPC,NCP1,RHO,U,T,YZ)
c         ENDIF
         NSIM_COUNT = NSIM_COUNT + 1
         PRINT *, NSIM_COUNT!, 'SIMULATIONS COMPLETED'
C FINISHED SIMULATION
      ENDDO
      PRINT *, 'FINISHED ALL ALLOCATED SIMULATIONS'
C END PROGRAM
      END
C INCLUDE ADDITIONAL LIBRARIES/FILES
        include 'BTriplet.f'
        include 'htlib.f'
        include 'cklib.f'
        include 'tranlib.f'
        include 'math.f'
        include 'i1mach.f'
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C                    PREMIX SUBROUTINES                           
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C--------------------------------------------------------------------
C
      SUBROUTINE HT_MTRNPR (nspc, nc, P, T, YZ,
     1 ICKWRK, RCKWRK, IMCWRK, RMCWRK, COND, D)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C
      DIMENSION YAV(nspc),XAV(nspc),T(*),YZ(nspc,*),
     1          ICKWRK(*), RCKWRK(*), IMCWRK(*), RMCWRK(*),
     2          COND(*),D(nspc,*)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   nspc   - NUMBER OF CHEMICAL SPECIES.
C   nc     - NUMBER OF MESH POINTS.
C   p      - PRESSURE.
C              CGS UNITS - DYNES/CM**2
C   T      - TEMPERATURE
C              UNITS - K
C   YZ     - SPECIES MASS FRACTIONS
C              (NUMBER OF SPECIES, NC+1 CELLS)
C WORK AND SCRATCH SPACE
C   YAV    - ARRAY OF MASS FRACTIONS AT MESH MIDPOINTS.
C              DIMENSION Y(*) AT LEAST KK.
C   XAV    - ARRAY OF MOLE FRACTIONS AT MESH MIDPOINTS.
C              DIMENSION X(*) AT LEAST KK
C   ICKWRK - INTEGER CHEMKIN WORK SPACE.
C              DIMENSIONING - SEE CHEMKIN DOCUMENTATION.
C   RCKWRK  - FLOATING POINT CHEMKIN WORK SPACE.
C              DIMENSIONING - SEE CHEMKIN DOCUMENTATION.
C   IMCWRK - INTEGER TRANSPORT PROPERTY WORK SPACE.
C              DIMENSIONING - SEE TRANSPORT DOCUMENTATION.
C   RMCWRK  - FLOATING POINT TRANSPORT PROPERTY WORK SPACE.
C              DIMENSIONING - SEE TRANSPORT DOCUMENTATION.
C OUTPUT-
C   COND   - ARRAY OF CONDUCTIVITIES AT THE MESH MID-POINTS.
C   D      - MATRIX OF DIFFUSION COEFFICIENTS AT THE MESH MID-POINTS.
C              DIMENSION D(nspc,*) EXACTLY nspc FOR THE FIRST DIMENSION,
C              AND AT LEAST nc+1 FOR THE SECOND.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C-----------------MIXTURE-AVERAGED FORMULAS---------------------
C
         DO 400 j = 1, nc
            TAV = 0.5 * (T(j) + T(j+1))
Ccc
C              DIMENSIONAL TEMPERATURE AT THE MID POINTS
Ccc
            DO 300 k = 1, nspc
               YAV(k) = 0.5 * (YZ(k,j) + YZ(k,j+1))
300         CONTINUE
            CALL CKYTX (YAV, ICKWRK, RCKWRK, XAV)
            CALL MCADIF(p, TAV, XAV, RMCWRK, D(1,j) )!PREMIX DIFFU
Ccc
C              DETERMINE THE MIXTURE CONDUCTIVITY AT J
Ccc
            CALL MCACON( TAV, XAV, RMCWRK, COND(j) )
C
400      CONTINUE
C         
      RETURN
      END
C
C--------------------------------------------------------------------
C 
      SUBROUTINE HT_MDIFV(nspc,nc,dx,YZ,T,WT,p,D,ICKWRK,RCKWRK,YV)
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
c      COMMON /LOCS/ NT, NM, NYS, NY, NTR
C
      DIMENSION YZ(nspc,*),T(*),WT(nspc),
     1     YAV(nspc),D(nspc,*),YV(nspc,*),
     2     ICKWRK(*),RCKWRK(*),
     3     XMFP(nspc),XMF(nspc)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   nspc   - NUMBER OF CHEMICAL SPECIES.
C   nc     - NUMBER OF MESH POINTS.
C   p      - PRESSURE.
C              CGS UNITS - DYNES/CM**2
C   T      - TEMPERATURE
C              UNITS - K
C   YZ     - SPECIES MASS FRACTIONS
C              (NUMBER OF SPECIES, NC+1 CELLS)
C   WT     - THE ARRAY OF SPECIES MOLECULAR WEIGHTS.
C              CGS UNITS - GM/MOLE
C              DIMENSION WT(*) AT LEAST KK.
C WORK AND SCRATCH SPACE
C   YAV    - ARRAY OF MASS FRACTIONS AT MESH MIDPOINTS.
C              DIMENSION Y(*) AT LEAST nspc.
C   XMF    - ARRAY OF MOLE FRACTIONS AT MESH POINT j.
C              DIMENSION XMF(*) AT LEAST nspc.
C   XMFP   - ARRAY OF MOLE FRACTIONS AT MESH POINT j+1.
C              DIMENSION XMFP(*) AT LEAST nspc.
C   D      - THESE ARE COMPUTED EACH TIME THE FUNCTION IS CALLED,
C            OTHERWISE THE STORED VALUES ARE USED.
C              CGS UNITS - CM**2/SEC
C              DIMENSION D(nspc,*) EXACTLY nspc FOR THE FIRST DIMENSION
C              AND AT LEAST nc FOR THE SECOND.
C   ICKWRK - INTEGER CHEMKIN WORK SPACE.
C              DIMENSIONING - SEE CHEMKIN DOCUMENTATION.
C   RCKWRK  - FLOATING POINT CHEMKIN WORK SPACE.
C              DIMENSIONING - SEE CHEMKIN DOCUMENTATION.
C   IMCWRK - INTEGER TRANSPORT PROPERTY WORK SPACE.
C              DIMENSIONING - SEE TRANSPORT DOCUMENTATION.
C   RMCWRK  - FLOATING POINT TRANSPORT PROPERTY WORK SPACE.
C              DIMENSIONING - SEE TRANSPORT DOCUMENTATION.
c
C OUTPUT-
C   YV     - MATRIX OF MASS FRACTIONS TIMES DIFFUSION VELOCITIES AT THE
C            MESH MIDPOINTS.  YV(k,j) IS THE FLUX OF KTH SPECIES BETWEEN
C            j AND j+1.
C              DIMENSION YV(nspc,*) EXACTLY nspc FOR THE FIRST DIMENSION
C              AND AT LEAST nc FOR THE SECOND.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      CALL CKYTX (YZ(1,1), ICKWRK, RCKWRK, XMFP)
C
C           LOOP OVER ALL MESH POINTS, COMPUTING THE DIFFUSION
C         VELOCITY AT THE MID POINTS.  THE INDEXING IS SUCH THAT
C         YV(K,J) IS THE DIFFUSION VELOCITY OF THE KTH SPECIES
C         MIDWAY BETWEEN NODES J AND J+1.
C
      DO 1000 j = 1, nc
C
         TAV = 0.5 * (T(j) + T(j+1))
C
         DO 300 k = 1, nspc
            YAV(K) = 0.5 * (YZ(k,j) + YZ(k,j+1))
300      CONTINUE
C
         DO 400 k = 1, nspc
            XMF(k) = XMFP(k)
400      CONTINUE
         CALL CKMMWY (YAV, ICKWRK, RCKWRK, WTM)
         CALL CKRHOY (p, TAV, YAV, ICKWRK, RCKWRK, RHOAV)
         CALL CKYTX (yz(1,j+1), ICKWRK, RCKWRK, XMFP)
C
C           USE MIXTURE-AVERAGED FORM FOR FICKIAN DIFFUSION,
C           WHETHER WE ARE USING THE MULTICOMPONENT FORMALISM
C           OR MIXTURE-AVERAGED
C
         DO 500 k = 1, nspc
            YV(k,j) = - D(k,j)*(WT(k)/WTM)*(XMFP(k)-XMF(k))/dx
500      CONTINUE
C    
C        COMPUTE AND ADD THE CORRECTION VELOCITY
C
         SUM = 0.0
         DO 700 k = 1, nspc
            SUM = SUM + YV(k,j)
700      CONTINUE
C
         VC = - SUM
C
         DO 800 k = 1, nspc
            YV(k,j) = YV(k,j) + YAV(k)*VC
800      CONTINUE
C
1000  CONTINUE
C
      RETURN
      END
c
c ----------------------------------------------------------
c
      SUBROUTINE PREMIXADV(i,nspc,nc,dx,dt,p,XMDOT,XMWT,GFAC,
     1  ICKWRK, RCKWRK, IMCWRK, RMCWRK, T, YZ, COND, D)
C      
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C
      DIMENSION YAV(nspc),XAV(nspc),H(nspc),WDOT(nspc),CP(nspc),
     1          XMWT(nspc),YV(nspc,nc),F(3+nspc,nc),
     2          ICKWRK(*), RCKWRK(*), IMCWRK(*), RMCWRK(*),
     3          T(*),YZ(nspc,*),COND(*),D(nspc,*)
C
C EVALUATE AND STORE THE TRANSPORT COEFFICIENTS, OBTAIN CONDUCTIVITIES
      IF ( mod(i,20) .eq. 0 .OR. i .eq. 1 ) THEN
         CALL HT_MTRNPR (nspc, nc, P, T, YZ,
     1        ICKWRK, RCKWRK, IMCWRK, RMCWRK, COND, D)
      ENDIF
C
C EVALUATE AND STORE THE DIFFUSION VELOCITIES
C
      CALL HT_MDIFV(nspc,nc,dx,YZ,T,XMWT,P,D,ICKWRK,RCKWRK,YV)
C
ccc   FIRST CELL
      TAV = 0.5 * ( T(1) + T(2) )
      DO 10 k = 1, nspc
         YAV(k) = 0.5 * (YZ(k,1) + YZ(k,2))
10    CONTINUE
      CALL CKRHOY (P, TAV, YAV, ICKWRK, RCKWRK, RHOP)
Ccc
C                   INTERIOR MESH POINTS
Ccc
CCC   INTERIOR CELLS
      DO 1000 j = 2, nc
Ccc
        TAV = 0.5 * ( T(j) + T(j+1) )
        DO 100 k = 1, nspc
           YAV(k) = 0.5 * (YZ(k,j) + YZ(k,j+1))
100     CONTINUE
Ccc
        RHOM = RHOP
        CALL CKRHOY (P, TAV, YAV, ICKWRK, RCKWRK, RHOP)
Ccc
C             FORM THE CHEMICAL RATE TERMS
Ccc
        CALL JYCWYP (P, T(j), YZ(1,j), ICKWRK, RCKWRK, WDOT)
        DO 35 KK = 1, NSPC
           WDOT(KK) = WDOT(KK)*GFAC
   35   CONTINUE
        CALL CKHML (T(j), ICKWRK, RCKWRK, H)
        CALL CKCPBS (T(j), YZ(1,j), ICKWRK, RCKWRK, CPB)
        CALL CKCPMS (T(j), ICKWRK, RCKWRK, CP)
Ccc
C               SPECIES CONSERVATION EQUATION
Ccc
        SUMYK = 0.
        XMDXM = XMDOT / dx
        DO 200 k = 1, nspc
           SUMYK = SUMYK + YZ(j,k)
           F(3+k,j) = - dt*(
     1            XMDXM * ( YZ(k,j)-YZ(k,j-1) ) 
     2          + (RHOP*YV(k,j) - RHOM*YV(k,j-1))/dx 
     3          - WDOT(k)*XMWT(k)
     4                       )/((RHOP + RHOM)/2.0)
200     CONTINUE
Ccc
C               ENERGY EQUATION
Ccc
        SUMX = 0.0
        TDOT = 0.0
        DO 400 k = 1, nspc
           TDOT = TDOT + WDOT(k)*H(k)
           SUMX = SUMX + 0.25 * (RHOP*YV(k,j) + RHOM*YV(k,j-1)) *
     1                CP(k)*(T(j+1)-T(j-1))/dx
400     CONTINUE
Ccc
        F(3,j) = - dt*(
     1            XMDOT*( T(j)-T(j-1) )/dx 
     2        - ( COND(j  )*( T(j+1)-T(j) )/dx - 
     3            COND(j-1)*( T(j)-T(j-1) )/dx  )/(CPB*dx)
     4        + (SUMX + TDOT) /CPB
     5                    )/((RHOP + RHOM)/2.0)
c
1000  CONTINUE
ccc
C                UPDATE ARRAYS
CCC
      DO j = 2,nc
         T(j) = T(j) + F(3,j)
         DO k = 1,nspc
            YZ(k,j) = YZ(k,j) + F(3+k,j)
            IF(YZ(k,j) .lt. 0.d0 ) YZ(k,j) = 0.d0
         ENDDO
      ENDDO
CCC
C              RIGHT BOUNDARY - ZERO GRADIENT
CCC
      T(nc+1) = T(nc)
      DO kk = 1, nspc
         YZ(kk,nc+1) = YZ(kk,nc)
      enddo
C
      RETURN
C
      END
c
c -------------------------------------------------------
c
      SUBROUTINE WRITEFILES(NSPC,NCP1,RHO,U,T,YZ)
C      
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C
      DIMENSION RHO(*), U(*), T(*), YZ(NSPC,*)
ccc
      do j=1, ncp1
         write(2,100) rho(j),u(j),t(j),yz(1,j),yz(2,j),
     &                 yz(3,j),yz(4,j),yz(5,j),yz(6,j),
     &                 yz(7,j),yz(8,j),yz(9,j),yz(10,j)
  100 format(F,F,F,F,F,F,F,F,F,F,F,F,F)
      enddo
ccc
      return
      end
c
c -------------------------------------------------------
c
      SUBROUTINE WRITEPV(ncp1,t)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C
      DIMENSION T(*)
ccc
      CALL XRecord(3, NCP1, T)
ccc
      RETURN
C
      END
c
c -------------------------------------------------------
c
      SUBROUTINE WRITEICFILES(NSPC,NCP1,X,RHO,U,T,YZ)
C      
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C
      DIMENSION RHO(*), U(*), T(*), YZ(NSPC,*), X(*)
c
      open(4, file="ic.csv", status="unknown")
ccc
      do j=1, ncp1
         write(4,100) x(j), rho(j), u(j), t(j),
     &                yz(1,j), yz(2,j), yz(3,j), yz(4,j), yz(5,j),
     &                yz(6,j), yz(7,j), yz(8,j), yz(9,j), yz(10,j)
  100 format(F,F,F,F,F,F,F,F,F,F,F,F,F,F)
      enddo
ccc close files
      close(4)
ccc
      return
      end
c
c -------------------------------------------------------
c
      SUBROUTINE SETIC(NSPC,NCP1,XMDOT,PARAMS,P,
     1           RHO,U,T,YZ,ICKWK,RCKWK)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C
      DIMENSION PARAMS(*),RHO(*),U(*),T(*),YZ(NSPC,*),
     1          ICKWK(*),RCKWK(*),Y(NSPC) 
CCC SET BOUNDARY Y
      YSUM = 0.d0
      DO J = 1, (NSPC-1)
         Y(J) = PARAMS(19+J)
         YSUM = YSUM + Y(J)
      ENDDO
      Y(NSPC) = 1.d0 - YSUM
ccc calculate gas constant using chemkin
      call CKMMWY(Y,ickwk,rckwk,WT)!weight in g/mole
      R = 8.314462d7/WT !ergs/mol-K
ccc set boundary conditions
      U(1)    = params(6)  !cm/s
      T(1)    = params(7)  !K
      P       = params(5)  !dyne/cm2
      RHO(1)  = P/(R*T(1)) !g/cm3
      XMDOT    = RHO(1)*U(1)!g/s (if A = 1cm2)
      DO J = 1, NSPC
         YZ(J,1)   = Y(J)
      ENDDO
ccc initialize rest of cells 
      DO J = 2, NCP1
         T(j)   = T(1)         
         U(j)   = U(1)
         RHO(j) = RHO(1) 
         DO K = 1, NSPC
            YZ(K,J)  = YZ(K,1)
         ENDDO
      ENDDO
ccc
      RETURN
C
      END
c
c -------------------------------------------------------------
c
      SUBROUTINE XRecord(ifile,ncp1,s)
C 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C
      DIMENSION S(*)
C
      DO j=1, NCP1
         WRITE(ifile,*) S(j)
      ENDDO
C  
      CALL FLUSH(ifile)
C
      RETURN
C
      END
c
c ---------------------------------------------------------------
c
      SUBROUTINE READIC(NSPC,NCP1,DX,XMDOT,X,RHO,U,T,YZ)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C
      DIMENSION X(*),RHO(*),U(*),T(*),YZ(NSPC,*)
ccc OPEN/READ INITIAL CONDITION FILE
      OPEN(4, file="ic.csv",status="unknown")
      DO j = 1, NCP1
         READ(4,*) X(j), RHO(j), U(j), T(j),
     &             YZ(1,j), YZ(2,j), YZ(3,j), YZ(4,j), YZ(5,j),
     &             YZ(6,j), YZ(7,j), YZ(8,j), YZ(9,j), YZ(10,j)
      ENDDO
CCC CHECK THAT DX MATCH
      IF( (DX - (X(2) - X(1))) .GT. 0.000000000001 ) THEN
         PRINT *, 'The dxs from input.dat file and x.ic do not match'
         PRINT *, 'check dx before continuing'
         STOP
      ENDIF
ccc close files
      CLOSE(4)
ccc set XMDOT from rho(1) and u(1)
      XMDOT = RHO(1)*U(1)
ccc
      RETURN
C
      END
c
c ---------------------------------------------------------------
c
      SUBROUTINE READINPUT(ifile,PARAMS)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C
      DIMENSION PARAMS(*) 
C     
      LOGICAL LREAD
C 
      LREAD = .TRUE.
C
      J = 1
      READ(ifile, *) 
      DO WHILE( LREAD )
         READ(ifile, *)
         READ(ifile, *) PARAMS(j)
         IF( PARAMS(j) .eq. 9999.d0) LREAD = .FALSE.
         J = J + 1
      ENDDO
c
      RETURN
C
      END
c
c --------------------------------------------------------      
c
      SUBROUTINE init_random_seed()
      INTEGER :: i, n, clock
      INTEGER, DIMENSION(:), ALLOCATABLE :: seed
          
      CALL RANDOM_SEED(size = n)
      ALLOCATE(seed(n))
         
      CALL SYSTEM_CLOCK(COUNT=clock)
          
      seed = clock + 37 * (/ (i - 1, i = 1, n) /)
      CALL RANDOM_SEED(PUT = seed)
       
      DEALLOCATE(seed)
      END SUBROUTINE
c
c ------------------------------------------------------
c
      SUBROUTINE EDDYLENGTH(N,dx,PDFA,PDFB,NLENGTH)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C
c RANDOM NUMBER FROM 0-1
      CALL random_number(r)
c MAKE SURE EDDY IS .GE. 6 CELLS LONG
      NSIZE = dint(  (((r-PDFA)/PDFB)**-0.6d0) /dx  )
      DO WHILE(NSIZE .le. 5)
         NSIZE = dint(  (((r-PDFA)/PDFB)**-0.6d0) /dx  )
      ENDDO
c MAKE SURE EDDY IS DIVISIBLE BY 3
      IF(     mod(NSIZE,3) .eq. 0) THEN
         NLENGTH = NSIZE
      ELSEIF( mod(NSIZE,3) .eq. 1) THEN
         NLENGTH = NSIZE-1
      ELSEIF( mod(NSIZE,3) .eq. 2) THEN
         NLENGTH = NSIZE+1
      ENDIF
c end subroutine
      END
c
c ----------------------------------------------------------------
c
      SUBROUTINE TM(NSPC,NC,NCP1,DX,PDFA,PDFB,T,YZ,MTS)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C
      DIMENSION T(*),YZ(NSPC,*),YY(NCP1)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C MULTIPLE MAPPINGS PER TIMESTEP
C 2013 DEC 10
C
      DO I = 1, MTS
C NEW MODULE 2013 APRIL 02
C
C GET POSITION(r) AND SIZE(L) of eddy      
         CALL RANDOM_NUMBER(r)
         M = dint(r*NC)
         DO WHILE( M .le. (NCP1/4) )  ! ******** 
            CALL RANDOM_NUMBER(r)
            M = dint(r*NC)
         ENDDO
         CALL EDDYLENGTH(NC,DX,PDFA,PDFB,L)
C CHECK EDDY DOES NOT EXCEED DOMAIN
         IF(  (L+M) .GT. NC  ) M = NC - L
         DO WHILE( (L + M) .gt. NC .OR. M .le. (NCP1/4) )!mod 2013/april/02
            CALL RANDOM_NUMBER(r)  
            M = dint(r*nc)
            CALL eddylength(nc,dx,PDFA,PDFB,L)
         ENDDO
C TEMPERATURE
         CALL BTriplet(nc,M,L,t ) 
C SPECIES MASS FRACTIONS
         DO k = 1, NSPC
            DO j = 1, NCP1
               YY(j) = YZ(k,j)
            ENDDO
            CALL BTriplet(NCP1,M,L,YY)
            DO j =2, NC
               YZ(k,j) = YY(j)
            ENDDO
            YZ(k,NCP1) = YZ(k,NC)
         ENDDO
C
      ENDDO
C
      END
c
c -------------------------------------------------------------
c
      SUBROUTINE CFUEL(NSPC,NCP1,NFL,DX,XMDOT,RHO,U,YZ)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C
      DIMENSION YZ(NSPC,*) 
CCC SET REGION TO STABILIZE FLAME
      NSTAB = NCP1/3
CCC MILD CORRECTION
      IF     (  YZ(NFL,NSTAB) .lt. 0.8d0*YZ(NFL,1) .AND.
     1          U         .lt. 30.d0            ) THEN
         U    = U + 0.001d0
      ELSEIF (  YZ(NFL,NSTAB+1) .gt. 0.8d0*YZ(NFL,1) .AND.
     1          U         .gt. 15.d0            ) THEN
         U    = U - 0.0010
      ENDIF
CCC AGRESSIVE CORRECTION
      IF     (  YZ(NFL,NSTAB-2) .lt. 0.8d0*YZ(NFL,1) .AND.
     1          U         .lt. 30.d0            ) THEN
         U    = U + 0.01d0
      ELSEIF (  YZ(NFL,NSTAB+3) .gt. 0.8d0*yz(NFL,1)  .AND.
     1          U         .gt. 15.d0             ) THEN
         U    = U - 0.01d0
      ENDIF
CCC MILD CORRECTION
c      IF     (  YZ(NFL,100) .lt. 0.8d0*YZ(NFL,1) .AND.
c     1          U           .lt. 30.d0            ) THEN
c         U    = U + 0.05d0
c      ELSEIF (  YZ(NFL,400) .gt. 0.8d0*YZ(NFL,1) .AND.
c     1          U           .gt. 10.d0             ) THEN
c         U    = U - 0.05d0
c      ENDIF
C 
      XMDOT = RHO*U
C
      END
c
c -------------------------------------------------------------
c
      SUBROUTINE SET_VELOCITY(NSPC,NCP1,NFL,DX,XMDOT,GFAC,P,
     1           ICKWRK,RCKWRK,RHO,U,XMWT,T,YZ)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C
      DIMENSION T(*),YZ(NSPC,*),XMWT(NSPC),
     1 ICKWRK(*),RCKWRK(*),WDOT(NSPC)
CCC CHECK DOMAIN FUEL CONSUMPTION
      CON = 0.d0
      DO J = 2, NCP1
         CALL JYCWYP (P,T(J),YZ(1,J),ICKWRK,RCKWRK,WDOT)
         CON = CON + WDOT(NFL)
      ENDDO
      CON = -CON*GFAC*DX*XMWT(NFL)
      UF = CON/(RHO*YZ(NFL,1))
      PRINT *, 'V_IN:', UF, 'cm/s'
CCC IF PARAMETERS MATCH, CHANGE THE DOMAIN VELOCITY
      IF(  YZ(NFL,299) .lt. 0.8d0*YZ(NFL,1) .OR.
     1     YZ(NFL,302) .gt. 0.8d0*YZ(NFL,1)      ) THEN
         RETURN
      ELSE
         U = 0.5d0*(UF+U)
         XMDOT = RHO*U
         PRINT *, 'VELOCITY MODIFIED'
      ENDIF
C
      END
C
c -------------------------------------------------------------
c
      SUBROUTINE CTEMP(nspc,nspcm1,ickwk,rckwk,params,P,
     1           RHO,U,T,XMDOT)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
ccc
      DIMENSION PARAMS(*),Y(nspc),ickwk(*),rckwk(*)
C
      IF( T .gt. 294.1d0 ) THEN
C GET BOUNDARY DENSITY
         YSUM = 0.d0 
         DO j = 1, NSPCM1
            Y(j) = PARAMS(19+j)
            YSUM = YSUM + Y(j)
         ENDDO
         Y(NSPC) = 1.d0 - YSUM   
C GET GAS WEIGHT
         CALL CKMMWY(Y,ickwk,rckwk,WT)!weight in g/mole
C SET BOUNDARY TEMPERATURE
         T = T - 0.5d0!K
         RHO = (P*WT)/(8.314462d7*T)  
         XMDOT = U*RHO !g/s
      ENDIF
C
      END
C
C -------------------------------------------------------------
C
      SUBROUTINE SETLOGICALS(PARAMS,LRIC,LTM,LWI,LOWIC,LCF,LCT)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C
      LOGICAL LRIC,LTM,LWI,LOWIC,LCF,LCT
C     
      DIMENSION PARAMS(*)
C
C     READ INITIAL CONDITIONS
      IF (params(12) .ne. 0) THEN
         LRIC = .TRUE.
      ELSE 
         LRIC = .FALSE.
      ENDIF
C     DO TRIPLET MAP
      IF (params(13) .ne. 0) THEN
         LTM = .TRUE.
      ELSE 
         LTM = .FALSE.
      ENDIF
C     WRITE DATA IN 1/10 SIMULATION INTERVALS
      IF (params(14) .ne. 0) THEN
         LWI = .TRUE.
      ELSE 
         LWI = .FALSE.
      ENDIF
C     OVERWRITE INITIAL CONDITIONS
      IF (params(15) .ne. 0) THEN
         LOWIC = .TRUE.
      ELSE 
         LOWIC = .FALSE.
      ENDIF
C     DO FUEL CORRECTION
      IF (params(16) .ne. 0) THEN
         LCF = .TRUE.
      ELSE 
         LCF = .FALSE.
      ENDIF
C     DO TEMPERATURE CORRECTION
      IF (params(17) .ne. 0) THEN
         LCT = .TRUE.
      ELSE 
         LCT = .FALSE.
      ENDIF
C
      END
C
C --------------------------------------------------
C
      SUBROUTINE INIT_ALL(NSPC,NCP1,NCM1,NFL,PARAMS,GFAC,
     1 LRIC,LTM,LWI,LOWIC,LCF,LCT,
     2 XMWT,RHO,U,T,YZ,XMDOT,P,
     3 lickwk,lrckwk,lcckwk,ickwk,rckwk,cckwk,
     4 lenimc,lenrmc,imcwrk,rmcwrk,
     5 TEND,DT,NTS,XMDT,DOM,DX,X,
     6 NSIM,NTSPSIM)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C
      LOGICAL LRIC,LTM,LWI,LOWIC,LCF,LCT
C     
      DIMENSION PARAMS(*),XMWT(*),ICKWK(*),RCKWK(*),CCKWK(*),
     1          IMCWRK(*),RMCWRK(*),
     2          X(*),RHO(*),U(*),T(*),YZ(NSPC,*)
ccc read input
      OPEN(99,file='input7.dat',status='old')
      CALL READINPUT(99,params)
ccc flags
      CALL SETLOGICALS(PARAMS,LRIC,LTM,LWI,LOWIC,LCF,LCT)
ccc retrieve molecular weights
      CALL CKWT(ickwk, rckwk, XMWT)
         OPEN(930,file='mwtht.d',form='formatted',status='unknown')
      DO i = 1, nspc
         WRITE(930,*) XMWT(i)
      ENDDO
      CLOSE(930)
ccc time parameters
      TEND = params(1)
      DT   = params(2)
      NTS  = tend/dt !number
ccc domain length parameters
      DOM = params(4)     !meters
      DX  = DOM/dfloat(NCM1) !meters
ccc lumping dx into dt, calling it mdt (modified dt)
      XMDT = DT/DX
ccc x coordinates
      DO j = 1, ncp1
         X(j) = (j-1)*DX
      ENDDO
CCC WDOT MULTIPLIER (SCALES REACTION SPEEDS)
      GFAC = PARAMS(18)
CCC FUEL ARRAY LOCATION
      NFL = PARAMS(19)
CCC NUMBER OF REALIZATIONS TO RECORD
      NSIM = PARAMS(31)
CCC NUMBER OF TIMESTEPS PER REALIZATION
      NTSPSIM = PARAMS(32)
ccc data files
      OPEN(1,file="x.csv", status="unknown")!spacing data
ccc write down the x spacing once
      CALL XRecord(1,NCP1,X)
ccc write input to z.ic as record for checking parameters later
      WRITE(9,*) 'input7.dat parameters'
      DO J = 1, 32
         WRITE(9,*) params(j)
      ENDDO
CCC
      IF(LRIC) then !read initial conditions from files OR
         CALL READIC(NSPC,NCP1,DX,XMDOT,X,RHO,U,T,YZ)
         P = params(5) !domain pressure set at input file
      ELSE !set initial conditions according to input.dat
         CALL SETIC(NSPC,NCP1,XMDOT,PARAMS,P,RHO,U,T,YZ,ICKWK,RCKWK)

      ENDIF
ccc
      END   
C
C -----------------------------------------------------------------------
C
      SUBROUTINE INIT_LEM(NC,PARAMS,DOM,TEND,DT,U,
     1           NTS_COUNT,NTS_PE,PDFA,PDFB,MTS)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C
      DIMENSION PARAMS(*),U(*)
C
      NTS_COUNT = 0
      XNU = params(9)!air at 300K in cm2/s
      Re  = params(3)!turbulent Reynolds 
      XLint = params(30)!Integral length
      XLk   = XLint/(Re**0.75) !Kolmogorov length
      C_lambda = 15.d0
c      Rate = DOM*(4.8*XNU/XLint**3)*( (XLint/XLk)**(4/3) )*
c     &       ( (XLint/XLk)**(5/3) - 1)/( 1 - (XLk/XLint)**(4/3) )
      Rate = DOM*(54.d0/5.d0)*
     & ( XNU*Re / (C_lambda*XLint**3) )*
     & ( (XLint/XLk)**(5/3) - 1)/( 1 - (XLk/XLint)**(4/3) )
      !eddies per second
      PDFA = (XLint**(5.d0/3.d0)) * (XLk**(-5.d0/3.d0)) / 
     &       ( (XLint/XLk)**(5.d0/3.d0) -1.d0 )
      PDFB = -(XLint**(5.d0/3.d0))/(  (XLint/XLk)**(5.d0/3.d0) -1.d0 )
      TAU = XLk*XLk/XNU !komogorov time
ccc calculate number of timesteps per eddy
      NTS_PE = (1.d0/Rate)/DT+1
CCC if (1/Rate)/DT is less than 1, need multiple triplet maps per timestep
      IF ( (1.d0/Rate)/DT < 1.d0 ) THEN
         MTS = Rate*DT+1
      ELSE
         MTS = 1
      ENDIF
      PRINT *, "Mappings per TimeStep: ", MTS
C
      write(9,*) 
      write(9,*) 'Number of cells: ', nc
      write(9,*) 'Domain length: ', dom, ' cm'
      write(9,*) 'Number of timesteps per eddy: ', NTS_PE
      write(9,*) 'Turbulent Reynolds: ', Re
      write(9,*) 'Boundary velocity: ', u(1)
      write(9,*) 'Sim time and timestep size: ', tend, dt
      write(9,*), 'PDFA: ', PDFA
      write(9,*), 'PDFB: ', PDFB
      write(9,*), 'Integral scale: ', XLint, ' cm'
      write(9,*), 'Kolmogorov scale: ', XLk, ' cm'
C
      END
C
C ----------------------------------------------------------------
C
      SUBROUTINE GET_RHO_U(NC,NCP1,NSPC,XMDOT,
     1                     T,YZ,P,ICKWK,RCKWK,RHO,U)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C
      DIMENSION RHO(*),U(*),T(*),ICKWK(*),RCKWK(*),YZ(NSPC,*)
C
      DO j = 2, NC
         call CKMMWY(yz(1,j),ickwk,rckwk,WT)!weight in g/mole
         RHO(j) = P*WT/(8.314462d7*T(j))
         U(j)   = XMDOT/RHO(j) 
      ENDDO 
      RHO(NCP1) = RHO(NC)
      U(NCP1)   = U(NC)
C
      END
C
C ----------------------------------------------------------------
C
      SUBROUTINE CHECKPOINT(XMDOT,I,NTS,DT,
     1 T,Y,P,ICKWK,RCKWK)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C
      DIMENSION ICKWK(*),RCKWK(*),Y(*)
C     CHECK FOR ERRORS AT LAST CELL
      CALL CKMMWY(Y,ickwk,rckwk,WT)!weight in g/mole
      U   = XMDOT*(8.314462d7*T)/(P*WT) 
C
      IF( .not. (U .le. 2.d4) .OR. (U .lt. 0.d0) ) THEN
C         PRINT *, 'Time of Death: ', dfloat((k-1)*nts/10+i)*dt
C         WRITE(9,*) 'Time of Death: ', dfloat((k-1)*nts/10+i)*dt
c         GO TO 79
      ENDIF
C 
      END  

           

