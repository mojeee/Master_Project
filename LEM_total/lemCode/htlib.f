C
C----------------------------------------------------------------------C
C
      SUBROUTINE HTW_L1  (P, T, Y, ICKWRK, RCKWRK, WDOT)
C
C----------------------------------------------------------------------C
C
C  This is from a 1-step mechanism
C
C  INPUT
C
C     P      - Pressure.
C                   cgs units - dynes/cm**2
C                   Data type - real scalar
C     T      - Temperature.
C                   cgs units - K
C                   Data type - real scalar
C     Y      - Mass fractions of the species.
C                   cgs units - none
C                   Data type - real array
C                   Dimension Y(*) at least KK, the total number of
C                   species.
C     ICKWRK - Array of integer workspace.
C                   Data type - integer array
C                   Dimension ICKWRK(*) at least LENIWK.
C     RCKWRK - Array of real work space.
C                   Data type - real array
C                   Dimension RCKWRK(*) at least LENRWK.
C
C  OUTPUT
C     WDOT   - Chemical molar production rates of the species.
C                   cgs units - kg/(m**3*sec)
C                   Data type - real array
C                   Dimension WDOT(*) at least KK, the total number of
C                   species.
C
C  END PROLOGUE
C
C*****double precision
      IMPLICIT REAL*8 (A-H, O-Z), INTEGER (I-N)
C*****END double precision
C*****single precision
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END single precision
C
      DIMENSION ICKWRK(*), RCKWRK(*), Y(*), WDOT(*)
      COMMON /CKSTRT/ NMM , NKK , NII , MXSP, MXTB, MXTP, NCP , NCP1,
     1                NCP2, NCP2T,NPAR, NLAR, NFAR, NLAN, NFAL, NREV,
     2                NTHB, NRLT, NWL,  IcMM, IcKK, IcNC, IcPH, IcCH,
     3                IcNT, IcNU, IcNK, IcNS, IcNR, IcLT, IcRL, IcRV,
     4                IcWL, IcFL, IcFO, IcKF, IcTB, IcKN, IcKT, NcAW,
     5                NcWT, NcTT, NcAA, NcCO, NcRV, NcLT, NcRL, NcFL,
     6                NcKT, NcWL, NcRU, NcRC, NcPA, NcK1, NcK2, NcK3,
     7                NcK4, NcI1, NcI2, NcI3, NcI4
C
C CKYTX = mass fraction to mole fraction
C
      CALL CKYTX (Y, ICKWRK, RCKWRK, RCKWRK(NcK1))
      DO 50 K=1,NKK
         WDOT(K) = 0.d0
   50 CONTINUE
C
C         SET LOCAL VALUES FOR THE CONCENTRATIONS
C
      Xch4 = RCKWRK(NcK1)
c      Xco2 = RCKWRK(NcK1+1)
c      Xh2o = RCKWRK(NcK1+2)
      Xo2  = RCKWRK(NcK1+3)
C
C         SET NEEDED LOCAL VARIABLES
C
      R=RCKWRK(NcRU)
      PBYRT2=(P/(R*T))**2.0 !units (mols/cc)^2
C
C         SET REACTION RATES
C
c      RK1 = 6.91d14*EXP(-1.590934d4/T)!original
c      RK1 = 4.d14*EXP(-1.590934d4/T)!edit1
      RK1 = 1.47d14*EXP(-1.6d4/T)!edit2 2013/02/08,22:00
c      RK1 = 2.119d11*EXP(-1.4906d4/T)
      RATE = RK1*Xch4*Xo2*PBYRT2
C         DETERMINE RATES OF COMBINED REACTANT/PRODUCT COMBINATIONS
C all in units of mol/cc-s
c      DXfDT = - RATE
c      DXoDT = - 2.d0*RATE 
c      DXiDT =   RATE
c      DXpDT =   2.d0*RATE
C         SET WDOT(I) FROM COMBINED REACTANT/PRODUCT RATES
c convert from mol/cc-s to kg/m3-s
c 1.d3 needed for unit conversion
      WDOT(1) = -     RATE*RCKWRK(NcWT  )*1.d3 
      WDOT(2) =       RATE*RCKWRK(NcWT+1)*1.d3
      WDOT(3) =   2.0*RATE*RCKWRK(NcWT+2)*1.d3
      WDOT(4) = - 2.0*RATE*RCKWRK(NcWT+3)*1.d3
C
      RETURN
      END
C
C----------------------------------------------------------------------C
C
      SUBROUTINE HTW_2dns(P, T, Y, ICKWRK, RCKWRK, WDOT)
C
C----------------------------------------------------------------------C
C
C     THIS DNS MECHANISM HAS BEEN MODIFIED TO FIT CAMBRIDGE 
C     EXPERIMENTS BETTER BY CHECKING/TESTING FLAME PROFILE 
C     FROM CANTERA FIRST!
C
C
C
C  INPUT
C
C     P      - Pressure.
C                   cgs units - dynes/cm**2
C                   Data type - real scalar
C     T      - Temperature.
C                   cgs units - K
C                   Data type - real scalar
C     Y      - Mass fractions of the species.
C                   cgs units - none
C                   Data type - real array
C                   Dimension Y(*) at least KK, the total number of
C                   species.
C     ICKWRK - Array of integer workspace.
C                   Data type - integer array
C                   Dimension ICKWRK(*) at least LENIWK.
C     RCKWRK - Array of real work space.
C                   Data type - real array
C                   Dimension RCKWRK(*) at least LENRWK.
C
C  OUTPUT

C     WDOT   - Chemical molar production rates of the species.
C                   cgs units - kg/(m**3*sec)
C                   Data type - real array
C                   Dimension WDOT(*) at least KK, the total number of
C                   species.
C     WR1    - rxn 1 rate
c                   units     - mols/cc/s
c                   data type - double precision
c     WR2    - rxn 2 rate
c                   units     - mols/cc/s
c                   data type  - double precision                 
C
C  END PROLOGUE
C
C*****double precision
      IMPLICIT REAL*8 (A-H, O-Z), INTEGER (I-N)
C*****END double precision
C*****single precision
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END single precision
C
      DIMENSION ICKWRK(*), RCKWRK(*), Y(*), WDOT(*)
      COMMON /CKSTRT/ NMM , NKK , NII , MXSP, MXTB, MXTP, NCP , NCP1,
     1                NCP2, NCP2T,NPAR, NLAR, NFAR, NLAN, NFAL, NREV,
     2                NTHB, NRLT, NWL,  IcMM, IcKK, IcNC, IcPH, IcCH,
     3                IcNT, IcNU, IcNK, IcNS, IcNR, IcLT, IcRL, IcRV,
     4                IcWL, IcFL, IcFO, IcKF, IcTB, IcKN, IcKT, NcAW,
     5                NcWT, NcTT, NcAA, NcCO, NcRV, NcLT, NcRL, NcFL,
     6                NcKT, NcWL, NcRU, NcRC, NcPA, NcK1, NcK2, NcK3,
     7                NcK4, NcI1, NcI2, NcI3, NcI4
C
C CKYTX = mass fraction to mole fraction
C
      CALL CKYTX (Y, ICKWRK, RCKWRK, RCKWRK(NcK1))
      DO 50 K=1,NKK
         WDOT(K) = 0.d0
   50 CONTINUE
C
C         SET LOCAL VALUES FOR THE CONCENTRATIONS
C
      Xf=RCKWRK(NcK1)
      Xo=RCKWRK(NcK1+1)
      Xi=RCKWRK(NcK1+2)
      Xp=RCKWRK(NcK1+3)
C
C         SET NEEDED LOCAL VARIABLES
C
      R=RCKWRK(NcRU)
      PBYRT=P/(R*T) !units mols/cc
      PBYRT3 = PBYRT**3 !will need this later
C
C         SET REACTION RATES
C         ORIGINAL
c      RK1=7.8d14*EXP(-7640.d0/T)
c      RK2=8.9d16*(T**-1.8)     
c      RK3=4.5d5*EXP(2300.d0/T)*(Xo**0.5d0)*(Xi**1.5d0)*
c     & 0.379748595d0*EXP(Xf/Xo)
c     The constant 0.3797.. is obtained from exp(-0.5*sqrt(15/4))
c
c         MODIFIED - working 02/23/22:51
c
c      RK1=13.d14*EXP(-8500.d0/T) 
c      RK2=8.5d16*(T**-1.78)   
c      RK3=4.5d5*EXP(2300.d0/T)*(Xo**0.5d0)*(Xi**1.5d0)*
c     & 0.017d0*EXP(Xf/Xo) !0.02
c
c         MODIFING - working 02/28/20:29 (closest to cantera phi = 0.61 values)
c
c      RK1=8.d14*EXP(-7640.d0/T)
c      RK2=8.5d16*(T**-1.78)   
c      RK3=4.5d5*EXP(2300.d0/T)*(Xo**0.5d0)*(Xi**1.5d0)*
c     & 0.017*EXP(Xf/Xo)
c
c         MODIFING - working 03/07/12:04 (closest to cantera phi = 1.0 values)

c      RK1=8.d14*EXP(-7640.d0/T)
c      RK2=20.d16*(T**-1.6)   
c      RK3=4.5d5*EXP(2300.d0/T)*(Xo**0.5d0)*(Xi**1.5d0)*
c     & 0.02*EXP(Xf/Xo)
c
c         MODIFING - working 03/14/02:46 (phi = 0.73)
c
      RK1=8.d14*EXP(-7640.d0/T)
      RK2=10.5d16*(T**-1.65)   
      RK3=4.5d5*EXP(2300.d0/T)*(Xo**0.5d0)*(Xi**1.5d0)*
     & 0.04*EXP(Xf/Xo)
ccc print
ccc print
C
C         DETERMINE RATES OF COMBINED REACTANT/PRODUCT COMBINATIONS
      RATE1 = RK1*Xf*RK3
      RATE2 = RK2*Xo*RK3
C all in units of mol/cc-s
      DXfDT = (-RATE1              )*PBYRT3
      DXoDT = (-RATE1 -      RATE2 )*PBYRT3
      DXiDT = ( RATE1 -      RATE2 )*PBYRT3
      DXpDT = ( RATE1 + 2.d0*RATE2 )*PBYRT3

C         SET WDOT(I) FROM COMBINED REACTANT/PRODUCT RATES
c convert from mol/cc-s to kg/m3-s
c 1.d3 needed for unit conversion
      WDOT(1) = DXfDT*RCKWRK(NcWT  )*1.d3 
      WDOT(2) = DXoDT*RCKWRK(NcWT+1)*1.d3
      WDOT(3) = DXiDT*RCKWRK(NcWT+2)*1.d3
      WDOT(4) = DXpDT*RCKWRK(NcWT+3)*1.d3
C
      RETURN
      END

C
C----------------------------------------------------------------------C
C
      SUBROUTINE HTW_CM2(P, T, Y, ICKWRK, RCKWRK, WDOT)
C
C----------------------------------------------------------------------C
C
C     standard 2S-CM2 mechanism
C
C
C
C  INPUT
C
C     P      - Pressure.
C                   cgs units - dynes/cm**2
C                   Data type - real scalar
C     T      - Temperature.
C                   cgs units - K
C                   Data type - real scalar
C     Y      - Mass fractions of the species.
C                   cgs units - none
C                   Data type - real array
C                   Dimension Y(*) at least KK, the total number of
C                   species.
C     ICKWRK - Array of integer workspace.
C                   Data type - integer array
C                   Dimension ICKWRK(*) at least LENIWK.
C     RCKWRK - Array of real work space.
C                   Data type - real array
C                   Dimension RCKWRK(*) at least LENRWK.
C
C  OUTPUT

C     WDOT   - Chemical molar production rates of the species.
C                   cgs units - kg/(m**3*sec)
C                   Data type - real array
C                   Dimension WDOT(*) at least KK, the total number of
C                   species.
C     WR1    - rxn 1 rate
c                   units     - mols/cc/s
c                   data type - double precision
c     WR2    - rxn 2 rate
c                   units     - mols/cc/s
c                   data type  - double precision                 
C
C  END PROLOGUE
C
C*****double precision
      IMPLICIT REAL*8 (A-H, O-Z), INTEGER (I-N)
C*****END double precision
C*****single precision
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END single precision
C
      DIMENSION ICKWRK(*), RCKWRK(*), Y(*), WDOT(*)
      COMMON /CKSTRT/ NMM , NKK , NII , MXSP, MXTB, MXTP, NCP , NCP1,
     1                NCP2, NCP2T,NPAR, NLAR, NFAR, NLAN, NFAL, NREV,
     2                NTHB, NRLT, NWL,  IcMM, IcKK, IcNC, IcPH, IcCH,
     3                IcNT, IcNU, IcNK, IcNS, IcNR, IcLT, IcRL, IcRV,
     4                IcWL, IcFL, IcFO, IcKF, IcTB, IcKN, IcKT, NcAW,
     5                NcWT, NcTT, NcAA, NcCO, NcRV, NcLT, NcRL, NcFL,
     6                NcKT, NcWL, NcRU, NcRC, NcPA, NcK1, NcK2, NcK3,
     7                NcK4, NcI1, NcI2, NcI3, NcI4
C
C CKYTX = mass fraction to mole fraction
C
      CALL CKYTX (Y, ICKWRK, RCKWRK, RCKWRK(NcK1))
      DO 50 K=1,NKK
         WDOT(K) = 0.d0
   50 CONTINUE
C
C         SET LOCAL VALUES FOR THE CONCENTRATIONS
C
      Xch4= RCKWRK(NcK1  )
      Xco = RCKWRK(NcK1+1)
      Xco2= RCKWRK(NcK1+2)
      Xh2 = RCKWRK(NcK1+3)
      Xo2 = RCKWRK(NcK1+4)
C
C         SET NEEDED LOCAL VARIABLES
C
      R=RCKWRK(NcRU)
      PBYRT=P/(R*T) !units mols/cc
c      write(8,*) PBYRT
      PBYRT2 = PBYRT**2 !will need this later
C
C         SET REACTION RATES
C
      RK1=2.0000d15*EXP(-17624.45d0/T)
      RK2=8.1104d10*EXP(-38871.48d0/T)   
      RK3=2.0000d9 *EXP(-6042.669d0/T)
ccc print
ccc print
C
C         DETERMINE RATES OF COMBINED REACTANT/PRODUCT COMBINATIONS
      RATE1 = RK1*(Xch4**0.9)*(Xo2**1.1)*PBYRT2 !the PBYRT accounts for 
      RATE2 = RK2*(Xco2)*PBYRT                  !the mol/cc of each 
      RATE3 = RK3*(Xco)*(Xo2**0.5) *PBYRT**1.5  !species mole fraction
C all in units of mol/cc-s
      DXch4DT = -RATE1
      DXcoDT  =  RATE1 + RATE2 - RATE3
      DXco2DT = -RATE2 + RATE3
      DXh2DT  = 2.0*RATE1 
      DXo2DT  = ( -RATE1 + RATE2 - RATE3 ) *0.5

C         SET WDOT(I) FROM COMBINED REACTANT/PRODUCT RATES
c convert from mol/cc-s to kg/m3-s
c 1.d3 needed for unit conversion
      WDOT(1) = DXch4DT*RCKWRK(NcWT  )*1.d3 
      WDOT(2) = DXcoDT *RCKWRK(NcWT+1)*1.d3
      WDOT(3) = DXco2DT*RCKWRK(NcWT+2)*1.d3
      WDOT(4) = DXh2DT *RCKWRK(NcWT+3)*1.d3
      WDOT(5) = DXo2DT *RCKWRK(NcWT+4)*1.d3
c      write(8,*) wdot(1), wdot(3), wdot(5)
c      write(8,*) RATE1, RATE2, RATE3
C
      RETURN
      END
C
C----------------------------------------------------------------------C
C
      SUBROUTINE HTW_WD(P, T, Y, ICKWRK, RCKWRK, WDOT)
C
C----------------------------------------------------------------------C
C
C     standard WD mechanism
C
C
C
C  INPUT
C
C     P      - Pressure.
C                   cgs units - dynes/cm**2
C                   Data type - real scalar
C     T      - Temperature.
C                   cgs units - K
C                   Data type - real scalar
C     Y      - Mass fractions of the species.
C                   cgs units - none
C                   Data type - real array
C                   Dimension Y(*) at least KK, the total number of
C                   species.
C     ICKWRK - Array of integer workspace.
C                   Data type - integer array
C                   Dimension ICKWRK(*) at least LENIWK.
C     RCKWRK - Array of real work space.
C                   Data type - real array
C                   Dimension RCKWRK(*) at least LENRWK.
C
C  OUTPUT

C     WDOT   - Chemical molar production rates of the species.
C                   cgs units - kg/(m**3*sec)
C                   Data type - real array
C                   Dimension WDOT(*) at least KK, the total number of
C                   species.
C     WR1    - rxn 1 rate
c                   units     - mols/cc/s
c                   data type - double precision
c     WR2    - rxn 2 rate
c                   units     - mols/cc/s
c                   data type  - double precision                 
C
C  END PROLOGUE
C
C*****double precision
      IMPLICIT REAL*8 (A-H, O-Z), INTEGER (I-N)
C*****END double precision
C*****single precision
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END single precision
C
      DIMENSION ICKWRK(*), RCKWRK(*), Y(*), WDOT(*)
      COMMON /CKSTRT/ NMM , NKK , NII , MXSP, MXTB, MXTP, NCP , NCP1,
     1                NCP2, NCP2T,NPAR, NLAR, NFAR, NLAN, NFAL, NREV,
     2                NTHB, NRLT, NWL,  IcMM, IcKK, IcNC, IcPH, IcCH,
     3                IcNT, IcNU, IcNK, IcNS, IcNR, IcLT, IcRL, IcRV,
     4                IcWL, IcFL, IcFO, IcKF, IcTB, IcKN, IcKT, NcAW,
     5                NcWT, NcTT, NcAA, NcCO, NcRV, NcLT, NcRL, NcFL,
     6                NcKT, NcWL, NcRU, NcRC, NcPA, NcK1, NcK2, NcK3,
     7                NcK4, NcI1, NcI2, NcI3, NcI4
C
C CKYTX = mass fraction to mole fraction
C
      CALL CKYTX (Y, ICKWRK, RCKWRK, RCKWRK(NcK1))
      DO 50 K=1,NKK
         WDOT(K) = 0.d0
   50 CONTINUE
C
C         SET LOCAL VALUES FOR THE CONCENTRATIONS
C
      Xch4= RCKWRK(NcK1  )
      Xco = RCKWRK(NcK1+1)
      Xco2= RCKWRK(NcK1+2)
      Xh2o= RCKWRK(NcK1+3)
      Xo2 = RCKWRK(NcK1+4)
C
C         SET NEEDED LOCAL VARIABLES
C
      R=RCKWRK(NcRU)
      PBYRT=P/(R*T) !units mols/cc
c      write(8,*) PBYRT
      PBYRT2 = PBYRT**2 !will need this later
C
C         SET REACTION RATES
C
      RK1=5.00d11*EXP(-24069.96d0/T)
      RK2=2.24d12*EXP(-20494.71d0/T)   
      RK3=5.00d8 *EXP(-20494.71d0/T)
ccc print
ccc print
C
C         DETERMINE RATES OF COMBINED REACTANT/PRODUCT COMBINATIONS
      RATE1 = RK1*(Xch4**0.7)*(Xo2**0.8)*PBYRT**1.5 !the PBYRT accounts for 
      RATE2 = RK2*(Xco)*(Xh2o)*PBYRT2               !the mol/cc of each 
      RATE3 = RK3*(Xco2)*PBYRT  !species mole fraction
C all in units of mol/cc-s
      DXch4DT = -RATE1
      DXcoDT  =  RATE1 - RATE2 + RATE3
      DXco2DT =  RATE2 - RATE3
      DXh2oDT = 2.0*RATE1 
      DXo2DT  = 0.5*( - 3.0*RATE1 - RATE2 + RATE3 )

C         SET WDOT(I) FROM COMBINED REACTANT/PRODUCT RATES
c convert from mol/cc-s to kg/m3-s
c 1.d3 needed for unit conversion
      WDOT(1) = DXch4DT*RCKWRK(NcWT  )*1.d3 
      WDOT(2) = DXcoDT *RCKWRK(NcWT+1)*1.d3
      WDOT(3) = DXco2DT*RCKWRK(NcWT+2)*1.d3
      WDOT(4) = DXh2oDT*RCKWRK(NcWT+3)*1.d3
      WDOT(5) = DXo2DT *RCKWRK(NcWT+4)*1.d3
c
c         write(8,*) WDOT(1), WDOT(2), WDOT(3)
c
c
c      if (WDOT(1) .lt. 0.d0 ) then
c         print *, 'wdot 1 error'
ccc      elseif (WDOT(2) .lt. 0.d0 ) then
c         print *, 'wdot 2 error'
c      elseif (WDOT(3) .lt. 0.d0 ) then
c         print *, 'wdot 3 error'
c      elseif (WDOT(4) .lt. 0.d0 ) then
c         print *, 'wdot 4 error'
c      elseif (WDOT(5) .lt. 0.d0 ) then
c         print *, 'wdot 5 error'
c      endif

c      write(8,*) wdot(1), wdot(3), wdot(5)
c      write(8,*) RATE1, RATE2, RATE3
C
      RETURN
      END


C
C----------------------------------------------------------------------C
C
      SUBROUTINE HTW_2a(P, T, Y, RHO, ICKWRK, RCKWRK, WDOT)
C
C----------------------------------------------------------------------C
C
C  INPUT
C
C     P      - Pressure.
C                   cgs units - dynes/cm**2
C                   Data type - real scalar
C     T      - Temperature.
C                   cgs units - K
C                   Data type - real scalar
C     Y      - Mass fractions of the species.
C                   cgs units - none
C                   Data type - real array
C                   Dimension Y(*) at least KK, the total number of
C                   species.
C     RHO    - Density
C                   cgs units - g/cm3
C                   Data type - real scalar
C
C
C     ICKWRK - Array of integer workspace.
C                   Data type - integer array
C                   Dimension ICKWRK(*) at least LENIWK.
C     RCKWRK - Array of real work space.
C                   Data type - real array
C                   Dimension RCKWRK(*) at least LENRWK.
C
C  OUTPUT

C     WDOT   - Chemical molar production rates of the species.
C                   cgs units - kg/(m**3*sec)
C                   Data type - real array
C                   Dimension WDOT(*) at least KK, the total number of
C                   species.
C
C  END PROLOGUE
C
C*****double precision
      IMPLICIT REAL*8 (A-H, O-Z), INTEGER (I-N)
C*****END double precision
C*****single precision
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END single precision
C
      DIMENSION ICKWRK(*), RCKWRK(*), Y(*), WDOT(*)
      COMMON /CKSTRT/ NMM , NKK , NII , MXSP, MXTB, MXTP, NCP , NCP1,
     1                NCP2, NCP2T,NPAR, NLAR, NFAR, NLAN, NFAL, NREV,

     2                NTHB, NRLT, NWL,  IcMM, IcKK, IcNC, IcPH, IcCH,
     3                IcNT, IcNU, IcNK, IcNS, IcNR, IcLT, IcRL, IcRV,
     4                IcWL, IcFL, IcFO, IcKF, IcTB, IcKN, IcKT, NcAW,
     5                NcWT, NcTT, NcAA, NcCO, NcRV, NcLT, NcRL, NcFL,
     6                NcKT, NcWL, NcRU, NcRC, NcPA, NcK1, NcK2, NcK3,
     7                NcK4, NcI1, NcI2, NcI3, NcI4
C
c      R=RCKWRK(NcRU)
c      PBYRT=P/(R*T) !units mols/cm3
C
C         SET REACTION RATES
C
C      THETA_ = EXP( (-4.26d-10*T**3.0)*Y(1)*EXP(4051.d0/T)/Y(2) )
c      GAMMA_ = 1.4d-2 + (3.92d-1)*Y(1) - (2.d-3)*Y(2) 
c     &                + (7.2d-2)*Y(3)  + (1.67d-1)*Y(4)
      HYDRA_ = (1.36d-1)*EXP(3045.d0/T)*
     &          EXP( (-4.26d-10*T**3.0)*Y(1)*EXP(4051.d0/T)/Y(2) )*
     &          ( 1.d0 - EXP( (-10.d-15)*T**5.0 ) )*
     &            ( (Y(2)*Y(3)**3.0)**0.5 )/Y(4)
C
C         units are mols/cm3/second 
C     1. extra factor of RHO in each RK expression for correct units
C     2. each pre-factor has been multiplied by 1000 so that I don't need
C     to convert from g/cm3 to kg/m3 at a later step
C - originals
c      RK1 = 1375.d3*(RHO**2.0)*(T**3.0)*EXP( -4404.d0/T )*Y(1)*HYDRA_
c      RK2 = 7.19d19*(RHO**3.0)*(T**-0.8)*Y(2)*
c     & (1.4d-2 + 3.92d-1*Y(1) - 2.d-3*Y(2) + 7.2d-2*Y(3) + 1.67d-1*Y(4))
c     & *HYDRA_     
c - changes
      RK1 = 7000.d3*(RHO**2.0)*(T**2.6)*EXP( -4600.d0/T )*Y(1)*HYDRA_
      RK2 = 6.5d19*(RHO**3.0)*(T**-0.72)*Y(2)*
     & (1.4d-2 + 3.92d-1*Y(1) - 2.d-3*Y(2) + 7.2d-2*Y(3) + 1.67d-1*Y(4))
     & *HYDRA_     
C
C         SET WDOT(I) FROM COMBINED REACTANT/PRODUCT RATES
C     convert from mol/cc-s to kg/m3-s with RCKWRK(NcWT ...)
C
      WDOT(1) = (-RK1           )*RCKWRK(NcWT  )
      WDOT(2) = (-RK1 -      RK2)*RCKWRK(NcWT+1)
      WDOT(3) = ( RK1 -      RK2)*RCKWRK(NcWT+2)
      WDOT(4) = ( RK1 + 2.d0*RK2)*RCKWRK(NcWT+3)
C
c      write(13,*) RK1, RK2, 'mol/cc/s x 1000'
c      write(13,*) RCKWRK(NcWT  ), RCKWRK(NcWT+1), 'g/mol'
c      write(13,*) RCKWRK(NcWT+2), RCKWRK(NcWT+3), 'g/mol'
c      write(13,*) WDOT(1), WDOT(2), WDOT(3)
c      write(13,*) WDOT(4), 'kg/m3/s'
c      write(13,*)
C
      RETURN
      END
C
C----------------------------------------------------------------------C
C
      SUBROUTINE HTW_2b(P, T, Y, RHO, ICKWRK, RCKWRK, WDOT)
C
C----------------------------------------------------------------------C
C
C  INPUT
C
C     P      - Pressure.
C                   cgs units - dynes/cm**2
C                   Data type - real scalar
C     T      - Temperature.
C                   cgs units - K
C                   Data type - real scalar
C     Y      - Mass fractions of the species.
C                   cgs units - none
C                   Data type - real array
C                   Dimension Y(*) at least KK, the total number of
C                   species.
C     RHO    - Density
C                   cgs units - g/cm3
C                   Data type - real scalar
C
C
C     ICKWRK - Array of integer workspace.
C                   Data type - integer array
C                   Dimension ICKWRK(*) at least LENIWK.
C     RCKWRK - Array of real work space.
C                   Data type - real array
C                   Dimension RCKWRK(*) at least LENRWK.
C
C  OUTPUT

C     WDOT   - Chemical molar production rates of the species.
C                   cgs units - kg/(m**3*sec)
C                   Data type - real array
C                   Dimension WDOT(*) at least KK, the total number of
C                   species.
C
C  END PROLOGUE
C
C*****double precision
      IMPLICIT REAL*8 (A-H, O-Z), INTEGER (I-N)
C*****END double precision
C*****single precision
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END single precision
C
      DIMENSION ICKWRK(*), RCKWRK(*), Y(*), WDOT(*)
      COMMON /CKSTRT/ NMM , NKK , NII , MXSP, MXTB, MXTP, NCP , NCP1,
     1                NCP2, NCP2T,NPAR, NLAR, NFAR, NLAN, NFAL, NREV,

     2                NTHB, NRLT, NWL,  IcMM, IcKK, IcNC, IcPH, IcCH,
     3                IcNT, IcNU, IcNK, IcNS, IcNR, IcLT, IcRL, IcRV,
     4                IcWL, IcFL, IcFO, IcKF, IcTB, IcKN, IcKT, NcAW,
     5                NcWT, NcTT, NcAA, NcCO, NcRV, NcLT, NcRL, NcFL,
     6                NcKT, NcWL, NcRU, NcRC, NcPA, NcK1, NcK2, NcK3,
     7                NcK4, NcI1, NcI2, NcI3, NcI4
C
c      R=RCKWRK(NcRU)
c      PBYRT=P/(R*T) !units mols/cm3
C
C         SET REACTION RATES
C
      THETA_ = EXP( -22.1d0*Y(1)/Y(2) )
      GAMMA_ = 1.4d-2 + (3.92d-1)*Y(1) - (2.d-3)*Y(2) 
     &                + (7.2d-2)*Y(3) + (1.67d-1)*Y(4)
      HYDRA_ = (1.36d-1)*EXP(3045.d0/T)*THETA_*
     &          ( 1.d0 - EXP( (-10.d-15)*T**5.0 ) )*
     &            ( (Y(2)*Y(3)**3.0)**0.5 )/Y(4)
C
C         units are mols/cm3/second 
C     1. extra factor of RHO in each RK expression for correct units
C     2. each pre-factor has been multiplied by 1000 so that I don't need
C     to convert from g/cm3 to kg/m3 at a later step
C
      RK1 = 1.57d17*(RHO**2.0)*EXP( -9624.d0/T )*Y(1)*HYDRA_
      RK2 = 7.19d19*(RHO**3.0)*(T**-0.8)*Y(2)*HYDRA_*GAMMA_     
C
C         SET WDOT(I) FROM COMBINED REACTANT/PRODUCT RATES
C     convert from mol/cc-s to kg/m3-s with RCKWRK(NcWT ...)
C
      WDOT(1) = (-RK1           )*RCKWRK(NcWT  )
      WDOT(2) = (-RK1 -      RK2)*RCKWRK(NcWT+1)
      WDOT(3) = ( RK1 -      RK2)*RCKWRK(NcWT+2)
      WDOT(4) = ( RK1 + 2.d0*RK2)*RCKWRK(NcWT+3)
C
c      write(13,*) RK1, RK2, 'mol/cc/s x 1000'
c      write(13,*) RCKWRK(NcWT  ), RCKWRK(NcWT+1), 'g/mol'
c      write(13,*) RCKWRK(NcWT+2), RCKWRK(NcWT+3), 'g/mol'
c      write(13,*) WDOT(1), WDOT(2), WDOT(3)
c      write(13,*) WDOT(4), 'kg/m3/s'
c      write(13,*)
C
      RETURN
      END
C
C ------------------------------------------------------------------- C
C
      SUBROUTINE HTH_EUL1  (T, CPM, ICKWRK, RCKWRK, ENTH)
C
C  START PROLOGUE
C 
C  This subroutine returns enthalpy of each species given temperature
C
C  Enth = enth_ref + cpavg*(T - T_ref)
C   
C  INPUT
C     T      - Temperature.
C                   cgs units - K
C                   Data type - real scalar
C     CPM    - heat capacity
C                   cgs units - J/kg-K
C                   data type - real scalar
C     ICKWRK - Array of integer workspace.
C                   Data type - integer array
C                   Dimension ICKWRK(*) at least LENIWK.
C     RCKWRK - Array of real work space.
C                   Data type - real array
C                   Dimension RCKWRK(*) at least LENRWK.
C
C  OUTPUT
C     ENTH   - Enthalpy
C                    units     - J/kg
c                    data type - real array           
C
C
C  END PROLOGUE
C
C*****double precision
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END double precision
C*****single precision
c        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END single precision
C
      DIMENSION ICKWRK(*), RCKWRK(*), ENTH(*)
      COMMON /CKSTRT/ NMM , NKK , NII , MXSP, MXTB, MXTP, NCP , NCP1,
     1                NCP2, NCP2T,NPAR, NLAR, NFAR, NLAN, NFAL, NREV,
     2                NTHB, NRLT, NWL,  IcMM, IcKK, IcNC, IcPH, IcCH,
     3                IcNT, IcNU, IcNK, IcNS, IcNR, IcLT, IcRL, IcRV,
     4                IcWL, IcFL, IcFO, IcKF, IcTB, IcKN, IcKT, NcAW,
     5                NcWT, NcTT, NcAA, NcCO, NcRV, NcLT, NcRL, NcFL,
     6                NcKT, NcWL, NcRU, NcRC, NcPA, NcK1, NcK2, NcK3,
     7                NcK4, NcI1, NcI2, NcI3, NcI4
C
      DELTAT = T - 298.15d0
      SUMH = CPM * DELTAT !J/kg
C Fuel
c using -74.85 kJ/mol and mass 16g/mol --> 4.678125d6 J/kg
      ENTH(1) = - 4.678125d6   + SUMH
C Oxid
c none
      ENTH(2) =                + SUMH 
C Intr
c using -393.5 kJ/mol and mass 44g/mol
      ENTH(3) = - 8.943182d6   + SUMH
C Prod
c using -241.82 kJ/mol and mass 18g/mol
      ENTH(4) = - 1.343444d7   + SUMH
C Nitro
c none
      ENTH(5) =                + SUMH
C
      RETURN
      END
C
C ------------------------------------------------------------------- C
C
      SUBROUTINE HTH_EUL2dns  (T, CPM, ICKWRK, RCKWRK, ENTH)
C
C  START PROLOGUE
C 
C  This subroutine returns enthalpy of each species given temperature
C
C  Enth = enth_ref + cpavg*(T - T_ref)
C   
C  INPUT
C     T      - Temperature.
C                   cgs units - K
C                   Data type - real scalar
C     CPM    - heat capacity
C                   cgs units - J/kg-K
C                   data type - real scalar
C     ICKWRK - Array of integer workspace.
C                   Data type - integer array
C                   Dimension ICKWRK(*) at least LENIWK.
C     RCKWRK - Array of real work space.
C                   Data type - real array
C                   Dimension RCKWRK(*) at least LENRWK.
C
C  OUTPUT
C     ENTH   - Enthalpy
C                    units     - J/kg
c                    data type - real array           
C
C
C  END PROLOGUE
C
C*****double precision
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END double precision
C*****single precision
c        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END single precision
C
      DIMENSION ICKWRK(*), RCKWRK(*), ENTH(*)
      COMMON /CKSTRT/ NMM , NKK , NII , MXSP, MXTB, MXTP, NCP , NCP1,
     1                NCP2, NCP2T,NPAR, NLAR, NFAR, NLAN, NFAL, NREV,
     2                NTHB, NRLT, NWL,  IcMM, IcKK, IcNC, IcPH, IcCH,
     3                IcNT, IcNU, IcNK, IcNS, IcNR, IcLT, IcRL, IcRV,
     4                IcWL, IcFL, IcFO, IcKF, IcTB, IcKN, IcKT, NcAW,
     5                NcWT, NcTT, NcAA, NcCO, NcRV, NcLT, NcRL, NcFL,
     6                NcKT, NcWL, NcRU, NcRC, NcPA, NcK1, NcK2, NcK3,
     7                NcK4, NcI1, NcI2, NcI3, NcI4
C
      DELTAT = T - 298.15d0
      SUMH = CPM * DELTAT !J/kg
C Fuel
c using -74.85 kJ/mol and mass 16g/mol --> 4.678125d6 J/kg
      ENTH(1) = - 4.678125d6   + SUMH
C Oxid
c none
      ENTH(2) =                + SUMH 
C Intr
c using -103.4 kJ/mol and mass 21g/mol
      ENTH(3) = - 4.92380952d6 + SUMH
C Prod
c using -292.37 kJ/mol and mass 26g/mol
      ENTH(4) = - 1.1245d7     + SUMH
C Nitro
c none
      ENTH(5) =                + SUMH
C
      RETURN
      END
C
C ------------------------------------------------------------------- C
C
      SUBROUTINE HTH_EUL2a  (T, CPM, ICKWRK, RCKWRK, ENTH)
C
C  START PROLOGUE
C 
C  This subroutine returns enthalpy of each species given temperature
C
C  Enth = enth_ref + cpavg*(T - T_ref)
C   
C  INPUT
C     T      - Temperature.
C                   cgs units - K
C                   Data type - real scalar
C     CPM    - heat capacity
C                   cgs units - J/kg-K
C                   data type - real scalar
C     ICKWRK - Array of integer workspace.
C                   Data type - integer array
C                   Dimension ICKWRK(*) at least LENIWK.
C     RCKWRK - Array of real work space.
C                   Data type - real array
C                   Dimension RCKWRK(*) at least LENRWK.
C
C  OUTPUT
C     ENTH   - Enthalpy
C                    units     - J/kg
c                    data type - real array           
C
C
C  END PROLOGUE
C
C*****double precision
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END double precision
C*****single precision
c        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END single precision
C
      DIMENSION ICKWRK(*), RCKWRK(*), ENTH(*)
      COMMON /CKSTRT/ NMM , NKK , NII , MXSP, MXTB, MXTP, NCP , NCP1,
     1                NCP2, NCP2T,NPAR, NLAR, NFAR, NLAN, NFAL, NREV,
     2                NTHB, NRLT, NWL,  IcMM, IcKK, IcNC, IcPH, IcCH,
     3                IcNT, IcNU, IcNK, IcNS, IcNR, IcLT, IcRL, IcRV,
     4                IcWL, IcFL, IcFO, IcKF, IcTB, IcKN, IcKT, NcAW,
     5                NcWT, NcTT, NcAA, NcCO, NcRV, NcLT, NcRL, NcFL,
     6                NcKT, NcWL, NcRU, NcRC, NcPA, NcK1, NcK2, NcK3,
     7                NcK4, NcI1, NcI2, NcI3, NcI4
C
      DELTAT = T - 298.15d0
      SUMH = CPM * DELTAT !J/kg
C Fuel
c using -74.0 kJ/mol and mass 16g/mol --> 4.625d6 J/kg
      ENTH(1) = - 4.625d6    + SUMH
C Oxid
c none
      ENTH(2) =              + SUMH 
C Intr
c using -73.0 kJ/mol and mass 64/3 g/mol
      ENTH(3) = - 3.421875d6 + SUMH
C Prod
c using -292.0 kJ/mol and mass 80/3 g/mol
      ENTH(4) = - 1.095d7    + SUMH
C Nitro
c none
      ENTH(5) =              + SUMH
C
      RETURN
      END
C
C ------------------------------------------------------------------- C
C
      SUBROUTINE HTH_CM2  (T, CPM, ICKWRK, RCKWRK, ENTH)
C
C  START PROLOGUE
C 
C  This subroutine returns enthalpy of each species given temperature
C
C  Enth = enth_ref + cpavg*(T - T_ref)
C   
C  INPUT
C     T      - Temperature.
C                   cgs units - K
C                   Data type - real scalar
C     CPM    - heat capacity
C                   cgs units - J/kg-K
C                   data type - real scalar
C     ICKWRK - Array of integer workspace.
C                   Data type - integer array
C                   Dimension ICKWRK(*) at least LENIWK.
C     RCKWRK - Array of real work space.
C                   Data type - real array
C                   Dimension RCKWRK(*) at least LENRWK.
C
C  OUTPUT
C     ENTH   - Enthalpy
C                    units     - J/kg
c                    data type - real array           
C
C
C  END PROLOGUE
C
C*****double precision
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END double precision
C*****single precision
c        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END single precision
C
      DIMENSION ICKWRK(*), RCKWRK(*), ENTH(*)
      COMMON /CKSTRT/ NMM , NKK , NII , MXSP, MXTB, MXTP, NCP , NCP1,
     1                NCP2, NCP2T,NPAR, NLAR, NFAR, NLAN, NFAL, NREV,
     2                NTHB, NRLT, NWL,  IcMM, IcKK, IcNC, IcPH, IcCH,
     3                IcNT, IcNU, IcNK, IcNS, IcNR, IcLT, IcRL, IcRV,
     4                IcWL, IcFL, IcFO, IcKF, IcTB, IcKN, IcKT, NcAW,
     5                NcWT, NcTT, NcAA, NcCO, NcRV, NcLT, NcRL, NcFL,
     6                NcKT, NcWL, NcRU, NcRC, NcPA, NcK1, NcK2, NcK3,
     7                NcK4, NcI1, NcI2, NcI3, NcI4
C
      DELTAT = T - 298.15d0
      SUMH = CPM * DELTAT !J/kg
C CH4
c using -74.87 kJ/mol and mass 16.04303g/mol --> 4.666824d6 J/kg
      ENTH(1) = - 4.666824d6 + SUMH
C CO
c using -110.53 and 28.01055
      ENTH(2) = - 3.946013d6 + SUMH 
C CO2
c using -393.52 and 44.00995
      ENTH(3) = - 8.941614d6 + SUMH
C H2
c none
      ENTH(4) =              + SUMH
C O2
c none
      ENTH(5) =              + SUMH
C N2
c none
      ENTH(6) =              + SUMH
C
      RETURN
      END
C
C ------------------------------------------------------------------- C
C
      SUBROUTINE HTH_WD  (T, CPM, ICKWRK, RCKWRK, ENTH)
C
C  START PROLOGUE
C 
C  This subroutine returns enthalpy of each species given temperature
C
C  Enth = enth_ref + cpavg*(T - T_ref)
C   
C  INPUT
C     T      - Temperature.
C                   cgs units - K
C                   Data type - real scalar
C     CPM    - heat capacity
C                   cgs units - J/kg-K
C                   data type - real scalar
C     ICKWRK - Array of integer workspace.
C                   Data type - integer array
C                   Dimension ICKWRK(*) at least LENIWK.
C     RCKWRK - Array of real work space.
C                   Data type - real array
C                   Dimension RCKWRK(*) at least LENRWK.
C
C  OUTPUT
C     ENTH   - Enthalpy
C                    units     - J/kg
c                    data type - real array           
C
C
C  END PROLOGUE
C
C*****double precision
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END double precision
C*****single precision
c        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END single precision
C
      DIMENSION ICKWRK(*), RCKWRK(*), ENTH(*)
      COMMON /CKSTRT/ NMM , NKK , NII , MXSP, MXTB, MXTP, NCP , NCP1,
     1                NCP2, NCP2T,NPAR, NLAR, NFAR, NLAN, NFAL, NREV,
     2                NTHB, NRLT, NWL,  IcMM, IcKK, IcNC, IcPH, IcCH,
     3                IcNT, IcNU, IcNK, IcNS, IcNR, IcLT, IcRL, IcRV,
     4                IcWL, IcFL, IcFO, IcKF, IcTB, IcKN, IcKT, NcAW,
     5                NcWT, NcTT, NcAA, NcCO, NcRV, NcLT, NcRL, NcFL,
     6                NcKT, NcWL, NcRU, NcRC, NcPA, NcK1, NcK2, NcK3,
     7                NcK4, NcI1, NcI2, NcI3, NcI4
C
      DELTAT = T - 298.15d0
      SUMH = CPM * DELTAT !J/kg
C CH4
c using -74.87 kJ/mol and mass 16.04303g/mol --> 4.666824d6 J/kg
      ENTH(1) = - 4.666824d6 + SUMH
C CO
c using -110.53 and 28.01055
      ENTH(2) = - 3.946013d6 + SUMH 
C CO2
c using -393.52 and 44.00995
      ENTH(3) = - 8.941614d6 + SUMH
C H2O
c using -241.82 and 18.01534 
      ENTH(4) = - 1.342301d7 + SUMH
C O2
c none
      ENTH(5) =              + SUMH
C N2
c none
      ENTH(6) =              + SUMH
C
      RETURN
      END
C
C ------------------------------------------------------------------- 
C
      SUBROUTINE HTH_L1  (T, CPM, ICKWRK, RCKWRK, ENTH)
C
C  START PROLOGUE
C 
C  This subroutine returns enthalpy of each species given temperature
C
C  Enth = enth_ref + cpavg*(T - T_ref)
C   
C  INPUT
C     T      - Temperature.
C                   cgs units - K
C                   Data type - real scalar
C     CPM    - heat capacity
C                   cgs units - J/kg-K
C                   data type - real scalar
C     ICKWRK - Array of integer workspace.
C                   Data type - integer array
C                   Dimension ICKWRK(*) at least LENIWK.
C     RCKWRK - Array of real work space.
C                   Data type - real array
C                   Dimension RCKWRK(*) at least LENRWK.
C
C  OUTPUT
C     ENTH   - Enthalpy
C                    units     - J/kg
c                    data type - real array           
C
C
C  END PROLOGUE
C
C*****double precision
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END double precision
C*****single precision
c        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END single precision
C
      DIMENSION ICKWRK(*), RCKWRK(*), ENTH(*)
      COMMON /CKSTRT/ NMM , NKK , NII , MXSP, MXTB, MXTP, NCP , NCP1,
     1                NCP2, NCP2T,NPAR, NLAR, NFAR, NLAN, NFAL, NREV,
     2                NTHB, NRLT, NWL,  IcMM, IcKK, IcNC, IcPH, IcCH,
     3                IcNT, IcNU, IcNK, IcNS, IcNR, IcLT, IcRL, IcRV,
     4                IcWL, IcFL, IcFO, IcKF, IcTB, IcKN, IcKT, NcAW,
     5                NcWT, NcTT, NcAA, NcCO, NcRV, NcLT, NcRL, NcFL,
     6                NcKT, NcWL, NcRU, NcRC, NcPA, NcK1, NcK2, NcK3,
     7                NcK4, NcI1, NcI2, NcI3, NcI4
C
      DELTAT = T - 298.15d0
      SUMH = CPM * DELTAT !J/kg
C CH4
c using -74.87 kJ/mol and mass 16.04303g/mol --> 4.666824d6 J/kg
      ENTH(1) = - 4.666824d6 + SUMH
C CO2
c using -393.52 and 44.00995
      ENTH(2) = - 8.941614d6 + SUMH
C H2O
c using -241.82 and 18.01534 
      ENTH(3) = - 1.342301d7 + SUMH
C O2
c none
      ENTH(4) =              + SUMH
C N2
c none
      ENTH(5) =              + SUMH
C
      RETURN
      END
C
C ------------------------------------------------------------------- 
C
      SUBROUTINE JYCH(T, CPM, ICKWRK, RCKWRK, ENTH)
C
C  START PROLOGUE
C 
C  This subroutine returns enthalpy of each species given temperature
C
C  Enth = enth_ref + cpavg*(T - T_ref)
C   
C  INPUT
C     T      - Temperature.
C                   cgs units - K
C                   Data type - real scalar
C     CPM    - heat capacity
C                   cgs units - J/kg-K
C                   data type - real scalar
C     ICKWRK - Array of integer workspace.
C                   Data type - integer array
C                   Dimension ICKWRK(*) at least LENIWK.
C     RCKWRK - Array of real work space.
C                   Data type - real array
C                   Dimension RCKWRK(*) at least LENRWK.
C
C  OUTPUT
C     ENTH   - Enthalpy
C                    units     - J/kg
c                    data type - real array           
C
C
C  END PROLOGUE
C
C*****double precision
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END double precision
C*****single precision
c        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END single precision
C
      DIMENSION ICKWRK(*), RCKWRK(*), ENTH(*)
      COMMON /CKSTRT/ NMM , NKK , NII , MXSP, MXTB, MXTP, NCP , NCP1,
     1                NCP2, NCP2T,NPAR, NLAR, NFAR, NLAN, NFAL, NREV,
     2                NTHB, NRLT, NWL,  IcMM, IcKK, IcNC, IcPH, IcCH,
     3                IcNT, IcNU, IcNK, IcNS, IcNR, IcLT, IcRL, IcRV,
     4                IcWL, IcFL, IcFO, IcKF, IcTB, IcKN, IcKT, NcAW,
     5                NcWT, NcTT, NcAA, NcCO, NcRV, NcLT, NcRL, NcFL,
     6                NcKT, NcWL, NcRU, NcRC, NcPA, NcK1, NcK2, NcK3,
     7                NcK4, NcI1, NcI2, NcI3, NcI4
C
      DELTAT = T - 298.15d0
      SUMH = CPM * DELTAT !J/kg
C *** 
C DATA TAKEN FROM 
C http://www.ohio.edu/mechanical/thermo/property_tables/combustion/Enth_Formation.html
C ***
C H2
C none 2.01594g/mol
      ENTH(1) =              + SUMH
C H
C none 1.00797g/mol
      ENTH(2) =              + SUMH
C O
C using +249.19 and 15.99940g/mol
      ENTH(3) = + 1.557496d7 + SUMH
C O2
c none 31.99880g/mol
      ENTH(4) =              + SUMH
C OH
C using +39.46 and 17.00737g/mol
      ENTH(5) = + 2.320171d6 + SUMH
C H2O
c using -241.82 and 18.01534 
      ENTH(6) = - 1.342301d7 + SUMH
C CH4
c using -74.85 kJ/mol and mass 16.04303g/mol --> 4.666824d6 J/kg
      ENTH(7) = - 4.665578d6 + SUMH
C CO
C using -110.53 and 28.01055g/mol
      ENTH(8) = - 3.946013d6 + SUMH
C CO2
c using -393.52 and 44.00995
      ENTH(9) = - 8.941614d6 + SUMH
C N2
c none
      ENTH(10) =             + SUMH
C
      RETURN
      END
C
C----------------------------------------------------------------------C
C
C
      SUBROUTINE JYCWYP (P, T, Y, ICKWRK, RCKWRK, WDOT)
C This routine can be changed to be compatible with CKWYR by
C replacing the above subroutine heading with the line below
C and by changing the section of computing XCON.
C      SUBROUTINE CKWYR (RHO, T, Y, ICKWRK, RCKWRK, WDOT)
C
C     This routine is automatically produced by CARM
C     for computing sources of reduced mechanism.
C     (Version: 1.0.13      Last updated on Dec-24-98  )
C
C last updated 3-15-99 by J-Y Chen                                                
C
C      SUMMARY OF REDUCED MECHANISM:
C      TOTAL NUMBER OF SPECIES= 10 WITH   6 STEPS
C ( 1)    2O = O2                                                               
C ( 2)   H + O = OH                                                             
C ( 3)   H2 + O = H + OH                                                        
C ( 4)   O +   .50 CH4 =   .50 H2 +   .50 H +   .50 OH +   .50 CO               
C ( 5)   O + CO = CO2                                                           
C ( 6)   O + H2O +   .25 CO =   .25 H2 +   .25 H + O2 +   .25 OH +              
C          .25 CH4                                                              
C
C  INPUT
C     P      - Pressure.
C                   cgs units - dynes/cm**2
C     (RHO   - Density.)
C                  (cgs units - gm/cm**3)
C                   Data type - real scalar
C     T      - Temperature.
C                   cgs units - K
C                   Data type - real scalar
C     Y      - Mass fractions of species.
C                   Data type - real array
C                   Dimension Y(*) at least  10 (total number of species)
C                   Species must be arranged in the order of:
C                   Y( 1)=    H2              
C                   Y( 2)=    H               
C                   Y( 3)=    O               
C                   Y( 4)=    O2              
C                   Y( 5)=    OH              
C                   Y( 6)=    H2O             
C                   Y( 7)=    CH4             
C                   Y( 8)=    CO              
C                   Y( 9)=    CO2             
C                   Y(10)=    N2              
C     ICKWRK - Dummy Array of integer workspace.
C                   Data type - integer array
C                   (Not used; simply to be Compatible with Chemkin)
C     RCKWRK - Dummy Array of real work space.
C                   Data type - real array
C                   (Not used; simply to be Compatible with Chemkin)
C
C  OUTPUT
C     WDOT   - Chemical molar production rates of the species.
C                   cgs units - moles/(cm**3*sec)
C                   Data type - real array
C                   Dimension WDOT(*) at least 10
C
C  END PROLOGUE
C
C*****double precision
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END double precision
C*****single precision
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END single precision
C
      PARAMETER (IGLO=6,IREAC=177,KK=10,KSS=21)                                 
      PARAMETER (NITER=50, RELACC=1.D-5)
      DIMENSION ICKWRK(*), RCKWRK(*), Y(KK), WDOT(KK)
      DIMENSION XM(IREAC), RF(IREAC), RB(IREAC), W(IREAC)
      DIMENSION RKF(IGLO), XCON(KK), WT(KK)
      DIMENSION B(KSS), ABV(KSS), DEN(KSS)
      LOGICAL LITER, LBOUND(KSS)
      SAVE RF,RB,TOLD
      DATA SMALL/1.D-50/,TOLD/0.D0/
      DATA RU/8.314510D7/
      DATA WT/  2.016,  1.008, 15.999, 31.999, 17.007, 18.015, 16.043,          
     &         28.011, 44.010, 28.013/                                          
C
C   compute concentrations from mass fractions 
C section for routine compatible with CKWYP
      SUMYOW = 0.0
      DO K = 1, KK
        SUMYOW = SUMYOW + Y(K)/WT(K)
      ENDDO
      SUMYOW = SUMYOW*T*RU
      DO K = 1, KK
        XCON(K) = P*Y(K) / ( SUMYOW*WT(K) )
      ENDDO
C Replacing the above section by following lines for CKWYR
C      DO K = 1, KK
C        XCON(K) = RHO*Y(K) / WT(K) 
C      ENDDO
C
C   SET LOCAL VALUES FOR THE CONCENTRATIONS
C
      XH2 = MAX( XCON(1), SMALL )                                               
      XH = MAX( XCON(2), SMALL )                                                
      XO = MAX( XCON(3), SMALL )                                                
      XO2 = MAX( XCON(4), SMALL )                                               
      XOH = MAX( XCON(5), SMALL )                                               
      XH2O = MAX( XCON(6), SMALL )                                              
      XCH4 = MAX( XCON(7), SMALL )                                              
      XCO = MAX( XCON(8), SMALL )                                               
      XCO2 = MAX( XCON(9), SMALL )                                              
      XN2 = MAX( XCON(10), SMALL )                                              

      BIG = 0.0
      DO 20 K = 1, KK
        IF( XCON(K) .GT. 0.0 ) THEN
          BIG = MAX( BIG, XCON(K) )
        ENDIF
        WDOT(K) = 0.0
 20   CONTINUE
C
C   EFFECTIVE THIRD BODY FOR ALL REACTIONS
C
      CALL THIRDB( XM, XH2, XH, XO, XO2, XOH, XH2O, XCH4, XCO, XCO2,            
     &             XN2 )                                                        
C
C   SET THE ELEMENTARY RATES
C
      CALL ELEMRATE( RF, RB, T, XM, TOLD )
C
C   EXPRESSIONS FOR STEADY-STATE SPECIES
C
      DO 40 K = 1, KSS
        B(K) = 0.0
        LBOUND(K) = .TRUE.
 40   CONTINUE
      CALL UPVALUE( 1, B, BIG, LBOUND, CONMAX, XHO2, XH2O2, XCH3, XCH2O,        
     &              XHCO, XCH2CO, XCH3OH, XCH3O, XCH2OH, XC2H6, XCH2S,          
     &              XCH2, XCH, XC2H4, XC2H3, XC2H2, XHCCO, XC2H, XHCCOH,        
     &              XC, XC2H5 )                                                 
      ADJ = 1.D0/BIG
      DO 30 N = 1, NITER

          LITER = .TRUE.
          CONMAX = 0.0
        CALL STEADYE( ABV, DEN, RF, RB, XM, ADJ, CONMAX, SMALL, LBOUND,
     &                LITER, XH2, XH, XO, XO2, XOH, XH2O, XHO2, XH2O2,          
     &                XC, XCH, XCH2, XCH2S, XCH3, XCH4, XCO, XCO2, XHCO,        
     &                XCH2O, XCH2OH, XCH3O, XCH3OH, XC2H, XC2H2, XC2H3,         
     &                XC2H4, XC2H5, XC2H6, XHCCO, XCH2CO, XHCCOH, XN2 )         

        IF( LITER .AND. CONMAX .LT. RELACC ) GO TO 35


 30   CONTINUE

 35   CONTINUE
C
C   NET PRODUCTION RATES FOR SKELETAL MECHANISM
C
      CALL NETRATE( W, RF, RB, XM, XH2, XH, XO, XO2, XOH, XH2O, XHO2,           
     &  XH2O2, XC, XCH, XCH2, XCH2S, XCH3, XCH4, XCO, XCO2, XHCO, XCH2O,        
     &  XCH2OH, XCH3O, XCH3OH, XC2H, XC2H2, XC2H3, XC2H4, XC2H5, XC2H6,         
     &  XHCCO, XCH2CO, XHCCOH, XN2 )                                            
C
C   GLOBAL RATES FOR REDUCED MECHANISM
C
      RKF(1) = +W(1) +W(7) +W(9) +W(10) -W(13) -W(14) +W(25) +W(26)             
     &         -W(31) -W(38) +W(43) -W(46) -W(55) -W(61) -W(66) +W(83)          
     &         +W(84) +W(85) +W(86) +W(87) +W(88) +W(89) +W(91) +W(92)          
     &         +W(93) +W(94) +W(95) +W(96) +W(97) +W(98) +W(101)                
     &         +W(102) +W(103) +W(104) +W(105) +W(109) +W(111) +W(112)          
     &         +W(113) +W(114) -W(120) -W(122) +W(132) -W(133) -W(144)          
     &         +W(153) -W(160) -W(166) -W(167) -W(168) +W(173) -W(176)          
      RKF(2) = +W(2) -W(7) +W(8) -W(9) -W(10) + 2.00*W(13) +W(14)               
     &         +W(22) +W(24) -W(25) -W(26) -W(30) +W(33) +W(34) +W(35)          
     &         +W(36) +W(38) +W(39) +W(40) +W(41) +W(42) +W(46) +W(49)          
     &         +W(50) +W(51) +W(52) +W(54) + 2.00*W(55) +W(56) +W(57)           
     &         +W(59) +W(61) +W(63) +W(66) +W(72) +W(73) +W(74) +W(76)          
     &         -W(83) -W(84) -W(86) -W(87) -W(88) -W(89) -W(90)                 
     &         - 2.00*W(91) -W(92) -W(94) -W(96) -W(97) -W(98) -W(99)           
     &         +W(100) -W(101) -W(102) -W(103) -W(104) -W(105) -W(106)          
     &         -W(112) -W(113) -W(114) -W(124) -W(125) -W(126) -W(127)          
     &         -W(129) -W(130) -W(132) +W(133) +W(137) +W(140) +W(144)          
     &         +W(145) +W(147) +W(158) + 2.00*W(160) +W(166) +W(167)            
     &         + 2.00*W(168) -W(171) -W(172) -W(173) +W(176)                    
      RKF(3) = +W(3) +  .50*W(6) -  .50*W(9) -  .50*W(10) + 1.50*W(13)          
     &         + 1.50*W(14) +  .50*W(20) +  .50*W(21) -  .50*W(23)              
     &         +W(24) -  .50*W(25) -  .50*W(26) +  .50*W(28)                    
     &         +  .50*W(29) -  .50*W(30) -W(39) -W(40) -W(41) -W(42)            
     &         -  .50*W(43) -  .50*W(44) -W(45) -W(47) -  .50*W(48)             
     &         +  .50*W(50) -  .50*W(51) -W(53) +  .50*W(55) -W(58)             
     &         -W(60) +  .50*W(61) -W(65) +  .50*W(66) -W(68) -W(69)            
     &         -W(71) -W(75) -W(77) -W(78) -  .50*W(79) -  .50*W(80)            
     &         -  .50*W(83) +  .50*W(84) -  .50*W(86) -  .50*W(87)              
     &         -  .50*W(88) -  .50*W(89) -  .50*W(90) -W(91)                    
     &         +  .50*W(93) -  .50*W(94) -  .50*W(95) -W(96)                    
     &         -  .50*W(97) -  .50*W(98) +W(100) -  .50*W(101)                  
     &         -  .50*W(102) -  .50*W(103) -  .50*W(104) -  .50*W(105)          
     &         +  .50*W(106) -  .50*W(109) +  .50*W(111) -  .50*W(112)          
     &         -  .50*W(113) -  .50*W(119) -  .50*W(122) -  .50*W(124)          
     &         -W(125) -  .50*W(127) +W(128) -  .50*W(129)                      
     &         -  .50*W(130) -W(132) +W(133) + 1.50*W(136) +W(137)              
     &         +  .50*W(138) +  .50*W(139) +  .50*W(140) -  .50*W(142)          
     &         +W(144) +  .50*W(145) +W(146) -  .50*W(148)                      
     &         -  .50*W(151) -  .50*W(152) -  .50*W(153) -  .50*W(155)          
     &         -  .50*W(156) + 1.50*W(160) + 1.50*W(166) + 1.50*W(167)          
     &         + 1.50*W(168) -  .50*W(171) +W(172) -W(173)                      
     &         +  .50*W(176)                                                    
      RKF(4) = +W(4) +W(5) +  .50*W(6) +  .50*W(9) +  .50*W(10) +W(11)          
     &         -  .50*W(13) -  .50*W(14) +W(15) +W(18) +W(19)                   
     &         +  .50*W(20) +  .50*W(21) + 1.50*W(23) -W(24)                    
     &         +  .50*W(25) +  .50*W(26) +W(27) +  .50*W(28)                    
     &         +  .50*W(29) + 1.50*W(30) -W(33) -W(34) -W(35) -W(36)            
     &         -  .50*W(43) +  .50*W(44) +W(45) +W(46) +W(47)                   
     &         + 1.50*W(48) -W(49) - 1.50*W(50) -  .50*W(51) -W(52)             
     &         +W(53) -W(54) -  .50*W(55) +W(58) -W(59) -  .50*W(61)            
     &         -W(63) -  .50*W(66) +W(68) +W(69) +W(71) -W(72) -W(73)           
     &         -W(74) +W(75) -W(76) +W(77) +W(78) +  .50*W(79)                  
     &         +  .50*W(80) -  .50*W(83) -  .50*W(84) - 2.00*W(85)              
     &         -  .50*W(86) +  .50*W(87) +  .50*W(88) +  .50*W(89)              
     &         + 1.50*W(90) +W(91) -W(92) - 1.50*W(93) -  .50*W(94)             
     &         -  .50*W(95) +W(96) +  .50*W(97) +  .50*W(98) -W(100)            
     &         +  .50*W(101) -  .50*W(102) -  .50*W(103) +  .50*W(104)          
     &         +  .50*W(105) +  .50*W(106) -  .50*W(109) - 1.50*W(111)          
     &         +  .50*W(112) +  .50*W(113) + 1.50*W(119) +W(120)                
     &         + 1.50*W(122) + 1.50*W(124) +W(125) +W(126)                      
     &         +  .50*W(127) -W(128) + 1.50*W(129) + 1.50*W(130)                
     &         +W(132) - 1.50*W(136) - 2.00*W(137) -  .50*W(138)                
     &         -  .50*W(139) - 1.50*W(140) +  .50*W(142) -  .50*W(145)          
     &         -W(146) -W(147) +  .50*W(148) +  .50*W(151)                      
     &         +  .50*W(152) -  .50*W(153) +  .50*W(155) +  .50*W(156)          
     &         -W(158) - 1.50*W(160) -  .50*W(166) -  .50*W(167)                
     &         - 1.50*W(168) -W(169) -W(170) + 1.50*W(171)                      
     &         +  .50*W(176)                                                    
      RKF(5) = +W(12) +W(14) +W(30) +W(31) +W(99) +W(120) -W(132)               
     &         -W(153)                                                          
      RKF(6) = +W(16) +W(17) -W(43) -W(44) -W(48) -W(56) -W(57) +W(60)          
     &         +W(61) +W(65) +W(66) -W(84) -W(86) -W(87) -W(88) -W(89)          
     &         -W(93) -W(95) -W(96) -W(97) -W(98) -W(100) -W(101)               
     &         -W(104) -W(105) -W(109) -W(111) -W(112) -W(113) -W(114)          
     &         -W(119) +W(127) -W(145) -W(155) +W(169) +W(170)                  
C
C   SPECIES PRODUCTION RATES
C
C H2                                                                            
      WDOT(1) = - RKF(3) +  .50*RKF(4) +  .25*RKF(6)                            
C H                                                                             
      WDOT(2) = - RKF(2) + RKF(3) +  .50*RKF(4) +  .25*RKF(6)                   
C O                                                                             
      WDOT(3) = - 2.00*RKF(1) - RKF(2) - RKF(3) - RKF(4) - RKF(5)               
     &          - RKF(6)                                                        
C O2                                                                            
      WDOT(4) = + RKF(1) + RKF(6)                                               
C OH                                                                            
      WDOT(5) = + RKF(2) + RKF(3) +  .50*RKF(4) +  .25*RKF(6)                   
C H2O                                                                           
      WDOT(6) = - RKF(6)                                                        
C CH4                                                                           
      WDOT(7) = -  .50*RKF(4) +  .25*RKF(6)                                     
C CO                                                                            
      WDOT(8) = +  .50*RKF(4) - RKF(5) -  .25*RKF(6)                            
C CO2                                                                           
      WDOT(9) = + RKF(5)                                                        
C N2                                                                            
      WDOT(10) = 0.0                                                            

      DO K = 1, KK
        IF( XCON(K) .LE. 0.D0 .AND. WDOT(K) .LT. 0.D0 ) WDOT(K) = 0.D0
      END DO

      RETURN
      END
C
C -------------------------------------------------------------------- C
C
      SUBROUTINE THIRDB( XM, XH2, XH, XO, XO2, XOH, XH2O, XCH4, XCO,            
     &                   XCO2, XN2 )                                            

C*****double precision
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END double precision
C*****single precision
C       IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END single precision

      DIMENSION XM(*)
C
C   EFFECTIVE THIRD BODY FOR ALL REACTIONS
C
      XM(1) = + 2.400*XH2+XH+XO+XO2+XOH+15.400*XH2O+ 2.000*XCH4                 
     &        + 1.750*XCO+ 3.600*XCO2+XN2                                       
      XM(2) = + 2.000*XH2+XH+XO+XO2+XOH+ 6.000*XH2O+ 2.000*XCH4                 
     &        + 1.500*XCO+ 2.000*XCO2+XN2                                       
      XM(12) = + 2.000*XH2+XH+XO+ 6.000*XO2+XOH+ 6.000*XH2O+ 2.000*XCH4         
     &         + 1.500*XCO+ 3.500*XCO2+XN2                                      
      XM(33) = +XH2+XH+XO+XOH+XCH4+  .750*XCO+ 1.500*XCO2                       
      XM(39) = +XH+XO+XO2+XOH+ 2.000*XCH4+XCO+XN2                               
      XM(43) = +  .730*XH2+XH+XO+XO2+XOH+ 3.650*XH2O+ 2.000*XCH4+XCO            
     &         +XCO2+XN2                                                        
      XM(50) = + 2.000*XH2+XH+XO+XO2+XOH+ 6.000*XH2O+ 2.000*XCH4                
     &         + 1.500*XCO+ 2.000*XCO2+XN2                                      
      XM(52) = + 2.000*XH2+XH+XO+XO2+XOH+ 6.000*XH2O+ 2.000*XCH4                
     &         + 1.500*XCO+ 2.000*XCO2+XN2                                      
      XM(54) = + 2.000*XH2+XH+XO+XO2+XOH+ 6.000*XH2O+ 2.000*XCH4                
     &         + 1.500*XCO+ 2.000*XCO2+XN2                                      
      XM(56) = + 2.000*XH2+XH+XO+XO2+XOH+ 6.000*XH2O+ 2.000*XCH4                
     &         + 1.500*XCO+ 2.000*XCO2+XN2                                      
      XM(57) = + 2.000*XH2+XH+XO+XO2+XOH+ 6.000*XH2O+ 2.000*XCH4                
     &         + 1.500*XCO+ 2.000*XCO2+XN2                                      
      XM(59) = + 2.000*XH2+XH+XO+XO2+XOH+ 6.000*XH2O+ 2.000*XCH4                
     &         + 1.500*XCO+ 2.000*XCO2+XN2                                      
      XM(63) = + 2.000*XH2+XH+XO+XO2+XOH+ 6.000*XH2O+ 2.000*XCH4                
     &         + 1.500*XCO+ 2.000*XCO2+XN2                                      
      XM(70) = + 2.000*XH2+XH+XO+XO2+XOH+ 6.000*XH2O+ 2.000*XCH4                
     &         + 1.500*XCO+ 2.000*XCO2+XN2                                      
      XM(71) = + 2.000*XH2+XH+XO+XO2+XOH+ 6.000*XH2O+ 2.000*XCH4                
     &         + 1.500*XCO+ 2.000*XCO2+XN2                                      
      XM(72) = + 2.000*XH2+XH+XO+XO2+XOH+ 6.000*XH2O+ 2.000*XCH4                
     &         + 1.500*XCO+ 2.000*XCO2+XN2                                      
      XM(74) = + 2.000*XH2+XH+XO+XO2+XOH+ 6.000*XH2O+ 2.000*XCH4                
     &         + 1.500*XCO+ 2.000*XCO2+XN2                                      
      XM(76) = + 2.000*XH2+XH+XO+XO2+XOH+ 6.000*XH2O+ 2.000*XCH4                
     &         + 1.500*XCO+ 2.000*XCO2+XN2                                      
      XM(83) = + 2.000*XH2+XH+XO+XO2+XOH+ 6.000*XH2O+ 2.000*XCH4                
     &         + 1.500*XCO+ 2.000*XCO2+XN2                                      
      XM(85) = + 2.000*XH2+XH+XO+XO2+XOH+ 6.000*XH2O+ 2.000*XCH4                
     &         + 1.500*XCO+ 2.000*XCO2+XN2                                      
      XM(95) = + 2.000*XH2+XH+XO+XO2+XOH+ 6.000*XH2O+ 2.000*XCH4                
     &         + 1.500*XCO+ 2.000*XCO2+XN2                                      
      XM(131) = + 2.000*XH2+XH+XO+XO2+XOH+ 6.000*XH2O+ 2.000*XCH4               
     &          + 1.500*XCO+ 2.000*XCO2+XN2                                     
      XM(140) = + 2.000*XH2+XH+XO+XO2+XOH+ 6.000*XH2O+ 2.000*XCH4               
     &          + 1.500*XCO+ 2.000*XCO2+XN2                                     
      XM(147) = + 2.000*XH2+XH+XO+XO2+XOH+ 6.000*XH2O+ 2.000*XCH4               
     &          + 1.500*XCO+ 2.000*XCO2+XN2                                     
      XM(158) = + 2.000*XH2+XH+XO+XO2+XOH+ 6.000*XH2O+ 2.000*XCH4               
     &          + 1.500*XCO+ 2.000*XCO2+XN2                                     
      XM(167) = + 2.000*XH2+XH+XO+XO2+XOH+ 2.000*XCH4+ 1.500*XCO                
     &          + 2.000*XCO2+XN2                                                
      XM(174) = + 2.000*XH2+XH+XO+XO2+XOH+ 6.000*XH2O+ 2.000*XCH4               
     &          + 1.500*XCO+ 2.000*XCO2+XN2                                     

      RETURN
      END
C
C -------------------------------------------------------------------- C
C
      SUBROUTINE ELEMRATE( RF, RB, T, XM, TOLD )

C*****double precision
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END double precision
C*****single precision
C       IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END single precision

      DIMENSION RF(*), RB(*), XM(*)
      DATA SMALL/1.D-50/

C   RATE CONSTANTS FOR SKELETAL MECHANSIM

      RUC = 8.31451D0/4.184D0
      ALOGT = DLOG(T)
      RTR = 1.0D3/(RUC*T)
      TM1 = 1.0/T                                                               
      TM2 = 1.0/T/T                                                             
      TM3 = 1.0/T/T/T                                                           
      TP1 = T                                                                   
      TP2 = T*T                                                                 
      TP3 = T*T*T                                                               

      IF(ABS(T-TOLD) .GT. 1.D-3) THEN
      RF(1) =  1.200E+17 * EXP( -1.000*ALOGT                )                   
      RB(1) = EXP(+( 1.59693E+10)*TM3+(-1.11204E+08)*TM2                        
     &        +( 2.30331E+05)*TM1+(-3.24217E+02)+( 2.25674E-01)*TP1             
     &        +(-6.88297E-05)*TP2+( 8.02658E-09)*TP3)                           
      RF(2) =  5.000E+17 * EXP( -1.000*ALOGT                )                   
      RB(2) = EXP(+( 2.42356E+06)*TM3+(-3.77760E+04)*TM2                        
     &        +(-5.07137E+04)*TM1+( 3.39494E+01)+(-2.73497E-04)*TP1             
     &        +(-2.08792E-08)*TP2+( 4.46337E-12)*TP3)                           
      RF(3) =  5.000E+04 * EXP(  2.670*ALOGT -  6.29068*RTR )                   
      RB(3) = EXP(+(-2.84005E+07)*TM3+( 3.35502E+05)*TM2                        
     &        +(-4.13521E+03)*TM1+( 2.85383E+01)+( 1.87759E-03)*TP1             
     &        +(-3.32578E-07)*TP2+( 3.02311E-11)*TP3)                           
      RF(4) =  2.000E+13                                                        
      RB(4) = EXP(+(-1.27958E+06)*TM3+(-1.70702E+04)*TM2                        
     &        +(-2.65672E+04)*TM1+( 3.03272E+01)+( 3.29220E-04)*TP1             
     &        +(-5.47023E-08)*TP2+( 4.96502E-12)*TP3)                           
      RF(5) =  9.630E+06 * EXP(  2.000*ALOGT -  4.00043*RTR )                   
      RB(5) = EXP(+(-1.75020E+07)*TM3+( 1.94410E+05)*TM2                        
     &        +(-1.05483E+04)*TM1+( 2.71317E+01)+( 2.13131E-03)*TP1             
     &        +(-3.90188E-07)*TP2+( 3.36872E-11)*TP3)                           
      RF(6) =  5.700E+13                                                        
      RB(6) = EXP(+(-2.60457E+10)*TM3+( 2.48595E+08)*TM2                        
     &        +(-8.38999E+05)*TM1+( 1.04210E+03)+(-6.59885E-01)*TP1             
     &        +( 2.05772E-04)*TP2+(-2.43373E-08)*TP3)                           
      RF(7) =  8.000E+13                                                        
      RB(7) = EXP(+(-7.62928E+06)*TM3+( 9.66395E+04)*TM2                        
     &        +(-4.63265E+04)*TM1+( 3.48539E+01)+(-6.64802E-04)*TP1             
     &        +( 1.64816E-07)*TP2+(-1.65701E-11)*TP3)                           
      RF(8) =  1.500E+13                                                        
      RB(8) = EXP(+(-3.55600E+10)*TM3+( 3.15390E+08)*TM2                        
     &        +(-9.92837E+05)*TM1+( 1.18310E+03)+(-7.27083E-01)*TP1             
     &        +( 2.20400E-04)*TP2+(-2.54940E-08)*TP3)                           
      RF(9) =  1.500E+13                                                        
      RB(9) = EXP(+( 1.28016E+06)*TM3+( 8.52154E+03)*TM2                        
     &        +(-5.05426E+04)*TM1+( 3.21800E+01)+(-5.91173E-04)*TP1             
     &        +( 1.51533E-07)*TP2+(-1.56605E-11)*TP3)                           
      RF(10) =  8.430E+13                                                       
      RB(10) = EXP(+(-1.21447E+07)*TM3+( 1.56172E+05)*TM2                       
     &         +(-3.51673E+04)*TM1+( 3.58716E+01)+(-6.92086E-04)*TP1            
     &         +( 1.66756E-07)*TP2+(-1.61858E-11)*TP3)                          
      RF(11) =  1.020E+09 * EXP(  1.500*ALOGT -  8.60093*RTR )                  
      RB(11) = EXP(+( 1.86553E+07)*TM3+(-2.10973E+05)*TM2                       
     &         +(-2.43174E+03)*TM1+( 2.51339E+01)+( 1.84137E-03)*TP1            
     &         +(-2.56572E-07)*TP2+( 2.02667E-11)*TP3)                          
      RF(12) =  6.020E+14 * EXP(              -  3.00033*RTR )                  
      RB(12) = EXP(+( 1.96931E+10)*TM3+(-1.31710E+08)*TM2                       
     &         +( 2.66739E+05)*TM1+(-3.61551E+02)+( 2.47612E-01)*TP1            
     &         +(-7.43956E-05)*TP2+( 8.57843E-09)*TP3)                          
      RF(13) =  3.000E+13                                                       
      RB(13) = EXP(+( 7.28675E+06)*TM3+(-1.04626E+05)*TM2                       
     &         +(-4.30633E+04)*TM1+( 3.04502E+01)+( 8.12574E-04)*TP1            
     &         +(-1.74620E-07)*TP2+( 1.61236E-11)*TP3)                          
      RF(14) =  3.000E+13                                                       
      RB(14) = EXP(+( 7.36101E+09)*TM3+(-5.12866E+07)*TM2                       
     &         +( 7.75036E+04)*TM1+(-1.27971E+02)+( 1.03700E-01)*TP1            
     &         +(-3.16211E-05)*TP2+( 3.68798E-09)*TP3)                          
      RF(15) =  3.900E+13 * EXP(              -  3.54038*RTR )                  
      RB(15) = EXP(+( 1.76497E+07)*TM3+(-2.22300E+05)*TM2                       
     &         +(-7.95215E+03)*TM1+( 2.63311E+01)+( 6.90719E-04)*TP1            
     &         +(-1.07102E-07)*TP2+( 7.61156E-12)*TP3)                          
      RF(16) =  1.000E+13                                                       
      RB(16) = EXP(+(-1.80690E+07)*TM3+( 1.79596E+05)*TM2                       
     &         +(-3.71047E+04)*TM1+( 3.06705E+01)+( 4.61263E-04)*TP1            
     &         +(-1.02063E-07)*TP2+( 1.00951E-11)*TP3)                          
      RF(17) =  1.000E+13                                                       
      RB(17) = EXP(+( 5.48786E+06)*TM3+(-1.06202E+05)*TM2                       
     &         +(-3.96655E+04)*TM1+( 2.71413E+01)+( 8.53658E-04)*TP1            
     &         +(-1.73931E-07)*TP2+( 1.75364E-11)*TP3)                          
      RF(18) =  3.880E+05 * EXP(  2.500*ALOGT -  3.10034*RTR )                  
      RB(18) = EXP(+( 4.30756E+06)*TM3+(-6.40031E+04)*TM2                       
     &         +(-4.65836E+03)*TM1+( 2.46529E+01)+( 2.46874E-03)*TP1            
     &         +(-4.11511E-07)*TP2+( 3.42466E-11)*TP3)                          
      RF(19) =  1.300E+05 * EXP(  2.500*ALOGT -  5.00054*RTR )                  
      RB(19) = EXP(+(-1.92493E+07)*TM3+( 2.21794E+05)*TM2                       
     &         +(-3.05373E+03)*TM1+( 2.70887E+01)+( 2.07635E-03)*TP1            
     &         +(-3.39643E-07)*TP2+( 2.68053E-11)*TP3)                          
      RF(20) =  5.000E+13                                                       
      RB(20) = EXP(+(-7.92400E+06)*TM3+( 8.75558E+04)*TM2                       
     &         +(-3.97788E+04)*TM1+( 3.07605E+01)+( 7.51599E-04)*TP1            
     &         +(-1.15653E-07)*TP2+( 7.46281E-12)*TP3)                          
      RF(21) =  1.020E+07 * EXP(  2.000*ALOGT -  1.90021*RTR )                  
      RB(21) = EXP(+(-1.39594E+07)*TM3+( 1.78089E+05)*TM2                       
     &         +(-1.20504E+04)*TM1+( 2.99686E+01)+( 1.03735E-03)*TP1            
     &         +(-1.01471E-07)*TP2+( 5.78045E-12)*TP3)                          
      RF(22) =  4.600E+19 * EXP( -1.410*ALOGT - 28.95314*RTR )                  
      RB(22) = EXP(+( 8.50122E+06)*TM3+(-1.78101E+05)*TM2                       
     &         +( 2.14760E+03)*TM1+( 3.04179E+01)+(-2.66299E-04)*TP1            
     &         +(-1.60415E-09)*TP2+( 3.75199E-12)*TP3)                          
      RF(23) =  1.020E+07 * EXP(  2.000*ALOGT -  1.90021*RTR )                  
      RB(23) = EXP(+(-3.62183E+07)*TM3+( 3.65800E+05)*TM2                       
     &         +(-2.61334E+04)*TM1+( 2.60321E+01)+( 2.41379E-03)*TP1            
     &         +(-4.08009E-07)*TP2+( 3.63589E-11)*TP3)                          
      RF(24) =  3.000E+13                                                       
      RB(24) = EXP(+( 1.34855E+07)*TM3+(-1.61609E+05)*TM2                       
     &         +(-4.49385E+04)*TM1+( 3.50709E+01)+(-4.50422E-04)*TP1            
     &         +( 1.18121E-07)*TP2+(-1.22866E-11)*TP3)                          
      RF(25) =  1.920E+07 * EXP(  1.830*ALOGT -   .22002*RTR )                  
      RB(25) = EXP(+(-1.74594E+07)*TM3+( 1.23640E+05)*TM2                       
     &         +(-1.42993E+04)*TM1+( 2.33531E+01)+( 2.26025E-03)*TP1            
     &         +(-3.79476E-07)*TP2+( 3.25662E-11)*TP3)                          
      RF(26) =  1.320E+14                                                       
      RB(26) = EXP(+(-8.50127E+06)*TM3+( 3.86564E+04)*TM2                       
     &         +(-3.95483E+04)*TM1+( 3.13334E+01)+( 8.34648E-04)*TP1            
     &         +(-1.55034E-07)*TP2+( 1.41614E-11)*TP3)                          
      RF(27) =  8.980E+07 * EXP(  1.920*ALOGT -  5.69062*RTR )                  
      RB(27) = EXP(+(-2.39104E+06)*TM3+( 2.05147E+04)*TM2                       
     &         +(-4.20434E+03)*TM1+( 2.51615E+01)+( 2.16905E-03)*TP1            
     &         +(-3.70410E-07)*TP2+( 3.11994E-11)*TP3)                          
      RF(28) =  1.000E+14                                                       
      RB(28) = EXP(+( 2.67056E+09)*TM3+(-1.86048E+07)*TM2                       
     &         +(-3.31199E+03)*TM1+(-3.01608E+01)+( 3.91930E-02)*TP1            
     &         +(-1.17713E-05)*TP2+( 1.36680E-09)*TP3)                          
      RF(29) =  1.000E+13 * EXP(              -  8.00087*RTR )                  
      RB(29) = EXP(+( 1.24589E+07)*TM3+(-1.58599E+05)*TM2                       
     &         +(-1.48132E+03)*TM1+( 2.53621E+01)+( 7.34187E-04)*TP1            
     &         +(-1.16887E-07)*TP2+( 8.55950E-12)*TP3)                          
      RF(30) =  1.750E+12 * EXP(              -  1.35015*RTR )                  
      RB(30) = EXP(+(-3.17357E+06)*TM3+( 1.64860E+04)*TM2                       
     &         +(-2.50023E+04)*TM1+( 2.71186E+01)+( 7.28210E-04)*TP1            
     &         +(-1.34682E-07)*TP2+( 1.21702E-11)*TP3)                          
      RF(31) =  2.500E+12 * EXP(              - 47.80519*RTR )                  
      RB(31) = EXP(+( 9.42359E+06)*TM3+(-6.78517E+04)*TM2                       
     &         +(-2.80333E+04)*TM1+( 3.22157E+01)+(-6.90284E-04)*TP1            
     &         +( 1.17961E-07)*TP2+(-9.32107E-12)*TP3)                          
      RF(32) =  1.000E+14 * EXP(              - 40.00434*RTR )                  
      RB(32) = EXP(+( 1.89293E+07)*TM3+(-2.05230E+05)*TM2                       
     &         +( 2.65775E+02)*TM1+( 2.75723E+01)+( 3.61499E-04)*TP1            
     &         +(-5.23996E-08)*TP2+( 2.64655E-12)*TP3)                          
      RF(33) =  2.800E+18 * EXP(  -.860*ALOGT                )                  
      RB(33) = EXP(+( 2.33761E+06)*TM3+(-4.37788E+03)*TM2                       
     &         +(-2.42443E+04)*TM1+( 3.69366E+01)+(-5.01661E-04)*TP1            
     &         +( 1.62723E-08)*TP2+( 1.04204E-12)*TP3)                          
      RF(34) =  3.000E+20 * EXP( -1.720*ALOGT                )                  
      RB(34) = EXP(+( 1.07259E+07)*TM3+(-1.04678E+05)*TM2                       
     &         +(-2.36435E+04)*TM1+( 3.56837E+01)+(-1.12244E-03)*TP1            
     &         +( 1.24084E-07)*TP2+(-8.44059E-12)*TP3)                          
      RF(35) =  9.380E+18 * EXP(  -.760*ALOGT                )                  
      RB(35) = EXP(+( 1.36223E+06)*TM3+( 7.28493E+03)*TM2                       
     &         +(-2.43141E+04)*TM1+( 3.88348E+01)+(-4.29478E-04)*TP1            
     &         +( 3.73602E-09)*TP2+( 2.14467E-12)*TP3)                          
      RF(36) =  3.750E+20 * EXP( -1.720*ALOGT                )                  
      RB(36) = EXP(+( 1.07259E+07)*TM3+(-1.04678E+05)*TM2                       
     &         +(-2.36435E+04)*TM1+( 3.59068E+01)+(-1.12244E-03)*TP1            
     &         +( 1.24084E-07)*TP2+(-8.44059E-12)*TP3)                          
      RF(38) =  8.300E+13 * EXP(              - 14.41456*RTR )                  
      RB(38) = EXP(+( 2.79711E+06)*TM3+(-5.52261E+04)*TM2                       
     &         +( 1.55414E+03)*TM1+( 2.82823E+01)+( 6.92137E-04)*TP1            
     &         +(-1.70781E-07)*TP2+( 1.76467E-11)*TP3)                          
      RF(39) =  1.000E+18 * EXP( -1.000*ALOGT                )                  
      RB(39) = EXP(+( 4.78151E+06)*TM3+(-6.18811E+04)*TM2                       
     &         +(-5.16092E+04)*TM1+( 3.53257E+01)+(-2.23802E-04)*TP1            
     &         +(-2.30197E-08)*TP2+( 3.67249E-12)*TP3)                          
      RF(40) =  9.000E+16 * EXP(  -.600*ALOGT                )                  
      RB(40) = EXP(+( 8.79998E+05)*TM3+(-1.52298E+04)*TM2                       
     &         +(-5.18886E+04)*TM1+( 3.56745E+01)+( 6.49300E-05)*TP1            
     &         +(-7.31647E-08)*TP2+( 8.08301E-12)*TP3)                          
      RF(41) =  6.000E+19 * EXP( -1.250*ALOGT                )                  
      RB(41) = EXP(+( 7.21995E+06)*TM3+(-9.10381E+04)*TM2                       
     &         +(-5.14346E+04)*TM1+( 3.76970E+01)+(-4.04260E-04)*TP1            
     &         +( 8.32098E-09)*TP2+( 9.15911E-13)*TP3)                          
      RF(42) =  5.500E+20 * EXP( -2.000*ALOGT                )                  
      RB(42) = EXP(+( 1.45353E+07)*TM3+(-1.78509E+05)*TM2                       
     &         +(-5.09106E+04)*TM1+( 3.47436E+01)+(-9.45633E-04)*TP1            
     &         +( 1.02343E-07)*TP2+(-7.35382E-12)*TP3)                          
      RF(43) =  2.200E+22 * EXP( -2.000*ALOGT                )                  
      RB(43) = EXP(+( 1.25196E+10)*TM3+(-8.72153E+07)*TM2                       
     &         +( 1.68125E+05)*TM1+(-2.42014E+02)+( 1.76470E-01)*TP1            
     &         +(-5.39470E-05)*TP2+( 6.29427E-09)*TP3)                          
      RF(44) =  3.970E+12 * EXP(              -   .67107*RTR )                  
      RB(44) = EXP(+(-6.88870E+06)*TM3+( 3.39387E+04)*TM2                       
     &         +(-2.71224E+04)*TM1+( 2.77176E+01)+( 8.63184E-04)*TP1            
     &         +(-2.15857E-07)*TP2+( 2.24823E-11)*TP3)                          
      RF(45) =  2.800E+13 * EXP(              -  1.06812*RTR )                  
      RB(45) = EXP(+( 1.07837E+06)*TM3+(-4.11753E+04)*TM2                       
     &         +(-2.80002E+04)*TM1+( 3.13468E+01)+( 3.78915E-04)*TP1            
     &         +(-5.68427E-08)*TP2+( 4.17414E-12)*TP3)                          
      RF(46) =  1.340E+14 * EXP(              -   .63507*RTR )                  
      RB(46) = EXP(+( 1.51753E+06)*TM3+(-7.22963E+04)*TM2                       
     &         +(-1.80790E+04)*TM1+( 2.84618E+01)+( 1.02136E-03)*TP1            
     &         +(-2.25484E-07)*TP2+( 2.26117E-11)*TP3)                          
      RF(47) =  1.210E+07 * EXP(  2.000*ALOGT -  5.20056*RTR )                  
      RB(47) = EXP(+(-1.51440E+07)*TM3+( 1.70305E+05)*TM2                       
     &         +(-1.20478E+04)*TM1+( 2.80431E+01)+( 2.18100E-03)*TP1            
     &         +(-3.92328E-07)*TP2+( 3.28963E-11)*TP3)                          
      RF(48) =  1.000E+13 * EXP(              -  3.60039*RTR )                  
      RB(48) = EXP(+(-4.88311E+06)*TM3+(-4.90772E+03)*TM2                       
     &         +(-3.57346E+04)*TM1+( 2.59087E+01)+( 1.55083E-03)*TP1            
     &         +(-3.55320E-07)*TP2+( 3.41168E-11)*TP3)                          
      RF(49) =  1.100E+14                                                       
      RB(49) = EXP(+( 5.62116E+06)*TM3+(-5.30958E+04)*TM2                       
     &         +(-1.16826E+04)*TM1+( 3.31389E+01)+( 1.42451E-04)*TP1            
     &         +( 1.39284E-08)*TP2+(-3.94506E-12)*TP3)                          
      RF(51) =  3.000E+13                                                       
      RB(51) = EXP(+( 1.48026E+07)*TM3+(-1.70866E+05)*TM2                       
     &         +(-5.37812E+03)*TM1+( 2.86830E+01)+( 6.56277E-04)*TP1            
     &         +(-1.52159E-07)*TP2+( 1.25744E-11)*TP3)                          
      RF(53) =  6.600E+08 * EXP(  1.620*ALOGT - 10.84118*RTR )                  
      RB(53) = EXP(+( 1.98428E+07)*TM3+(-2.21082E+05)*TM2                       
     &         +(-4.53842E+03)*TM1+( 2.62087E+01)+( 1.97768E-03)*TP1            
     &         +(-2.73756E-07)*TP2+( 2.07990E-11)*TP3)                          
      RF(55) =  7.340E+13                                                       
      RB(55) = EXP(+( 9.64470E+06)*TM3+(-1.28731E+05)*TM2                       
     &         +(-4.39588E+04)*TM1+( 3.20280E+01)+( 8.62269E-04)*TP1            
     &         +(-1.76760E-07)*TP2+( 1.53327E-11)*TP3)                          
      RF(58) =  2.300E+10 * EXP(  1.050*ALOGT -  3.27536*RTR )                  
      RB(58) = EXP(+( 9.76617E+06)*TM3+(-1.23946E+05)*TM2                       
     &         +(-9.44779E+03)*TM1+( 2.68150E+01)+( 1.49834E-03)*TP1            
     &         +(-2.40873E-07)*TP2+( 1.83983E-11)*TP3)                          
      RF(60) =  2.000E+13                                                       
      RB(60) = EXP(+(-1.57111E+07)*TM3+( 1.55491E+05)*TM2                       
     &         +(-3.80002E+04)*TM1+( 3.20467E+01)+( 5.10958E-04)*TP1            
     &         +(-1.04203E-07)*TP2+( 9.30426E-12)*TP3)                          
      RF(61) =  1.200E+13                                                       
      RB(61) = EXP(+(-5.92436E+06)*TM3+( 2.34233E+04)*TM2                       
     &         +(-1.93737E+03)*TM1+( 2.70466E+01)+( 1.15335E-03)*TP1            
     &         +(-2.68819E-07)*TP2+( 2.62810E-11)*TP3)                          
      RF(62) =  6.000E+12                                                       
      RB(62) = EXP(+(-1.01057E+07)*TM3+( 5.50086E+04)*TM2                       
     &         +(-1.75797E+03)*TM1+( 2.61302E+01)+( 1.58498E-03)*TP1            
     &         +(-3.51071E-07)*TP2+( 3.32378E-11)*TP3)                          
      RF(64) =  3.400E+06 * EXP(  1.600*ALOGT                )                  
      RB(64) = EXP(+( 7.95084E+06)*TM3+(-9.91922E+04)*TM2                       
     &         +(-3.67854E+03)*TM1+( 2.25372E+01)+( 1.54732E-03)*TP1            
     &         +(-2.72448E-07)*TP2+( 2.50834E-11)*TP3)                          
      RF(65) =  2.000E+13                                                       
      RB(65) = EXP(+( 7.84581E+06)*TM3+(-1.30307E+05)*TM2                       
     &         +(-4.05610E+04)*TM1+( 2.85175E+01)+( 9.03353E-04)*TP1            
     &         +(-1.76072E-07)*TP2+( 1.67455E-11)*TP3)                          
      RF(66) =  3.200E+13                                                       
      RB(66) = EXP(+( 1.76325E+07)*TM3+(-2.62374E+05)*TM2                       
     &         +(-4.49821E+03)*TM1+( 2.44982E+01)+( 1.54574E-03)*TP1            
     &         +(-3.40688E-07)*TP2+( 3.37223E-11)*TP3)                          
      RF(67) =  1.600E+13                                                       
      RB(67) = EXP(+( 1.34511E+07)*TM3+(-2.30789E+05)*TM2                       
     &         +(-4.31882E+03)*TM1+( 2.35818E+01)+( 1.97738E-03)*TP1            
     &         +(-4.22939E-07)*TP2+( 4.06791E-11)*TP3)                          
      RF(68) =  1.700E+07 * EXP(  2.100*ALOGT -  4.87053*RTR )                  
      RB(68) = EXP(+( 1.05670E+07)*TM3+(-1.34759E+05)*TM2                       
     &         +(-6.16524E+03)*TM1+( 2.63592E+01)+( 2.22970E-03)*TP1            
     &         +(-3.63506E-07)*TP2+( 2.90452E-11)*TP3)                          
      RF(69) =  4.200E+06 * EXP(  2.100*ALOGT -  4.87053*RTR )                  
      RB(69) = EXP(+(-1.29899E+07)*TM3+( 1.51038E+05)*TM2                       
     &         +(-3.60440E+03)*TM1+( 2.84902E+01)+( 1.83731E-03)*TP1            
     &         +(-2.91638E-07)*TP2+( 2.16039E-11)*TP3)                          
      RF(73) =  3.000E+13                                                       
      RB(73) = EXP(+( 2.27542E+07)*TM3+(-2.89145E+05)*TM2                       
     &         +(-3.35921E+04)*TM1+( 3.11358E+01)+( 7.39776E-04)*TP1            
     &         +(-1.50161E-07)*TP2+( 1.17542E-11)*TP3)                          
      RF(75) =  1.325E+06 * EXP(  2.530*ALOGT - 12.24133*RTR )                  
      RB(75) = EXP(+(-5.19346E+06)*TM3+( 4.62776E+04)*TM2                       
     &         +(-3.33116E+03)*TM1+( 2.60762E+01)+( 2.63102E-03)*TP1            
     &         +(-4.43848E-07)*TP2+( 3.67619E-11)*TP3)                          
      RF(77) =  2.000E+12                                                       
      RB(77) = EXP(+( 1.11164E+07)*TM3+(-1.17959E+05)*TM2                       
     &         +(-3.37042E+04)*TM1+( 2.88930E+01)+( 6.35767E-04)*TP1            
     &         +(-1.14214E-07)*TP2+( 8.59401E-12)*TP3)                          
      RF(78) =  1.150E+08 * EXP(  1.900*ALOGT -  7.53082*RTR )                  
      RB(78) = EXP(+( 1.61981E+05)*TM3+(-5.92296E+03)*TM2                       
     &         +(-6.01191E+03)*TM1+( 2.59541E+01)+( 2.20431E-03)*TP1            
     &         +(-3.70043E-07)*TP2+( 3.01880E-11)*TP3)                          
      RF(79) =  1.000E+14                                                       
      RB(79) = EXP(+(-3.11684E+07)*TM3+( 2.75828E+05)*TM2                       
     &         +(-9.86681E+03)*TM1+( 2.92997E+01)+( 1.30281E-03)*TP1            
     &         +(-2.93255E-07)*TP2+( 2.96688E-11)*TP3)                          
      RF(80) =  5.000E+13 * EXP(              -  8.00087*RTR )                  
      RB(80) = EXP(+( 1.48168E+07)*TM3+(-1.82704E+05)*TM2                       
     &         +(-2.37684E+03)*TM1+( 2.76546E+01)+( 7.83882E-04)*TP1            
     &         +(-1.19028E-07)*TP2+( 7.76862E-12)*TP3)                          
      RF(81) =  1.130E+13 * EXP(              -  3.42837*RTR )                  
      RB(81) = EXP(+(-2.29343E+07)*TM3+( 1.91879E+05)*TM2                       
     &         +(-1.82519E+04)*TM1+( 2.55460E+01)+( 1.44720E-03)*TP1            
     &         +(-3.18264E-07)*TP2+( 3.11420E-11)*TP3)                          
      RF(82) =  1.000E+13                                                       
      RB(82) = EXP(+(-1.15433E+07)*TM3+( 1.46248E+05)*TM2                       
     &         +(-1.57368E+04)*TM1+( 2.98964E+01)+( 4.78691E-05)*TP1            
     &         +( 1.31837E-08)*TP2+(-1.74954E-12)*TP3)                          
      RF(84) =  2.160E+08 * EXP(  1.510*ALOGT -  3.43037*RTR )                  
      RB(84) = EXP(+(-2.54924E+07)*TM3+( 3.06449E+05)*TM2                       
     &         +(-1.09109E+04)*TM1+( 3.16895E+01)+( 8.82096E-04)*TP1            
     &         +(-1.77530E-07)*TP2+( 1.73112E-11)*TP3)                          
      RF(86) =  3.570E+04 * EXP(  2.400*ALOGT +  2.11023*RTR )                  
      RB(86) = EXP(+(-3.18153E+07)*TM3+( 3.86143E+05)*TM2                       
     &         +(-9.63997E+03)*TM1+( 2.97985E+01)+( 1.57422E-03)*TP1            
     &         +(-2.91243E-07)*TP2+( 2.63337E-11)*TP3)                          
      RF(87) =  2.900E+13 * EXP(              +   .50005*RTR )                  
      RB(87) = EXP(+(-9.68582E+06)*TM3+( 8.91648E+04)*TM2                       
     &         +(-3.53409E+04)*TM1+( 3.34736E+01)+( 1.71047E-04)*TP1            
     &         +(-4.50753E-08)*TP2+( 4.83558E-12)*TP3)                          
      RF(88) =  1.750E+12 * EXP(              -   .32003*RTR )                  
      RB(88) = EXP(+(-6.40064E+06)*TM3+( 6.73886E+04)*TM2                       
     &         +(-1.63245E+04)*TM1+( 2.82328E+01)+( 5.29473E-04)*TP1            
     &         +(-1.29836E-07)*TP2+( 1.15051E-11)*TP3)                          
      RF(89) =  5.800E+14 * EXP(              -  9.56104*RTR )                  
      RB(89) = EXP(+(-6.40064E+06)*TM3+( 6.73886E+04)*TM2                       
     &         +(-2.09747E+04)*TM1+( 3.40362E+01)+( 5.29473E-04)*TP1            
     &         +(-1.29836E-07)*TP2+( 1.15051E-11)*TP3)                          
      RF(90) =  5.000E+13                                                       
      RB(90) = EXP(+(-3.24478E+09)*TM3+( 6.90299E+07)*TM2                       
     &         +(-3.43763E+05)*TM1+( 4.40613E+02)+(-2.88005E-01)*TP1            
     &         +( 9.47620E-05)*TP2+(-1.16471E-08)*TP3)                          
      RF(91) =  3.000E+13                                                       
      RB(91) = EXP(+(-1.11645E+07)*TM3+( 1.55282E+05)*TM2                       
     &         +(-4.60600E+04)*TM1+( 3.59054E+01)+(-1.19775E-03)*TP1            
     &         +( 3.01551E-07)*TP2+(-2.90258E-11)*TP3)                          
      RF(92) =  2.000E+13                                                       
      RB(92) = EXP(+(-2.52790E+07)*TM3+( 3.18940E+05)*TM2                       
     &         +(-4.01559E+04)*TM1+( 3.84311E+01)+(-1.35552E-03)*TP1            
     &         +( 2.71918E-07)*TP2+(-2.41817E-11)*TP3)                          
      RF(93) =  1.130E+07 * EXP(  2.000*ALOGT -  3.00033*RTR )                  
      RB(93) = EXP(+(-2.43786E+07)*TM3+( 2.80849E+05)*TM2                       
     &         +(-1.21987E+04)*TM1+( 3.07668E+01)+( 1.81844E-03)*TP1            
     &         +(-3.77833E-07)*TP2+( 3.43788E-11)*TP3)                          
      RF(94) =  3.000E+13                                                       
      RB(94) = EXP(+(-1.63695E+07)*TM3+( 2.30822E+05)*TM2                       
     &         +(-4.43720E+04)*TM1+( 3.78366E+01)+(-1.28189E-03)*TP1            
     &         +( 2.58635E-07)*TP2+(-2.32721E-11)*TP3)                          
      RF(96) =  5.600E+07 * EXP(  1.600*ALOGT -  5.42059*RTR )                  
      RB(96) = EXP(+(-1.08780E+07)*TM3+( 1.30072E+05)*TM2                       
     &         +(-7.88219E+03)*TM1+( 2.76448E+01)+( 1.66019E-03)*TP1            
     &         +(-2.96115E-07)*TP2+( 2.55085E-11)*TP3)                          
      RF(97) =  2.501E+13                                                       
      RB(97) = EXP(+(-4.18137E+06)*TM3+( 3.15853E+04)*TM2                       
     &         +( 1.79394E+02)*TM1+( 3.06270E+01)+( 4.31632E-04)*TP1            
     &         +(-8.22514E-08)*TP2+( 6.95683E-12)*TP3)                          
      RF(98) =  1.000E+08 * EXP(  1.600*ALOGT -  3.12034*RTR )                  
      RB(98) = EXP(+( 9.27365E+06)*TM3+(-9.30748E+04)*TM2                       
     &         +(-8.76901E+03)*TM1+( 2.62756E+01)+( 1.75538E-03)*TP1            
     &         +(-2.59481E-07)*TP2+( 2.12399E-11)*TP3)                          
      RF(99) =  4.760E+07 * EXP(  1.228*ALOGT -   .07001*RTR )                  
      RB(99) = EXP(+(-5.35116E+06)*TM3+( 1.30594E+05)*TM2                       
     &         +(-1.36778E+04)*TM1+( 3.35776E+01)+(-4.96013E-04)*TP1            
     &         +( 1.34797E-07)*TP2+(-1.34274E-11)*TP3)                          
      RF(100) =  5.000E+13                                                      
      RB(100) = EXP(+( 1.04729E+09)*TM3+(-7.30012E+06)*TM2                      
     &          +(-3.30597E+04)*TM1+( 1.00685E+01)+( 1.55334E-02)*TP1           
     &          +(-4.69353E-06)*TP2+( 5.43807E-10)*TP3)                         
      RF(101) =  3.430E+09 * EXP(  1.180*ALOGT +   .44705*RTR )                 
      RB(101) = EXP(+(-2.26600E+06)*TM3+( 2.15558E+04)*TM2                      
     &          +(-1.57952E+04)*TM1+( 2.78998E+01)+( 1.38431E-03)*TP1           
     &          +(-2.45403E-07)*TP2+( 2.04932E-11)*TP3)                         
      RF(102) =  5.000E+12                                                      
      RB(102) = EXP(+(-2.64752E+07)*TM3+( 2.85831E+05)*TM2                      
     &          +(-4.61300E+04)*TM1+( 3.27522E+01)+( 3.03090E-04)*TP1           
     &          +(-9.24361E-08)*TP2+( 9.96570E-12)*TP3)                         
      RF(103) =  5.000E+12                                                      
      RB(103) = EXP(+(-2.91837E+06)*TM3+( 3.33709E+01)*TM2                      
     &          +(-4.86909E+04)*TM1+( 2.92230E+01)+( 6.95485E-04)*TP1           
     &          +(-1.64304E-07)*TP2+( 1.74070E-11)*TP3)                         
      RF(104) =  1.440E+06 * EXP(  2.000*ALOGT +   .84009*RTR )                 
      RB(104) = EXP(+( 7.78212E+05)*TM3+(-1.60821E+04)*TM2                      
     &          +(-1.13515E+04)*TM1+( 2.52932E+01)+( 1.94965E-03)*TP1           
     &          +(-3.39203E-07)*TP2+( 2.86040E-11)*TP3)                         
      RF(105) =  6.300E+06 * EXP(  2.000*ALOGT -  1.50016*RTR )                 
      RB(105) = EXP(+(-2.27787E+07)*TM3+( 2.69715E+05)*TM2                      
     &          +(-9.96833E+03)*TM1+( 3.02983E+01)+( 1.55726E-03)*TP1           
     &          +(-2.67334E-07)*TP2+( 2.11627E-11)*TP3)                         
      RF(106) =  2.000E+13                                                      
      RB(106) = EXP(+( 1.07997E+07)*TM3+(-4.15119E+04)*TM2                      
     &          +(-2.54294E+04)*TM1+( 3.58130E+01)+(-1.15780E-03)*TP1           
     &          +( 3.27620E-07)*TP2+(-3.55713E-11)*TP3)                         
      RF(107) =  2.180E-04 * EXP(  4.500*ALOGT +  1.00011*RTR )                 
      RB(107) = EXP(+(-5.08027E+07)*TM3+( 6.28259E+05)*TM2                      
     &          +(-1.48822E+04)*TM1+( 2.72011E+01)+( 2.10774E-03)*TP1           
     &          +(-2.97990E-07)*TP2+( 2.47867E-11)*TP3)                         
      RF(108) =  5.040E+05 * EXP(  2.300*ALOGT - 13.50147*RTR )                 
      RB(108) = EXP(+(-1.78011E+07)*TM3+( 2.25429E+05)*TM2                      
     &          +(-4.90595E+03)*TM1+( 3.36373E+01)+( 4.71839E-04)*TP1           
     &          +(-3.53759E-08)*TP2+( 2.27839E-12)*TP3)                         
      RF(109) =  3.370E+07 * EXP(  2.000*ALOGT - 14.00152*RTR )                 
      RB(109) = EXP(+(-3.31654E+07)*TM3+( 3.25836E+05)*TM2                      
     &          +(-1.73592E+03)*TM1+( 2.87523E+01)+( 2.03697E-03)*TP1           
     &          +(-4.19464E-07)*TP2+( 4.12223E-11)*TP3)                         
      RF(110) =  4.830E-04 * EXP(  4.000*ALOGT +  2.00022*RTR )                 
      RB(110) = EXP(+(-6.88602E+07)*TM3+( 7.61824E+05)*TM2                      
     &          +(-3.05563E+04)*TM1+( 2.00408E+01)+( 3.19402E-03)*TP1           
     &          +(-5.53572E-07)*TP2+( 5.04156E-11)*TP3)                         
      RF(111) =  5.000E+12                                                      
      RB(111) = EXP(+( 1.19900E+07)*TM3+(-1.58805E+05)*TM2                      
     &          +(-4.17219E+04)*TM1+( 3.14358E+01)+( 5.31907E-04)*TP1           
     &          +(-1.38393E-07)*TP2+( 1.24157E-11)*TP3)                         
      RF(112) =  3.600E+06 * EXP(  2.000*ALOGT -  2.50027*RTR )                 
      RB(112) = EXP(+(-1.07881E+07)*TM3+( 1.14805E+05)*TM2                      
     &          +(-6.18888E+03)*TM1+( 2.55148E+01)+( 2.04058E-03)*TP1           
     &          +(-3.65638E-07)*TP2+( 3.15794E-11)*TP3)                         
      RF(113) =  3.540E+06 * EXP(  2.120*ALOGT -   .87009*RTR )                 
      RB(113) = EXP(+(-1.27480E+07)*TM3+( 1.50075E+05)*TM2                      
     &          +(-1.09436E+04)*TM1+( 2.60813E+01)+( 2.15524E-03)*TP1           
     &          +(-3.85855E-07)*TP2+( 3.32753E-11)*TP3)                         
      RF(114) =  7.500E+12 * EXP(              -  2.00022*RTR )                 
      RB(114) = EXP(+( 4.05265E+06)*TM3+(-5.23637E+04)*TM2                      
     &          +(-7.48703E+03)*TM1+( 2.78493E+01)+( 5.76014E-04)*TP1           
     &          +(-1.07260E-07)*TP2+( 8.43006E-12)*TP3)                         
      RF(115) =  1.300E+11 * EXP(              +  1.63018*RTR )                 
      RB(115) = EXP(+(-3.28517E+06)*TM3+( 2.17763E+04)*TM2                      
     &          +(-1.86088E+04)*TM1+( 2.80240E+01)+(-3.58426E-04)*TP1           
     &          +( 8.47605E-08)*TP2+(-6.66952E-12)*TP3)                         
      RF(116) =  4.200E+14 * EXP(              - 12.00130*RTR )                 
      RB(116) = EXP(+(-3.28517E+06)*TM3+( 2.17763E+04)*TM2                      
     &          +(-2.54683E+04)*TM1+( 3.61044E+01)+(-3.58426E-04)*TP1           
     &          +( 8.47605E-08)*TP2+(-6.66952E-12)*TP3)                         
      RF(117) =  2.000E+13                                                      
      RB(117) = EXP(+( 1.24347E+10)*TM3+(-8.65210E+07)*TM2                      
     &          +( 1.68209E+05)*TM1+(-2.46879E+02)+( 1.76475E-01)*TP1           
     &          +(-5.37669E-05)*TP2+( 6.27052E-09)*TP3)                         
      RF(118) =  1.000E+12                                                      
      RB(118) = EXP(+(-3.45655E+07)*TM3+( 3.68845E+05)*TM2                      
     &          +(-2.95114E+04)*TM1+( 3.32786E+01)+(-4.29403E-04)*TP1           
     &          +( 1.38259E-08)*TP2+( 1.23780E-12)*TP3)                         
      RF(119) =  2.000E+13                                                      
      RB(119) = EXP(+(-1.61150E+07)*TM3+( 1.90078E+05)*TM2                      
     &          +(-1.32612E+04)*TM1+( 3.31582E+01)+(-5.24386E-04)*TP1           
     &          +( 1.15204E-07)*TP2+(-1.11106E-11)*TP3)                         
      RF(120) =  1.500E+14 * EXP(              - 23.60256*RTR )                 
      RB(120) = EXP(+( 8.14400E+06)*TM3+(-8.49219E+04)*TM2                      
     &          +(-4.24214E+04)*TM1+( 3.60105E+01)+(-3.61064E-04)*TP1           
     &          +( 6.32591E-08)*TP2+(-4.35605E-12)*TP3)                         
      RF(121) =  1.000E+12 * EXP(              -  8.00087*RTR )                 
      RB(121) = EXP(+( 1.56441E+07)*TM3+(-1.83454E+05)*TM2                      
     &          +(-3.05863E+03)*TM1+( 2.54003E+01)+( 3.07240E-06)*TP1           
     &          +( 3.23609E-08)*TP2+(-4.02297E-12)*TP3)                         
      RF(122) =  5.800E+13 * EXP(              -   .57606*RTR )                 
      RB(122) = EXP(+( 1.22142E+10)*TM3+(-5.96839E+07)*TM2                      
     &          +( 3.43789E+04)*TM1+(-4.91468E+01)+( 2.86505E-02)*TP1           
     &          +(-3.58138E-06)*TP2+(-3.31173E-11)*TP3)                         
      RF(123) =  5.000E+13                                                      
      RB(123) = EXP(+( 4.31826E+06)*TM3+(-6.65517E+04)*TM2                      
     &          +(-3.88238E+04)*TM1+( 3.44601E+01)+(-6.96583E-04)*TP1           
     &          +( 8.97804E-08)*TP2+(-4.75515E-12)*TP3)                         
      RF(124) =  5.000E+13                                                      
      RB(124) = EXP(+( 2.27042E+07)*TM3+(-2.15664E+05)*TM2                      
     &          +(-4.95676E+04)*TM1+( 3.56016E+01)+(-7.84630E-04)*TP1           
     &          +( 1.62984E-07)*TP2+(-1.60584E-11)*TP3)                         
      RF(125) =  3.300E+13                                                      
      RB(125) = EXP(+(-8.36735E+06)*TM3+( 1.00056E+05)*TM2                      
     &          +(-3.72522E+04)*TM1+( 3.22332E+01)+(-5.05617E-04)*TP1           
     &          +( 1.30770E-07)*TP2+(-1.13791E-11)*TP3)                         
      RF(126) =  1.107E+08 * EXP(  1.790*ALOGT -  1.67018*RTR )                 
      RB(126) = EXP(+(-2.33524E+07)*TM3+( 2.91512E+05)*TM2                      
     &          +(-9.28919E+02)*TM1+( 3.22083E+01)+( 7.09430E-04)*TP1           
     &          +(-8.55237E-08)*TP2+( 8.07232E-12)*TP3)                         
      RF(127) =  5.710E+12 * EXP(              +   .75508*RTR )                 
      RB(127) = EXP(+(-2.04079E+07)*TM3+( 2.71348E+05)*TM2                      
     &          +(-3.04841E+04)*TM1+( 3.64350E+01)+(-1.73030E-03)*TP1           
     &          +( 3.99026E-07)*TP2+(-3.65079E-11)*TP3)                         
      RF(128) =  4.000E+13                                                      
      RB(128) = EXP(+( 1.89203E+10)*TM3+(-1.24460E+08)*TM2                      
     &          +( 2.43907E+05)*TM1+(-3.31051E+02)+( 2.24984E-01)*TP1           
     &          +(-6.70994E-05)*TP2+( 7.69394E-09)*TP3)                         
      RF(129) =  3.000E+13                                                      
      RB(129) = EXP(+( 5.57118E+06)*TM3+( 2.03853E+04)*TM2                      
     &          +(-2.76581E+04)*TM1+( 3.57946E+01)+(-1.38196E-03)*TP1           
     &          +( 3.27073E-07)*TP2+(-3.17576E-11)*TP3)                         
      RF(130) =  6.000E+13                                                      
      RB(130) = EXP(+( 2.17315E+07)*TM3+(-1.40843E+05)*TM2                      
     &          +(-3.02056E+04)*TM1+( 3.66811E+01)+(-1.37843E-03)*TP1           
     &          +( 3.83085E-07)*TP2+(-3.76867E-11)*TP3)                         
      RF(132) =  3.400E+12 * EXP(              -   .69007*RTR )                 
      RB(132) = EXP(+(-1.77909E+07)*TM3+( 1.67908E+05)*TM2                      
     &          +(-3.36225E+04)*TM1+( 2.62921E+01)+( 1.84667E-04)*TP1           
     &          +( 1.28082E-08)*TP2+(-2.05805E-12)*TP3)                         
      RF(133) =  9.460E+13 * EXP(              +   .51506*RTR )                 
      RB(133) = EXP(+( 3.12013E+07)*TM3+(-2.97396E+05)*TM2                      
     &          +(-3.71701E+04)*TM1+( 3.71755E+01)+(-1.14029E-03)*TP1           
     &          +( 2.78439E-07)*TP2+(-2.78584E-11)*TP3)                         
      RF(134) =  5.000E+13                                                      
      RB(134) = EXP(+(-5.31373E+09)*TM3+( 8.62306E+07)*TM2                      
     &          +(-3.94150E+05)*TM1+( 5.05494E+02)+(-3.29611E-01)*TP1           
     &          +( 1.07708E-04)*TP2+(-1.31787E-08)*TP3)                         
      RF(135) =  1.320E+13 * EXP(              -  1.50016*RTR )                 
      RB(135) = EXP(+(-4.83217E+06)*TM3+( 4.14134E+04)*TM2                      
     &          +(-3.82736E+04)*TM1+( 2.92846E+01)+( 2.73353E-05)*TP1           
     &          +(-5.96537E-09)*TP2+( 1.07654E-12)*TP3)                         
      RF(136) =  5.000E+05 * EXP(  2.000*ALOGT -  7.23078*RTR )                 
      RB(136) = EXP(+(-3.49998E+07)*TM3+( 4.20129E+05)*TM2                      
     &          +(-9.12882E+03)*TM1+( 3.02213E+01)+( 7.30532E-04)*TP1           
     &          +(-1.43423E-07)*TP2+( 1.48476E-11)*TP3)                         
      RF(137) =  3.200E+13                                                      
      RB(137) = EXP(+( 1.73347E+10)*TM3+(-1.08558E+08)*TM2                      
     &          +( 1.91092E+05)*TM1+(-2.58931E+02)+( 1.74799E-01)*TP1           
     &          +(-5.07768E-05)*TP2+( 5.69799E-09)*TP3)                         
      RF(138) =  4.000E+13                                                      
      RB(138) = EXP(+(-8.01929E+06)*TM3+( 1.86429E+05)*TM2                      
     &          +(-3.34163E+04)*TM1+( 4.01904E+01)+(-1.60410E-03)*TP1           
     &          +( 3.14879E-07)*TP2+(-2.89582E-11)*TP3)                         
      RF(139) =  2.460E+06 * EXP(  2.000*ALOGT -  8.27090*RTR )                 
      RB(139) = EXP(+( 6.44065E+05)*TM3+( 1.01091E+04)*TM2                      
     &          +(-7.60352E+03)*TM1+( 2.65506E+01)+( 1.53885E-03)*TP1           
     &          +(-2.14092E-07)*TP2+( 1.77840E-11)*TP3)                         
      RF(141) =  3.000E+13                                                      
      RB(141) = EXP(+(-2.62869E+07)*TM3+( 3.12221E+05)*TM2                      
     &          +(-4.69961E+04)*TM1+( 3.38238E+01)+(-1.35994E-04)*TP1           
     &          +(-1.10378E-08)*TP2+( 3.28056E-12)*TP3)                         
      RF(142) =  1.500E+13 * EXP(              -   .60007*RTR )                 
      RB(142) = EXP(+( 8.90944E+06)*TM3+(-8.81179E+04)*TM2                      
     &          +(-4.51812E+03)*TM1+( 2.93391E+01)+( 7.36294E-05)*TP1           
     &          +(-1.32833E-08)*TP2+( 9.09603E-13)*TP3)                         
      RF(144) =  2.800E+13                                                      
      RB(144) = EXP(+( 1.86942E+07)*TM3+(-2.30183E+05)*TM2                      
     &          +(-3.33859E+04)*TM1+( 2.83665E+01)+( 4.65205E-04)*TP1           
     &          +(-4.76266E-08)*TP2+( 2.62009E-12)*TP3)                         
      RF(145) =  1.200E+13                                                      
      RB(145) = EXP(+(-3.54045E+10)*TM3+( 3.13634E+08)*TM2                      
     &          +(-9.86040E+05)*TM1+( 1.17239E+03)+(-7.20529E-01)*TP1           
     &          +( 2.18278E-04)*TP2+(-2.52372E-08)*TP3)                         
      RF(146) =  7.000E+13                                                      
      RB(146) = EXP(+(-6.58281E+06)*TM3+( 9.87548E+04)*TM2                      
     &          +(-8.30921E+03)*TM1+( 3.41946E+01)+(-6.39501E-04)*TP1           
     &          +( 9.40188E-08)*TP2+(-6.29539E-12)*TP3)                         
      RF(148) =  3.000E+13                                                      
      RB(148) = EXP(+( 8.90944E+06)*TM3+(-8.81179E+04)*TM2                      
     &          +(-4.21616E+03)*TM1+( 3.00323E+01)+( 7.36294E-05)*TP1           
     &          +(-1.32833E-08)*TP2+( 9.09603E-13)*TP3)                         
      RF(149) =  1.200E+13 * EXP(              +   .57006*RTR )                 
      RB(149) = EXP(+( 8.90157E+05)*TM3+( 9.83113E+04)*TM2                      
     &          +(-3.73456E+04)*TM1+( 3.79865E+01)+(-1.53047E-03)*TP1           
     &          +( 3.01595E-07)*TP2+(-2.80486E-11)*TP3)                         
      RF(150) =  1.600E+13 * EXP(              +   .57006*RTR )                 
      RB(150) = EXP(+( 2.90611E+07)*TM3+(-3.11265E+05)*TM2                      
     &          +(-5.97364E+03)*TM1+( 2.74546E+01)+( 1.68818E-04)*TP1           
     &          +( 2.33502E-08)*TP2+(-3.35905E-12)*TP3)                         
      RF(151) =  9.000E+12                                                      
      RB(151) = EXP(+( 8.90944E+06)*TM3+(-8.81179E+04)*TM2                      
     &          +(-4.21616E+03)*TM1+( 2.88283E+01)+( 7.36294E-05)*TP1           
     &          +(-1.32833E-08)*TP2+( 9.09603E-13)*TP3)                         
      RF(152) =  7.000E+12                                                      
      RB(152) = EXP(+( 8.90944E+06)*TM3+(-8.81179E+04)*TM2                      
     &          +(-4.21616E+03)*TM1+( 2.85770E+01)+( 7.36294E-05)*TP1           
     &          +(-1.32833E-08)*TP2+( 9.09603E-13)*TP3)                         
      RF(153) =  1.400E+13                                                      
      RB(153) = EXP(+(-2.29960E+07)*TM3+( 2.43448E+05)*TM2                      
     &          +(-3.15873E+04)*TM1+( 2.96385E+01)+( 1.00530E-04)*TP1           
     &          +(-3.01081E-08)*TP2+( 3.69565E-12)*TP3)                         
      RF(154) =  4.000E+13 * EXP(              +   .55006*RTR )                 
      RB(154) = EXP(+( 1.21113E+07)*TM3+(-1.28762E+05)*TM2                      
     &          +(-8.92743E+03)*TM1+( 2.79339E+01)+( 1.93331E-04)*TP1           
     &          +(-3.78352E-08)*TP2+( 2.94265E-12)*TP3)                         
      RF(155) =  2.675E+13 * EXP(              - 28.80313*RTR )                 
      RB(155) = EXP(+(-1.48354E+07)*TM3+( 2.07148E+05)*TM2                      
     &          +(-1.18821E+03)*TM1+( 3.37485E+01)+(-8.53606E-04)*TP1           
     &          +( 1.69906E-07)*TP2+(-1.60756E-11)*TP3)                         
      RF(156) =  3.600E+10 * EXP(              -  8.94097*RTR )                 
      RB(156) = EXP(+(-9.34754E+06)*TM3+( 1.00946E+05)*TM2                      
     &          +(-3.08588E+04)*TM1+( 2.43455E+01)+( 5.15687E-08)*TP1           
     &          +(-4.02517E-09)*TP2+( 1.46085E-12)*TP3)                         
      RF(157) =  2.450E+04 * EXP(  2.470*ALOGT -  5.18056*RTR )                 
      RB(157) = EXP(+(-5.53721E+07)*TM3+( 6.35140E+05)*TM2                      
     &          +(-1.44147E+04)*TM1+( 3.03441E+01)+( 1.71195E-03)*TP1           
     &          +(-3.80580E-07)*TP2+( 3.51423E-11)*TP3)                         
      RF(159) =  4.990E+12 * EXP(   .100*ALOGT - 10.60115*RTR )                 
      RB(159) = EXP(+(-4.61876E+06)*TM3+( 1.29179E+05)*TM2                      
     &          +(-1.02355E+03)*TM1+( 3.49143E+01)+(-1.45455E-03)*TP1           
     &          +( 3.09254E-07)*TP2+(-2.92446E-11)*TP3)                         
      RF(160) =  2.648E+13                                                      
      RB(160) = EXP(+(-2.59992E+07)*TM3+( 2.81289E+05)*TM2                      
     &          +(-4.60075E+04)*TM1+( 3.62725E+01)+( 5.39505E-05)*TP1           
     &          +(-1.06092E-07)*TP2+( 1.23964E-11)*TP3)                         
      RF(161) =  3.320E+03 * EXP(  2.810*ALOGT -  5.86064*RTR )                 
      RB(161) = EXP(+(-4.30443E+07)*TM3+( 4.91340E+05)*TM2                      
     &          +(-1.40269E+04)*TM1+( 2.84579E+01)+( 1.96044E-03)*TP1           
     &          +(-3.90843E-07)*TP2+( 3.48683E-11)*TP3)                         
      RF(162) =  3.000E+07 * EXP(  1.500*ALOGT -  9.94108*RTR )                 
      RB(162) = EXP(+(-1.92246E+07)*TM3+( 2.05284E+05)*TM2                      
     &          +(-1.03464E+04)*TM1+( 2.80560E+01)+( 9.88287E-04)*TP1           
     &          +(-2.17620E-07)*TP2+( 1.94931E-11)*TP3)                         
      RF(163) =  1.000E+07 * EXP(  1.500*ALOGT -  9.94108*RTR )                 
      RB(163) = EXP(+(-4.27815E+07)*TM3+( 4.91081E+05)*TM2                      
     &          +(-7.78556E+03)*TM1+( 3.04866E+01)+( 5.95892E-04)*TP1           
     &          +(-1.45752E-07)*TP2+( 1.20518E-11)*TP3)                         
      RF(164) =  2.270E+05 * EXP(  2.000*ALOGT -  9.20100*RTR )                 
      RB(164) = EXP(+(-3.56678E+07)*TM3+( 3.94485E+05)*TM2                      
     &          +(-3.47969E+03)*TM1+( 2.59233E+01)+( 1.44013E-03)*TP1           
     &          +(-3.06737E-07)*TP2+( 2.79817E-11)*TP3)                         
      RF(165) =  6.140E+06 * EXP(  1.740*ALOGT - 10.45113*RTR )                 
      RB(165) = EXP(+(-3.39213E+07)*TM3+( 3.85436E+05)*TM2                      
     &          +(-9.41840E+03)*TM1+( 2.71854E+01)+( 1.28050E-03)*TP1           
     &          +(-2.79316E-07)*TP2+( 2.54875E-11)*TP3)                         
      RF(166) =  2.244E+18 * EXP( -1.000*ALOGT - 17.00185*RTR )                 
      RB(166) = EXP(+( 2.43707E+07)*TM3+(-3.00106E+05)*TM2                      
     &          +( 4.91899E+02)*TM1+( 3.46927E+01)+(-3.57591E-04)*TP1           
     &          +( 9.69846E-08)*TP2+(-1.03924E-11)*TP3)                         
      RF(167) =  1.870E+17 * EXP( -1.000*ALOGT - 17.00185*RTR )                 
      RB(167) = EXP(+( 2.43707E+07)*TM3+(-3.00106E+05)*TM2                      
     &          +( 4.91899E+02)*TM1+( 3.22078E+01)+(-3.57591E-04)*TP1           
     &          +( 9.69846E-08)*TP2+(-1.03924E-11)*TP3)                         
      RF(168) =  7.600E+12 * EXP(              -   .40004*RTR )                 
      RB(168) = EXP(+( 8.56633E+06)*TM3+(-8.75559E+04)*TM2                      
     &          +(-1.66974E+04)*TM1+( 2.93767E+01)+( 4.83354E-04)*TP1           
     &          +(-1.19918E-07)*TP2+( 1.11586E-11)*TP3)                         
      RF(169) =  1.800E+13 * EXP(              -   .90010*RTR )                 
      RB(169) = EXP(+(-1.67894E+07)*TM3+( 1.96666E+05)*TM2                      
     &          +(-1.09904E+04)*TM1+( 3.15578E+01)+( 1.32043E-04)*TP1           
     &          +(-4.73608E-08)*TP2+( 5.13012E-12)*TP3)                         
      RF(170) =  4.280E-13 * EXP(  7.600*ALOGT +  3.53038*RTR )                 
      RB(170) = EXP(+(-6.73612E+07)*TM3+( 7.97243E+05)*TM2                      
     &          +(-1.66308E+04)*TM1+( 2.14066E+01)+( 6.01035E-03)*TP1           
     &          +(-1.07198E-06)*TP2+( 9.63714E-11)*TP3)                         
      RF(171) =  5.000E+13 * EXP(              -  1.50016*RTR )                 
      RB(171) = EXP(+(-5.10208E+09)*TM3+( 8.35052E+07)*TM2                      
     &          +(-3.82863E+05)*TM1+( 4.87577E+02)+(-3.19719E-01)*TP1           
     &          +( 1.04617E-04)*TP2+(-1.28058E-08)*TP3)                         
      RF(172) =  4.070E+05 * EXP(  2.400*ALOGT -   .20002*RTR )                 
      RB(172) = EXP(+(-2.05154E+07)*TM3+( 3.17668E+05)*TM2                      
     &          +(-1.66140E+04)*TM1+( 3.39138E+01)+( 9.31217E-04)*TP1           
     &          +(-1.20364E-07)*TP2+( 7.95494E-12)*TP3)                         
      RF(173) =  3.980E+12 * EXP(              +   .24003*RTR )                 
      RB(173) = EXP(+(-2.60832E+07)*TM3+( 2.35843E+05)*TM2                      
     &          +(-4.46407E+04)*TM1+( 2.91618E+01)+( 1.84252E-04)*TP1           
     &          +(-2.95475E-08)*TP2+( 4.19268E-12)*TP3)                         
      RF(175) =  8.400E+11 * EXP(              -  3.87542*RTR )                 
      RB(175) = EXP(+( 1.00380E+07)*TM3+(-7.67841E+04)*TM2                      
     &          +(-8.19167E+03)*TM1+( 2.76420E+01)+( 2.56852E-04)*TP1           
     &          +(-5.73709E-08)*TP2+( 4.41987E-12)*TP3)                         
      RF(176) =  1.600E+12 * EXP(              -   .85409*RTR )                 
      RB(176) = EXP(+(-1.24741E+07)*TM3+( 4.56455E+04)*TM2                      
     &          +(-4.36825E+04)*TM1+( 2.25678E+01)+( 1.76802E-03)*TP1           
     &          +(-3.40881E-07)*TP2+( 3.22889E-11)*TP3)                         
      RF(177) =  1.000E+13                                                      
      RB(177) = EXP(+(-2.08194E+07)*TM3+( 1.56039E+05)*TM2                      
     &          +(-4.23634E+04)*TM1+( 2.81212E+01)+( 1.48220E-03)*TP1           
     &          +(-3.19354E-07)*TP2+( 3.09144E-11)*TP3)                         
      TOLD=T
      END IF
      RFH50 =  2.500E+16 * EXP(  -.800*ALOGT                )                   
      RFL50 =  3.200E+27 * EXP( -3.140*ALOGT -  1.23013*RTR )                   
      PR = RFL50*XM(50)/RFH50                                                   
      PRLOG = DLOG10(MAX(PR,SMALL))
      F4 =     6.80000E-01
      F5 =     7.80000E+01
      F6 =     1.99500E+03
      F7 =     5.59000E+03
      FC = (1.0D0-F4)*EXP(-T/F5)+F4*EXP(-T/F6)+EXP(-F7/T)
      FCLOG = DLOG10(MAX(FC,SMALL))
      CPRLOG = PRLOG - (0.4D0 + 0.67D0*FCLOG)
      X = CPRLOG/(0.75D0-1.27D0*FCLOG-0.14D0*CPRLOG)
      F = 10.0D0**( FCLOG/(1.0D0+X*X) )
      PCOR = PR*F/(1.0D0+PR)
      RF(50) = RFH50*PCOR                                                       
      RB(50) = EXP(+( 7.09562E+09)*TM3+(-4.93578E+07)*TM2                       
     &         +( 7.31751E+04)*TM1+(-1.24136E+02)+( 1.00088E-01)*TP1            
     &         +(-3.06445E-05)*TP2+( 3.57727E-09)*TP3)                          
      RB(50) = RB(50)*PCOR                                                      
      RFH52 =  1.270E+16 * EXP(  -.630*ALOGT -   .38304*RTR )                   
      RFL52 =  2.477E+33 * EXP( -4.760*ALOGT -  2.44026*RTR )                   
      PR = RFL52*XM(52)/RFH52                                                   
      PRLOG = DLOG10(MAX(PR,SMALL))
      F4 =     7.83000E-01
      F5 =     7.40000E+01
      F6 =     2.94100E+03
      F7 =     6.96400E+03
      FC = (1.0D0-F4)*EXP(-T/F5)+F4*EXP(-T/F6)+EXP(-F7/T)
      FCLOG = DLOG10(MAX(FC,SMALL))
      CPRLOG = PRLOG - (0.4D0 + 0.67D0*FCLOG)
      X = CPRLOG/(0.75D0-1.27D0*FCLOG-0.14D0*CPRLOG)
      F = 10.0D0**( FCLOG/(1.0D0+X*X) )
      PCOR = PR*F/(1.0D0+PR)
      RF(52) = RFH52*PCOR                                                       
      RB(52) = EXP(+( 3.59162E+08)*TM3+(-2.35019E+06)*TM2                       
     &         +(-4.69646E+04)*TM1+( 2.98875E+01)+( 4.82138E-03)*TP1            
     &         +(-1.69901E-06)*TP2+( 2.02987E-10)*TP3)                          
      RB(52) = RB(52)*PCOR                                                      
      RFH54 =  1.090E+12 * EXP(   .480*ALOGT +   .26003*RTR )                   
      RFL54 =  1.350E+24 * EXP( -2.570*ALOGT -  1.42515*RTR )                   
      PR = RFL54*XM(54)/RFH54                                                   
      PRLOG = DLOG10(MAX(PR,SMALL))
      F4 =     7.82400E-01
      F5 =     2.71000E+02
      F6 =     2.75500E+03
      F7 =     6.57000E+03
      FC = (1.0D0-F4)*EXP(-T/F5)+F4*EXP(-T/F6)+EXP(-F7/T)
      FCLOG = DLOG10(MAX(FC,SMALL))
      CPRLOG = PRLOG - (0.4D0 + 0.67D0*FCLOG)
      X = CPRLOG/(0.75D0-1.27D0*FCLOG-0.14D0*CPRLOG)
      F = 10.0D0**( FCLOG/(1.0D0+X*X) )
      PCOR = PR*F/(1.0D0+PR)
      RF(54) = RFH54*PCOR                                                       
      RB(54) = EXP(+(-2.96617E+07)*TM3+( 3.57134E+05)*TM2                       
     &         +(-4.54461E+04)*TM1+( 3.60768E+01)+( 1.04094E-04)*TP1            
     &         +(-9.93140E-08)*TP2+( 1.31708E-11)*TP3)                          
      RB(54) = RB(54)*PCOR                                                      
      RFH56 =  5.400E+11 * EXP(   .454*ALOGT -  3.60039*RTR )                   
      RFL56 =  1.270E+32 * EXP( -4.820*ALOGT -  6.53071*RTR )                   
      PR = RFL56*XM(56)/RFH56                                                   
      PRLOG = DLOG10(MAX(PR,SMALL))
      F4 =     7.18700E-01
      F5 =     1.03000E+02
      F6 =     1.29100E+03
      F7 =     4.16000E+03
      FC = (1.0D0-F4)*EXP(-T/F5)+F4*EXP(-T/F6)+EXP(-F7/T)
      FCLOG = DLOG10(MAX(FC,SMALL))
      CPRLOG = PRLOG - (0.4D0 + 0.67D0*FCLOG)
      X = CPRLOG/(0.75D0-1.27D0*FCLOG-0.14D0*CPRLOG)
      F = 10.0D0**( FCLOG/(1.0D0+X*X) )
      PCOR = PR*F/(1.0D0+PR)
      RF(56) = RFH56*PCOR                                                       
      RB(56) = EXP(+( 6.31058E+06)*TM3+(-4.77943E+04)*TM2                       
     &         +(-1.64365E+04)*TM1+( 2.94949E+01)+( 3.14782E-04)*TP1            
     &         +(-1.01093E-07)*TP2+( 1.04005E-11)*TP3)                          
      RB(56) = RB(56)*PCOR                                                      
      RFH57 =  5.400E+11 * EXP(   .454*ALOGT -  2.60028*RTR )                   
      RFL57 =  2.200E+30 * EXP( -4.800*ALOGT -  5.56060*RTR )                   
      PR = RFL57*XM(57)/RFH57                                                   
      PRLOG = DLOG10(MAX(PR,SMALL))
      F4 =     7.58000E-01
      F5 =     9.40000E+01
      F6 =     1.55500E+03
      F7 =     4.20000E+03
      FC = (1.0D0-F4)*EXP(-T/F5)+F4*EXP(-T/F6)+EXP(-F7/T)
      FCLOG = DLOG10(MAX(FC,SMALL))
      CPRLOG = PRLOG - (0.4D0 + 0.67D0*FCLOG)
      X = CPRLOG/(0.75D0-1.27D0*FCLOG-0.14D0*CPRLOG)
      F = 10.0D0**( FCLOG/(1.0D0+X*X) )
      PCOR = PR*F/(1.0D0+PR)
      RF(57) = RFH57*PCOR                                                       
      RB(57) = EXP(+(-1.72463E+07)*TM3+( 2.38003E+05)*TM2                       
     &         +(-1.33724E+04)*TM1+( 3.30241E+01)+(-7.76132E-05)*TP1            
     &         +(-2.92251E-08)*TP2+( 2.95920E-12)*TP3)                          
      RB(57) = RB(57)*PCOR                                                      
      RFH59 =  1.800E+13                                                        
      RFL59 =  3.000E+31 * EXP( -4.800*ALOGT -  3.30036*RTR )                   
      PR = RFL59*XM(59)/RFH59                                                   
      PRLOG = DLOG10(MAX(PR,SMALL))
      F4 =     7.67900E-01
      F5 =     3.38000E+02
      F6 =     1.81200E+03
      F7 =     5.08100E+03
      FC = (1.0D0-F4)*EXP(-T/F5)+F4*EXP(-T/F6)+EXP(-F7/T)
      FCLOG = DLOG10(MAX(FC,SMALL))
      CPRLOG = PRLOG - (0.4D0 + 0.67D0*FCLOG)
      X = CPRLOG/(0.75D0-1.27D0*FCLOG-0.14D0*CPRLOG)
      F = 10.0D0**( FCLOG/(1.0D0+X*X) )
      PCOR = PR*F/(1.0D0+PR)
      RF(59) = RFH59*PCOR                                                       
      RB(59) = EXP(+(-3.60222E+07)*TM3+( 4.34426E+05)*TM2                       
     &         +(-5.00604E+04)*TM1+( 3.60552E+01)+(-2.15830E-04)*TP1            
     &         +(-4.81374E-08)*TP2+( 8.80885E-12)*TP3)                          
      RB(59) = RB(59)*PCOR                                                      
      RFH63 =  5.000E+13                                                        
      RFL63 =  8.600E+28 * EXP( -4.000*ALOGT -  3.02533*RTR )                   
      PR = RFL63*XM(63)/RFH63                                                   
      PRLOG = DLOG10(MAX(PR,SMALL))
      F4 =     8.90200E-01
      F5 =     1.44000E+02
      F6 =     2.83800E+03
      F7 =     4.55690E+04
      FC = (1.0D0-F4)*EXP(-T/F5)+F4*EXP(-T/F6)+EXP(-F7/T)
      FCLOG = DLOG10(MAX(FC,SMALL))
      CPRLOG = PRLOG - (0.4D0 + 0.67D0*FCLOG)
      X = CPRLOG/(0.75D0-1.27D0*FCLOG-0.14D0*CPRLOG)
      F = 10.0D0**( FCLOG/(1.0D0+X*X) )
      PCOR = PR*F/(1.0D0+PR)
      RF(63) = RFH63*PCOR                                                       
      RB(63) = EXP(+( 1.62639E+09)*TM3+(-1.12653E+07)*TM2                       
     &         +(-2.28756E+04)*TM1+(-3.44874E+00)+( 2.34351E-02)*TP1            
     &         +(-7.19893E-06)*TP2+( 8.41318E-10)*TP3)                          
      RB(63) = RB(63)*PCOR                                                      
      RFH70 =  1.000E+17 * EXP( -1.000*ALOGT                )                   
      RFL70 =  3.750E+33 * EXP( -4.800*ALOGT -  1.90021*RTR )                   
      PR = RFL70*XM(70)/RFH70                                                   
      PRLOG = DLOG10(MAX(PR,SMALL))
      F4 =     6.46400E-01
      F5 =     1.32000E+02
      F6 =     1.31500E+03
      F7 =     5.56600E+03
      FC = (1.0D0-F4)*EXP(-T/F5)+F4*EXP(-T/F6)+EXP(-F7/T)
      FCLOG = DLOG10(MAX(FC,SMALL))
      CPRLOG = PRLOG - (0.4D0 + 0.67D0*FCLOG)
      X = CPRLOG/(0.75D0-1.27D0*FCLOG-0.14D0*CPRLOG)
      F = 10.0D0**( FCLOG/(1.0D0+X*X) )
      PCOR = PR*F/(1.0D0+PR)
      RF(70) = RFH70*PCOR                                                       
      RB(70) = EXP(+( 1.73139E+10)*TM3+(-1.08823E+08)*TM2                       
     &         +( 1.93177E+05)*TM1+(-2.62814E+02)+( 1.77135E-01)*TP1            
     &         +(-5.16503E-05)*TP2+( 5.81053E-09)*TP3)                          
      RB(70) = RB(70)*PCOR                                                      
      RFH71 =  5.600E+12 * EXP(              -  2.40026*RTR )                   
      RFL71 =  3.800d+40 * EXP( -7.270*ALOGT -  7.22078*RTR )!mod                   
      PR = RFL71*XM(71)/RFH71                                                   
      PRLOG = DLOG10(MAX(PR,SMALL))
      F4 =     7.50700E-01
      F5 =     9.85000E+01
      F6 =     1.30200E+03
      F7 =     4.16700E+03
      FC = (1.0D0-F4)*EXP(-T/F5)+F4*EXP(-T/F6)+EXP(-F7/T)
      FCLOG = DLOG10(MAX(FC,SMALL))
      CPRLOG = PRLOG - (0.4D0 + 0.67D0*FCLOG)
      X = CPRLOG/(0.75D0-1.27D0*FCLOG-0.14D0*CPRLOG)
      F = 10.0D0**( FCLOG/(1.0D0+X*X) )
      PCOR = PR*F/(1.0D0+PR)
      RF(71) = RFH71*PCOR                                                       
      RB(71) = EXP(+(-2.77264E+07)*TM3+( 3.43893E+05)*TM2                       
     &         +(-1.99235E+04)*TM1+( 3.00214E+01)+(-2.41747E-04)*TP1            
     &         +( 1.77849E-09)*TP2+( 2.94459E-12)*TP3)                          
      RB(71) = RB(71)*PCOR                                                      
      RFH72 =  6.080E+12 * EXP(   .270*ALOGT -   .28003*RTR )                   
      RFL72 =  1.400E+30 * EXP( -3.860*ALOGT -  3.32036*RTR )                   
      PR = RFL72*XM(72)/RFH72                                                   
      PRLOG = DLOG10(MAX(PR,SMALL))
      F4 =     7.82000E-01
      F5 =     2.07500E+02
      F6 =     2.66300E+03
      F7 =     6.09500E+03
      FC = (1.0D0-F4)*EXP(-T/F5)+F4*EXP(-T/F6)+EXP(-F7/T)
      FCLOG = DLOG10(MAX(FC,SMALL))
      CPRLOG = PRLOG - (0.4D0 + 0.67D0*FCLOG)
      X = CPRLOG/(0.75D0-1.27D0*FCLOG-0.14D0*CPRLOG)
      F = 10.0D0**( FCLOG/(1.0D0+X*X) )
      PCOR = PR*F/(1.0D0+PR)
      RF(72) = RFH72*PCOR                                                       
      RB(72) = EXP(+( 8.21559E+09)*TM3+(-5.70716E+07)*TM2                       
     &         +( 9.23731E+04)*TM1+(-1.48549E+02)+( 1.16868E-01)*TP1            
     &         +(-3.56592E-05)*TP2+( 4.15851E-09)*TP3)                          
      RB(72) = RB(72)*PCOR                                                      
      RFH74 =  1.080E+12 * EXP(   .454*ALOGT -  1.82020*RTR )                   
      RFL74 =  1.200d+42 * EXP( -7.620*ALOGT -  6.97076*RTR )!mod                     
      PR = RFL74*XM(74)/RFH74                                                   
      PRLOG = DLOG10(MAX(PR,SMALL))
      F4 =     9.75300E-01
      F5 =     2.10000E+02
      F6 =     9.84000E+02
      F7 =     4.37400E+03
      FC = (1.0D0-F4)*EXP(-T/F5)+F4*EXP(-T/F6)+EXP(-F7/T)
      FCLOG = DLOG10(MAX(FC,SMALL))
      CPRLOG = PRLOG - (0.4D0 + 0.67D0*FCLOG)
      X = CPRLOG/(0.75D0-1.27D0*FCLOG-0.14D0*CPRLOG)
      F = 10.0D0**( FCLOG/(1.0D0+X*X) )
      PCOR = PR*F/(1.0D0+PR)
      RF(74) = RFH74*PCOR                                                       
      RB(74) = EXP(+(-2.05168E+07)*TM3+( 2.25656E+05)*TM2                       
     &         +(-1.98366E+04)*TM1+( 3.10392E+01)+( 1.89973E-04)*TP1            
     &         +(-9.10833E-08)*TP2+( 1.11107E-11)*TP3)                          
      RB(74) = RB(74)*PCOR                                                      
      RFH76 =  5.210E+17 * EXP(  -.990*ALOGT -  1.58017*RTR )                   
      RFL76 =  1.990d+41 * EXP( -7.080*ALOGT -  6.68573*RTR )!mod                     
      PR = RFL76*XM(76)/RFH76                                                   
      PRLOG = DLOG10(MAX(PR,SMALL))
      F4 =     8.42200E-01
      F5 =     1.25000E+02
      F6 =     2.21900E+03
      F7 =     6.88200E+03
      FC = (1.0D0-F4)*EXP(-T/F5)+F4*EXP(-T/F6)+EXP(-F7/T)
      FCLOG = DLOG10(MAX(FC,SMALL))
      CPRLOG = PRLOG - (0.4D0 + 0.67D0*FCLOG)
      X = CPRLOG/(0.75D0-1.27D0*FCLOG-0.14D0*CPRLOG)
      F = 10.0D0**( FCLOG/(1.0D0+X*X) )
      PCOR = PR*F/(1.0D0+PR)
      RF(76) = RFH76*PCOR                                                       
      RB(76) = EXP(+(-1.40102E+07)*TM3+( 1.66802E+05)*TM2                       
     &         +(-5.15163E+04)*TM1+( 4.04437E+01)+(-1.04942E-03)*TP1            
     &         +( 1.07581E-07)*TP2+(-5.45529E-12)*TP3)                          
      RB(76) = RB(76)*PCOR                                                      
      RFH83 =  4.300E+07 * EXP(  1.500*ALOGT - 79.60864*RTR )                   
      RFL83 =  5.070E+27 * EXP( -3.420*ALOGT - 84.35915*RTR )                   
      PR = RFL83*XM(83)/RFH83                                                   
      PRLOG = DLOG10(MAX(PR,SMALL))
      F4 =     9.32000E-01
      F5 =     1.97000E+02
      F6 =     1.54000E+03
      F7 =     1.03000E+04
      FC = (1.0D0-F4)*EXP(-T/F5)+F4*EXP(-T/F6)+EXP(-F7/T)
      FCLOG = DLOG10(MAX(FC,SMALL))
      CPRLOG = PRLOG - (0.4D0 + 0.67D0*FCLOG)
      X = CPRLOG/(0.75D0-1.27D0*FCLOG-0.14D0*CPRLOG)
      F = 10.0D0**( FCLOG/(1.0D0+X*X) )
      PCOR = PR*F/(1.0D0+PR)
      RF(83) = RFH83*PCOR                                                       
      RB(83) = EXP(+(-4.92553E+07)*TM3+( 6.04826E+05)*TM2                       
     &         +(-4.23911E+04)*TM1+( 3.28651E+01)+(-2.19077E-05)*TP1            
     &         +(-5.04236E-08)*TP2+( 9.08484E-12)*TP3)                          
      RB(83) = RB(83)*PCOR                                                      
      RFH85 =  7.400E+13 * EXP(  -.370*ALOGT                )                   
      RFL85 =  2.300E+18 * EXP(  -.900*ALOGT +  1.70018*RTR )                   
      PR = RFL85*XM(85)/RFH85                                                   
      PRLOG = DLOG10(MAX(PR,SMALL))
      F4 =     7.34600E-01
      F5 =     9.40000E+01
      F6 =     1.75600E+03
      F7 =     5.18200E+03
      FC = (1.0D0-F4)*EXP(-T/F5)+F4*EXP(-T/F6)+EXP(-F7/T)
      FCLOG = DLOG10(MAX(FC,SMALL))
      CPRLOG = PRLOG - (0.4D0 + 0.67D0*FCLOG)
      X = CPRLOG/(0.75D0-1.27D0*FCLOG-0.14D0*CPRLOG)
      F = 10.0D0**( FCLOG/(1.0D0+X*X) )
      PCOR = PR*F/(1.0D0+PR)
      RF(85) = RFH85*PCOR                                                       
      RB(85) = EXP(+(-7.24444E+06)*TM3+( 1.46842E+05)*TM2                       
     &         +(-2.62562E+04)*TM1+( 3.62729E+01)+(-1.52775E-03)*TP1            
     &         +( 2.65089E-07)*TP2+(-2.28363E-11)*TP3)                          
      RB(85) = RB(85)*PCOR                                                      
      RFH95 =  6.300E+13                                                        
      RFL95 =  2.700E+38 * EXP( -6.300*ALOGT -  3.10034*RTR )                   
      PR = RFL95*XM(95)/RFH95                                                   
      PRLOG = DLOG10(MAX(PR,SMALL))
      F4 =     2.10500E-01
      F5 =     8.35000E+01
      F6 =     5.39800E+03
      F7 =     8.37000E+03
      FC = (1.0D0-F4)*EXP(-T/F5)+F4*EXP(-T/F6)+EXP(-F7/T)
      FCLOG = DLOG10(MAX(FC,SMALL))
      CPRLOG = PRLOG - (0.4D0 + 0.67D0*FCLOG)
      X = CPRLOG/(0.75D0-1.27D0*FCLOG-0.14D0*CPRLOG)
      F = 10.0D0**( FCLOG/(1.0D0+X*X) )
      PCOR = PR*F/(1.0D0+PR)
      RF(95) = RFH95*PCOR                                                       
      RB(95) = EXP(+(-3.00978E+07)*TM3+( 4.11002E+05)*TM2                       
     &         +(-4.81230E+04)*TM1+( 4.03773E+01)+(-1.36918E-03)*TP1            
     &         +( 2.20682E-07)*TP2+(-1.74721E-11)*TP3)                          
      RB(95) = RB(95)*PCOR                                                      
      RFH131 =  5.000E+13                                                       
      RFL131 =  2.690E+28 * EXP( -3.740*ALOGT -  1.93621*RTR )                  
      PR = RFL131*XM(131)/RFH131                                                
      PRLOG = DLOG10(MAX(PR,SMALL))
      F4 =     5.75700E-01
      F5 =     2.37000E+02
      F6 =     1.65200E+03
      F7 =     5.06900E+03
      FC = (1.0D0-F4)*EXP(-T/F5)+F4*EXP(-T/F6)+EXP(-F7/T)
      FCLOG = DLOG10(MAX(FC,SMALL))
      CPRLOG = PRLOG - (0.4D0 + 0.67D0*FCLOG)
      X = CPRLOG/(0.75D0-1.27D0*FCLOG-0.14D0*CPRLOG)
      F = 10.0D0**( FCLOG/(1.0D0+X*X) )
      PCOR = PR*F/(1.0D0+PR)
      RF(131) = RFH131*PCOR                                                     
      RB(131) = EXP(+( 1.13935E+07)*TM3+(-5.02156E+04)*TM2                      
     &          +(-3.70628E+04)*TM1+( 3.75999E+01)+(-1.46106E-03)*TP1           
     &          +( 2.97031E-07)*TP2+(-2.75444E-11)*TP3)                         
      RB(131) = RB(131)*PCOR                                                    
      RFH140 =  8.100E+11 * EXP(   .500*ALOGT -  4.51049*RTR )                  
      RFL140 =  2.690E+33 * EXP( -5.110*ALOGT -  7.09577*RTR )                  
      PR = RFL140*XM(140)/RFH140                                                
      PRLOG = DLOG10(MAX(PR,SMALL))
      F4 =     5.90700E-01
      F5 =     2.75000E+02
      F6 =     1.22600E+03
      F7 =     5.18500E+03
      FC = (1.0D0-F4)*EXP(-T/F5)+F4*EXP(-T/F6)+EXP(-F7/T)
      FCLOG = DLOG10(MAX(FC,SMALL))
      CPRLOG = PRLOG - (0.4D0 + 0.67D0*FCLOG)
      X = CPRLOG/(0.75D0-1.27D0*FCLOG-0.14D0*CPRLOG)
      F = 10.0D0**( FCLOG/(1.0D0+X*X) )
      PCOR = PR*F/(1.0D0+PR)
      RF(140) = RFH140*PCOR                                                     
      RB(140) = EXP(+(-2.40706E+06)*TM3+( 1.08055E+05)*TM2                      
     &          +(-4.24931E+04)*TM1+( 3.94623E+01)+(-1.30138E-03)*TP1           
     &          +( 2.14502E-07)*TP2+(-1.81351E-11)*TP3)                         
      RB(140) = RB(140)*PCOR                                                    
      RFH147 =  2.000E+13                                                       
      RFL147 =  2.700E+38 * EXP( -6.300*ALOGT -  3.10034*RTR )                  
      PR = RFL147*XM(147)/RFH147                                                
      PRLOG = DLOG10(MAX(PR,SMALL))
      F4 =     1.50700E-01
      F5 =     1.34000E+02
      F6 =     2.38300E+03
      F7 =     7.26500E+03
      FC = (1.0D0-F4)*EXP(-T/F5)+F4*EXP(-T/F6)+EXP(-F7/T)
      FCLOG = DLOG10(MAX(FC,SMALL))
      CPRLOG = PRLOG - (0.4D0 + 0.67D0*FCLOG)
      X = CPRLOG/(0.75D0-1.27D0*FCLOG-0.14D0*CPRLOG)
      F = 10.0D0**( FCLOG/(1.0D0+X*X) )
      PCOR = PR*F/(1.0D0+PR)
      RF(147) = RFH147*PCOR                                                     
      RB(147) = EXP(+(-2.59165E+07)*TM3+( 3.79417E+05)*TM2                      
     &          +(-4.83024E+04)*TM1+( 3.94532E+01)+(-1.80081E-03)*TP1           
     &          +( 3.02933E-07)*TP2+(-2.44290E-11)*TP3)                         
      RB(147) = RB(147)*PCOR                                                    
      RFH158 =  2.120E+16 * EXP(  -.970*ALOGT -   .62007*RTR )                  
      RFL158 =  1.770d+50 * EXP( -9.670*ALOGT -  6.22068*RTR )!mod                    
      PR = RFL158*XM(158)/RFH158                                                
      PRLOG = DLOG10(MAX(PR,SMALL))
      F4 =     5.32500E-01
      F5 =     1.51000E+02
      F6 =     1.03800E+03
      F7 =     4.97000E+03
      FC = (1.0D0-F4)*EXP(-T/F5)+F4*EXP(-T/F6)+EXP(-F7/T)
      FCLOG = DLOG10(MAX(FC,SMALL))
      CPRLOG = PRLOG - (0.4D0 + 0.67D0*FCLOG)
      X = CPRLOG/(0.75D0-1.27D0*FCLOG-0.14D0*CPRLOG)
      F = 10.0D0**( FCLOG/(1.0D0+X*X) )
      PCOR = PR*F/(1.0D0+PR)
      RF(158) = RFH158*PCOR                                                     
      RB(158) = EXP(+(-1.78486E+07)*TM3+( 2.86650E+05)*TM2                      
     &          +(-4.66662E+04)*TM1+( 4.23664E+01)+(-2.56171E-03)*TP1           
     &          +( 4.26863E-07)*TP2+(-3.55820E-11)*TP3)                         
      RB(158) = RB(158)*PCOR                                                    
      RFH174 =  8.000E+12 * EXP(   .440*ALOGT - 88.77963*RTR )                  
      RFL174 =  7.000d+50 * EXP( -9.310*ALOGT - 99.87084*RTR )                  
      PR = RFL174*XM(174)/RFH174                                                
      PRLOG = DLOG10(MAX(PR,SMALL))
      F4 =     7.34500E-01
      F5 =     1.80000E+02
      F6 =     1.03500E+03
      F7 =     5.41700E+03
      FC = (1.0D0-F4)*EXP(-T/F5)+F4*EXP(-T/F6)+EXP(-F7/T)
      FCLOG = DLOG10(MAX(FC,SMALL))
      CPRLOG = PRLOG - (0.4D0 + 0.67D0*FCLOG)
      X = CPRLOG/(0.75D0-1.27D0*FCLOG-0.14D0*CPRLOG)
      F = 10.0D0**( FCLOG/(1.0D0+X*X) )
      PCOR = PR*F/(1.0D0+PR)
      RF(174) = RFH174*PCOR                                                     
      RB(174) = EXP(+( 4.29184E+07)*TM3+(-5.41368E+05)*TM2                      
     &          +(-2.16709E+04)*TM1+( 2.66179E+01)+( 1.36414E-03)*TP1           
     &          +(-1.83619E-07)*TP2+( 1.07724E-11)*TP3)                         
      RB(174) = RB(174)*PCOR                                                    

      RETURN
      END
C
C -------------------------------------------------------------------- C
C
      SUBROUTINE STEADYE( ABV, DEN, RF, RB, XM, ADJ, CONMX, SMALL,
     &  LBOUND, LITER,  XH2, XH, XO, XO2, XOH, XH2O, XHO2, XH2O2, XC,           
     &  XCH, XCH2, XCH2S, XCH3, XCH4, XCO, XCO2, XHCO, XCH2O, XCH2OH,           
     &  XCH3O, XCH3OH, XC2H, XC2H2, XC2H3, XC2H4, XC2H5, XC2H6, XHCCO,          
     &  XCH2CO, XHCCOH, XN2 )                                                   

C*****double precision
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END double precision
C*****single precision
C       IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END single precision

      DIMENSION ABV(*), DEN(*), RF(*), RB(*), XM(*)
      LOGICAL LBOUND(*), LITER

C     STEADY-STATE EXPRESSION FOR HO2       

      ABV(1) = RB(115)*XO2*XH2O2 +RF(175)*XO2*XC2H5 +RB(121)*XH2O2*XHCO         
     &         +RF(157)*XH2O2*XCH3 +RB(118)*XO2*XCH4 +RB(117)*XOH*XCH2O         
     &         +RB(116)*XO2*XH2O2 +RF(32)*XO2*XCH2O +RB(119)*XOH*XCH3O +        
     &         RF(170)*XO2*XCH3O +RF(34)*XH*XO2*XO2 +RF(88)*XOH*XH2O2 +         
     &         RF(5)*XO*XH2O2 +RB(120)*XOH*XCO2 +RF(47)*XH*XH2O2 +              
     &         RF(169)*XO2*XCH2OH +RB(44)*XO*XH2O +RB(4)*XO2*XOH +RF(36)        
     &         *XH*XO2*XN2 +RF(89)*XOH*XH2O2 +RF(33)*XH*XO2*XM(33) +            
     &         RF(168)*XO2*XHCO +RB(45)*XH2*XO2 +RB(87)*XO2*XH2O +RB(46)        
     &         *XOH*XOH +RF(35)*XH*XO2*XH2O                                     
      DEN(1) = RF(115)*XHO2 +RB(175)*XC2H4 +RF(121)*XCH2O +RB(157)*XCH4         
     &         +RF(118)*XCH3 +RF(117)*XCH2 +RF(116)*XHO2 +RB(32)*XHCO +         
     &         RF(119)*XCH3 +RB(170)*XCH2O +RB(34)*XO2 +RB(88)*XH2O +           
     &         RB(5)*XOH +RF(120)*XCO +RB(47)*XH2 +RB(169)*XCH2O +RF(44)        
     &         *XH +RF(4)*XO +RB(36)*XN2 +RB(89)*XH2O +RB(33)*XM(33) +          
     &         RB(168)*XCO +RF(45)*XH +RF(87)*XOH +RF(46)*XH +RB(35)            
     &         *XH2O                                                            
      IF( LITER .OR. .NOT. LBOUND(1) ) THEN                                     
      IF(DEN(1).LT.1.0)DEN(1)=MAX(ADJ*ABV(1),DEN(1),SMALL)                      
      VOLD = XHO2                                                               
      AHO2 = RF(115) +RF(116)                                                   
      BHO2 = DEN(1) - AHO2*XHO2                                                 
      ZHO2 = ABV(1)/BHO2                                                        
      IF( AHO2 .GT. SMALL ) THEN                                                
        B2P4AC = BHO2*BHO2 + 4.0*AHO2*ABV(1)                                    
        XHO2 = 0.5*( SQRT(B2P4AC) - BHO2 )/AHO2                                 
        IF( XHO2 .LT. SMALL ) XHO2 = ZHO2                                       
      ELSE
        XHO2 = ZHO2                                                             
      ENDIF
      DIFF = ABS( (XHO2-VOLD)/MAX(XHO2,VOLD,SMALL))                             
      CONMX = MAX( CONMX, DIFF )
      ENDIF

C     STEADY-STATE EXPRESSION FOR H2O2      

      ABV(2) = RF(115)*XHO2*XHO2 +RF(121)*XHO2*XCH2O +RB(157)*XHO2*XCH4         
     &         +RF(116)*XHO2*XHO2 +RB(88)*XH2O*XHO2 +RB(48)*XOH*XH2O +          
     &         RB(5)*XOH*XHO2 +RB(47)*XH2*XHO2 +RB(89)*XH2O*XHO2 +RF(85)        
     &         *XOH*XOH                                                         
      DEN(2) = RB(115)*XO2 +RB(121)*XHCO +RF(157)*XCH3 +RB(116)*XO2 +           
     &         RF(88)*XOH +RF(48)*XH +RF(5)*XO +RF(47)*XH +RF(89)*XOH +         
     &         RB(85)                                                           
      IF( LITER .OR. .NOT. LBOUND(2) ) THEN                                     
      IF(DEN(2).LT.1.0)DEN(2)=MAX(ADJ*ABV(2),DEN(2),SMALL)                      
      VOLD = XH2O2                                                              
      XH2O2 = ABV(2)/DEN(2)                                                     
      DIFF = ABS( (XH2O2-VOLD)/MAX(XH2O2,VOLD,SMALL))                           
      CONMX = MAX( CONMX, DIFF )
      ENDIF

C     STEADY-STATE EXPRESSION FOR CH3       

      ABV(3) = RB(164)*XCH4*XC2H3 +RB(165)*XCH4*XC2H5 +RF(154)*XCH2S            
     &         *XC2H6 +RB(163)*XCH4*XCH3O +RB(162)*XCH4*XCH2OH +RF(110)         
     &         *XOH*XC2H2 +RB(157)*XHO2*XCH4 +RB(118)*XO2*XCH4 +RB(149)         
     &         *XH*XC2H4 +RF(26)*XO*XC2H5 +RF(139)*XCH2*XCH4 +RB(159)*XH        
     &         *XC2H5 +RB(124)*XH*XC2H2 +RB(129)*XH*XC2H3 +RB(161)*XCH4         
     &         *XHCO +RF(150)*XCH2S*XCH4 +RF(81)*XH*XCH2CO +RB(158)             
     &         *XC2H6 +RB(156)*XOH*XCH2O +RB(160)*XCH4*XCO +RB(119)*XOH         
     &         *XCH3O +RF(25)*XO*XC2H4 +RF(66)*XH*XCH3O +RB(138)*XH             
     &         *XC2H4 +RB(155)*XO*XCH3O +RF(50)*XH*XCH2 +RF(136)*XH2            
     &         *XCH2 +RB(95)*XCH3OH +RF(61)*XH*XCH2OH +RB(52)*XCH4 +            
     &         RF(146)*XH2*XCH2S +RB(96)*XH2O*XCH2 +RF(11)*XO*XCH4 +            
     &         RB(97)*XH2O*XCH2S +RF(53)*XH*XCH4 +RF(98)*XOH*XCH4 +             
     &         RB(10)*XH*XCH2O                                                  
      DEN(3) = RF(164)*XC2H4 +RF(165)*XC2H6 +RB(154)*XC2H5 +RF(163)             
     &         *XCH3OH +RF(162)*XCH3OH +RB(110)*XCO +RF(157)*XH2O2 +            
     &         RF(118)*XHO2 +RF(149)*XCH2S +RB(26)*XCH2O +RB(139)*XCH3 +        
     &         RF(159)*XCH3 +RF(124)*XC +RF(129)*XCH +RF(161)*XCH2O +           
     &         RB(150)*XCH3 +RB(81)*XCO +RF(158)*XCH3 +RF(156)*XO2 +            
     &         RF(160)*XHCO +RF(119)*XHO2 +RB(25)*XHCO +RB(66)*XOH +            
     &         RF(138)*XCH2 +RF(155)*XO2 +RB(50) +RB(136)*XH +RF(95)*XOH        
     &         +RB(61)*XOH +RF(52)*XH +RB(146)*XH +RF(96)*XOH +RB(11)           
     &         *XOH +RF(97)*XOH +RB(53)*XH2 +RB(98)*XH2O +RF(10)*XO             
      IF( LITER .OR. .NOT. LBOUND(3) ) THEN                                     
      IF(DEN(3).LT.1.0)DEN(3)=MAX(ADJ*ABV(3),DEN(3),SMALL)                      
      VOLD = XCH3                                                               
      ACH3 = RB(139) +RF(159) +RB(150) +RF(158)                                 
      BCH3 = DEN(3) - ACH3*XCH3                                                 
      ZCH3 = ABV(3)/BCH3                                                        
      IF( ACH3 .GT. SMALL ) THEN                                                
        B2P4AC = BCH3*BCH3 + 4.0*ACH3*ABV(3)                                    
        XCH3 = 0.5*( SQRT(B2P4AC) - BCH3 )/ACH3                                 
        IF( XCH3 .LT. SMALL ) XCH3 = ZCH3                                       
      ELSE
        XCH3 = ZCH3                                                             
      ENDIF
      DIFF = ABS( (XCH3-VOLD)/MAX(XCH3,VOLD,SMALL))                             
      CONMX = MAX( CONMX, DIFF )
      ENDIF

C     STEADY-STATE EXPRESSION FOR CH2O      

      ABV(4) = RB(121)*XH2O2*XHCO +RF(117)*XHO2*XCH2 +RF(26)*XO*XC2H5 +         
     &         RF(17)*XO*XCH3O +RF(173)*XO2*XC2H3 +RF(103)*XOH*XCH3O +          
     &         RF(54)*XH*XHCO +RB(161)*XCH4*XHCO +RB(32)*XHO2*XHCO +            
     &         RF(156)*XO2*XCH3 +RF(65)*XH*XCH3O +RB(133)*XH*XCH2CO +           
     &         RF(83)*XH2*XCO +RF(170)*XO2*XCH3O +RF(16)*XO*XCH2OH +            
     &         RF(102)*XOH*XCH2OH +RF(60)*XH*XCH2OH +RF(169)*XO2*XCH2OH         
     &         +RB(57)*XCH3O +RF(94)*XOH*XCH2S +RB(56)*XCH2OH +RF(153)          
     &         *XCH2S*XCO2 +RF(92)*XOH*XCH2 +RB(15)*XOH*XHCO +RF(127)           
     &         *XH2O*XCH +RB(58)*XH2*XHCO +RB(101)*XH2O*XHCO +RF(10)*XO         
     &         *XCH3                                                            
      DEN(4) = RF(121)*XHO2 +RB(117)*XOH +RB(26)*XCH3 +RB(17)*XOH +             
     &         RB(173)*XHCO +RB(103)*XH2O +RB(54) +RF(161)*XCH3 +RF(32)         
     &         *XO2 +RB(156)*XOH +RB(65)*XH2 +RF(133)*XCH +RB(83) +             
     &         RB(170)*XHO2 +RB(16)*XOH +RB(102)*XH2O +RB(60)*XH2 +             
     &         RB(169)*XHO2 +RF(57)*XH +RB(94)*XH +RF(56)*XH +RB(153)           
     &         *XCO +RB(92)*XH +RF(15)*XO +RB(127)*XH +RF(58)*XH +              
     &         RF(101)*XOH +RB(10)*XH                                           
      IF( LITER .OR. .NOT. LBOUND(4) ) THEN                                     
      IF(DEN(4).LT.1.0)DEN(4)=MAX(ADJ*ABV(4),DEN(4),SMALL)                      
      VOLD = XCH2O                                                              
      XCH2O = ABV(4)/DEN(4)                                                     
      DIFF = ABS( (XCH2O-VOLD)/MAX(XCH2O,VOLD,SMALL))                           
      CONMX = MAX( CONMX, DIFF )
      ENDIF

C     STEADY-STATE EXPRESSION FOR HCO       

      ABV(5) = RF(121)*XHO2*XCH2O +RF(171)*XO2*XC2H +RF(173)*XO2*XC2H3 +        
     &         RB(54)*XCH2O +RF(161)*XCH3*XCH2O +RF(32)*XO2*XCH2O +             
     &         RB(160)*XCH4*XCO +RF(25)*XO*XC2H4 +RF(9)*XO*XCH2S +              
     &         RF(132)*XCH*XCO2 +RB(13)*XOH*XCO +RB(14)*XH*XCO2 +RF(91)         
     &         *XOH*XCH +RB(168)*XHO2*XCO +RF(135)*XO2*XCH2 +RF(125)*XO2        
     &         *XCH +RF(7)*XO*XCH2 +RB(100)*XH2O*XCO +RF(15)*XO*XCH2O +         
     &         RB(55)*XH2*XCO +RB(167)*XH*XCO*XM(167) +RF(58)*XH*XCH2O +        
     &         RF(101)*XOH*XCH2O +RB(166)*XH*XCO*XH2O                           
      DEN(5) = RB(121)*XH2O2 +RB(171)*XCO +RB(173)*XCH2O +RF(54)*XH +           
     &         RB(161)*XCH4 +RB(32)*XHO2 +RF(160)*XCH3 +RB(25)*XCH3 +           
     &         RB(9)*XH +RB(132)*XCO +RF(13)*XO +RF(14)*XO +RB(91)*XH +         
     &         RF(168)*XO2 +RB(135)*XOH +RB(125)*XO +RB(7)*XH +RF(100)          
     &         *XOH +RB(15)*XOH +RF(55)*XH +RF(167)*XM(167) +RB(58)*XH2         
     &         +RB(101)*XH2O +RF(166)*XH2O                                      
      IF( LITER .OR. .NOT. LBOUND(5) ) THEN                                     
      IF(DEN(5).LT.1.0)DEN(5)=MAX(ADJ*ABV(5),DEN(5),SMALL)                      
      VOLD = XHCO                                                               
      XHCO = ABV(5)/DEN(5)                                                      
      DIFF = ABS( (XHCO-VOLD)/MAX(XHCO,VOLD,SMALL))                             
      CONMX = MAX( CONMX, DIFF )
      ENDIF

C     STEADY-STATE EXPRESSION FOR CH2CO     

      ABV(6) = RF(82)*XHCCOH*XH +RF(107)*XOH*XC2H2 +RB(29)*XOH*XHCCO +          
     &         RB(30)*XCH2*XCO2 +RF(24)*XO*XC2H3 +RB(81)*XCH3*XCO +             
     &         RB(114)*XH2O*XHCCO +RB(80)*XH2*XHCCO +RF(140)*XCH2*XCO +         
     &         RF(133)*XCH*XCH2O                                                
      DEN(6) = RB(82)*XH +RB(107)*XH +RF(29)*XO +RF(30)*XO +RB(24)*XH +         
     &         RF(81)*XH +RF(114)*XOH +RF(80)*XH +RB(140) +RB(133)*XH           
      IF( LITER .OR. .NOT. LBOUND(6) ) THEN                                     
      IF(DEN(6).LT.1.0)DEN(6)=MAX(ADJ*ABV(6),DEN(6),SMALL)                      
      VOLD = XCH2CO                                                             
      XCH2CO = ABV(6)/DEN(6)                                                    
      DIFF = ABS( (XCH2CO-VOLD)/MAX(XCH2CO,VOLD,SMALL))                         
      CONMX = MAX( CONMX, DIFF )
      ENDIF

C     STEADY-STATE EXPRESSION FOR CH3OH     

      ABV(7) = RB(163)*XCH4*XCH3O +RB(162)*XCH4*XCH2OH +RF(63)*XH*XCH3O         
     &         +RF(59)*XH*XCH2OH +RB(19)*XOH*XCH3O +RF(147)*XH2O*XCH2S +        
     &         RB(104)*XH2O*XCH2OH +RB(69)*XH2*XCH3O +RB(18)*XOH*XCH2OH         
     &         +RB(105)*XH2O*XCH3O +RB(68)*XH2*XCH2OH +RF(95)*XOH*XCH3          
      DEN(7) = RF(163)*XCH3 +RF(162)*XCH3 +RB(63) +RB(59) +RF(19)*XO +          
     &         RB(147) +RF(104)*XOH +RF(69)*XH +RF(18)*XO +RF(105)*XOH +        
     &         RF(68)*XH +RB(95)                                                
      IF( LITER .OR. .NOT. LBOUND(7) ) THEN                                     
      IF(DEN(7).LT.1.0)DEN(7)=MAX(ADJ*ABV(7),DEN(7),SMALL)                      
      VOLD = XCH3OH                                                             
      XCH3OH = ABV(7)/DEN(7)                                                    
      DIFF = ABS( (XCH3OH-VOLD)/MAX(XCH3OH,VOLD,SMALL))                         
      CONMX = MAX( CONMX, DIFF )
      ENDIF

C     STEADY-STATE EXPRESSION FOR CH3O      

      ABV(8) = RF(163)*XCH3*XCH3OH +RB(63)*XCH3OH +RB(64)*XCH2OH*XH +           
     &         RB(67)*XH2O*XCH2S +RB(17)*XOH*XCH2O +RB(103)*XH2O*XCH2O +        
     &         RB(65)*XH2*XCH2O +RF(119)*XHO2*XCH3 +RB(66)*XOH*XCH3 +           
     &         RB(170)*XHO2*XCH2O +RF(19)*XO*XCH3OH +RF(155)*XO2*XCH3 +         
     &         RF(69)*XH*XCH3OH +RF(105)*XOH*XCH3OH +RF(57)*XH*XCH2O            
      DEN(8) = RB(163)*XCH4 +RF(63)*XH +RF(64)*XH +RF(67)*XH +RF(17)*XO         
     &         +RF(103)*XOH +RF(65)*XH +RB(119)*XOH +RF(66)*XH +RF(170)         
     &         *XO2 +RB(19)*XOH +RB(155)*XO +RB(69)*XH2 +RB(105)*XH2O +         
     &         RB(57)                                                           
      IF( LITER .OR. .NOT. LBOUND(8) ) THEN                                     
      IF(DEN(8).LT.1.0)DEN(8)=MAX(ADJ*ABV(8),DEN(8),SMALL)                      
      VOLD = XCH3O                                                              
      XCH3O = ABV(8)/DEN(8)                                                     
      DIFF = ABS( (XCH3O-VOLD)/MAX(XCH3O,VOLD,SMALL))                           
      CONMX = MAX( CONMX, DIFF )
      ENDIF

C     STEADY-STATE EXPRESSION FOR CH2OH     

      ABV(9) = RF(162)*XCH3*XCH3OH +RB(59)*XCH3OH +RF(64)*XCH3O*XH +            
     &         RB(16)*XOH*XCH2O +RB(102)*XH2O*XCH2O +RF(104)*XOH*XCH3OH         
     &         +RB(60)*XH2*XCH2O +RF(18)*XO*XCH3OH +RB(62)*XH2O*XCH2S +         
     &         RB(169)*XHO2*XCH2O +RF(68)*XH*XCH3OH +RB(61)*XOH*XCH3 +          
     &         RF(56)*XH*XCH2O                                                  
      DEN(9) = RB(162)*XCH4 +RF(59)*XH +RB(64)*XH +RF(16)*XO +RF(102)           
     &         *XOH +RB(104)*XH2O +RF(60)*XH +RB(18)*XOH +RF(62)*XH +           
     &         RF(169)*XO2 +RB(68)*XH2 +RF(61)*XH +RB(56)                       
      IF( LITER .OR. .NOT. LBOUND(9) ) THEN                                     
      IF(DEN(9).LT.1.0)DEN(9)=MAX(ADJ*ABV(9),DEN(9),SMALL)                      
      VOLD = XCH2OH                                                             
      XCH2OH = ABV(9)/DEN(9)                                                    
      DIFF = ABS( (XCH2OH-VOLD)/MAX(XCH2OH,VOLD,SMALL))                         
      CONMX = MAX( CONMX, DIFF )
      ENDIF

C     STEADY-STATE EXPRESSION FOR C2H6      

      ABV(10) = RB(165)*XCH4*XC2H5 +RB(154)*XCH3*XC2H5 +RF(76)*XH*XC2H5         
     &          +RB(27)*XOH*XC2H5 +RB(78)*XH2*XC2H5 +RB(113)*XH2O*XC2H5         
     &          +RF(158)*XCH3*XCH3                                              
      DEN(10) = RF(165)*XCH3 +RF(154)*XCH2S +RB(76) +RF(27)*XO +RF(78)          
     &          *XH +RF(113)*XOH +RB(158)                                       
      IF( LITER .OR. .NOT. LBOUND(10) ) THEN                                    
      IF(DEN(10).LT.1.0)DEN(10)=MAX(ADJ*ABV(10),DEN(10),SMALL)                  
      VOLD = XC2H6                                                              
      XC2H6 = ABV(10)/DEN(10)                                                   
      DIFF = ABS( (XC2H6-VOLD)/MAX(XC2H6,VOLD,SMALL))                           
      CONMX = MAX( CONMX, DIFF )
      ENDIF

C     STEADY-STATE EXPRESSION FOR CH2S      

      ABV(11) = RB(154)*XCH3*XC2H5 +RB(149)*XH*XC2H4 +RF(67)*XH*XCH3O +         
     &          RB(150)*XCH3*XCH3 +RF(79)*XH*XHCCO +RB(9)*XH*XHCO +RB(8)        
     &          *XH2*XCO +RB(147)*XCH3OH +RF(62)*XH*XCH2OH +RB(151)*XCH2        
     &          *XCO +RB(51)*XH2*XCH +RB(145)*XH2O*XCO +RB(94)*XH*XCH2O         
     &          +RB(152)*XCH2*XCO2 +RB(144)*XH*XOH*XCO +RB(153)*XCO             
     &          *XCH2O +RB(146)*XH*XCH3 +RB(148)*XCH2*XH2O +RB(142)*XCH2        
     &          *XN2 +RF(97)*XOH*XCH3                                           
      DEN(11) = RF(154)*XC2H6 +RF(149)*XCH3 +RB(67)*XH2O +RF(150)*XCH4 +        
     &          RB(79)*XCO +RF(9)*XO +RF(8)*XO +RF(147)*XH2O +RB(62)            
     &          *XH2O +RF(151)*XCO +RF(51)*XH +RF(145)*XO2 +RF(94)*XOH +        
     &          RF(152)*XCO2 +RF(144)*XO2 +RF(153)*XCO2 +RF(146)*XH2 +          
     &          RF(148)*XH2O +RF(142)*XN2 +RB(97)*XH2O                          
      IF( LITER .OR. .NOT. LBOUND(11) ) THEN                                    
      IF(DEN(11).LT.1.0)DEN(11)=MAX(ADJ*ABV(11),DEN(11),SMALL)                  
      VOLD = XCH2S                                                              
      XCH2S = ABV(11)/DEN(11)                                                   
      DIFF = ABS( (XCH2S-VOLD)/MAX(XCH2S,VOLD,SMALL))                           
      CONMX = MAX( CONMX, DIFF )
      ENDIF

C     STEADY-STATE EXPRESSION FOR CH2       

      ABV(12) = RB(141)*XCO*XC2H3 +RB(123)*XH*XC2H +RB(128)*XH*XC2H2 +          
     &          RF(30)*XO*XCH2CO +RB(117)*XOH*XCH2O +RB(139)*XCH3*XCH3 +        
     &          RB(137)*XH2*XC2H2 +RB(140)*XCH2CO +RB(138)*XH*XC2H4 +           
     &          RF(23)*XO*XC2H2 +RB(50)*XCH3 +RF(151)*XCH2S*XCO +RB(136)        
     &          *XH*XCH3 +RF(152)*XCH2S*XCO2 +RB(93)*XH2O*XCH +RB(92)*XH        
     &          *XCH2O +RB(135)*XOH*XHCO +RB(7)*XH*XHCO +RF(96)*XOH*XCH3        
     &          +RF(148)*XCH2S*XH2O +RF(142)*XCH2S*XN2 +RF(126)*XH2*XCH         
      DEN(12) = RF(141)*XHCCO +RF(123)*XC +RF(128)*XCH +RB(30)*XCO2 +           
     &          RF(117)*XHO2 +RF(139)*XCH4 +RF(137)*XCH2 +RF(140)*XCO +         
     &          RF(138)*XCH3 +RB(23)*XCO +RF(50)*XH +RB(151)*XCO +              
     &          RF(136)*XH2 +RB(152)*XCO2 +RF(93)*XOH +RF(92)*XOH +             
     &          RF(135)*XO2 +RF(7)*XO +RB(96)*XH2O +RB(148)*XH2O +              
     &          RB(142)*XN2 +RB(126)*XH                                         
      IF( LITER .OR. .NOT. LBOUND(12) ) THEN                                    
      IF(DEN(12).LT.1.0)DEN(12)=MAX(ADJ*ABV(12),DEN(12),SMALL)                  
      VOLD = XCH2                                                               
      ACH2 = RF(137)                                                            
      BCH2 = DEN(12) - ACH2*XCH2                                                
      ZCH2 = ABV(12)/BCH2                                                       
      IF( ACH2 .GT. SMALL ) THEN                                                
        B2P4AC = BCH2*BCH2 + 4.0*ACH2*ABV(12)                                   
        XCH2 = 0.5*( SQRT(B2P4AC) - BCH2 )/ACH2                                 
        IF( XCH2 .LT. SMALL ) XCH2 = ZCH2                                       
      ELSE
        XCH2 = ZCH2                                                             
      ENDIF
      DIFF = ABS( (XCH2-VOLD)/MAX(XCH2,VOLD,SMALL))                             
      CONMX = MAX( CONMX, DIFF )
      ENDIF

C     STEADY-STATE EXPRESSION FOR CH        

      ABV(13) = RB(134)*XCO*XC2H2 +RF(20)*XO*XC2H +RB(128)*XH*XC2H2 +           
     &          RB(129)*XH*XC2H3 +RB(131)*XHCCO +RB(133)*XH*XCH2CO +            
     &          RB(130)*XH*XC2H4 +RF(51)*XH*XCH2S +RB(6)*XH*XCO +RB(132)        
     &          *XCO*XHCO +RB(91)*XH*XHCO +RF(93)*XOH*XCH2 +RB(125)*XO          
     &          *XHCO +RB(49)*XH2*XC +RB(127)*XH*XCH2O +RB(126)*XH*XCH2         
      DEN(13) = RF(134)*XHCCO +RB(20)*XCO +RF(128)*XCH2 +RF(129)*XCH3 +         
     &          RF(131)*XCO +RF(133)*XCH2O +RF(130)*XCH4 +RB(51)*XH2 +          
     &          RF(6)*XO +RF(132)*XCO2 +RF(91)*XOH +RB(93)*XH2O +RF(125)        
     &          *XO2 +RF(49)*XH +RF(127)*XH2O +RF(126)*XH2                      
      IF( LITER .OR. .NOT. LBOUND(13) ) THEN                                    
      IF(DEN(13).LT.1.0)DEN(13)=MAX(ADJ*ABV(13),DEN(13),SMALL)                  
      VOLD = XCH                                                                
      XCH = ABV(13)/DEN(13)                                                     
      DIFF = ABS( (XCH-VOLD)/MAX(XCH,VOLD,SMALL))                               
      CONMX = MAX( CONMX, DIFF )
      ENDIF

C     STEADY-STATE EXPRESSION FOR C2H4      

      ABV(14) = RB(164)*XCH4*XC2H3 +RF(175)*XO2*XC2H5 +RF(77)*XH*XC2H5 +        
     &          RB(174)*XH2*XC2H2 +RF(72)*XH*XC2H3 +RF(149)*XCH2S*XCH3 +        
     &          RB(74)*XC2H5 +RB(25)*XCH3*XHCO +RB(112)*XH2O*XC2H3 +            
     &          RF(138)*XCH2*XCH3 +RB(75)*XH2*XC2H3 +RF(130)*XCH*XCH4           
      DEN(14) = RF(164)*XCH3 +RB(175)*XHO2 +RB(77)*XH2 +RF(174) +RB(72)         
     &          +RB(149)*XH +RF(74)*XH +RF(25)*XO +RF(112)*XOH +RB(138)         
     &          *XH +RF(75)*XH +RB(130)*XH                                      
      IF( LITER .OR. .NOT. LBOUND(14) ) THEN                                    
      IF(DEN(14).LT.1.0)DEN(14)=MAX(ADJ*ABV(14),DEN(14),SMALL)                  
      VOLD = XC2H4                                                              
      XC2H4 = ABV(14)/DEN(14)                                                   
      DIFF = ABS( (XC2H4-VOLD)/MAX(XC2H4,VOLD,SMALL))                           
      CONMX = MAX( CONMX, DIFF )
      ENDIF

C     STEADY-STATE EXPRESSION FOR C2H3      

      ABV(15) = RF(164)*XCH3*XC2H4 +RF(141)*XCH2*XHCCO +RB(72)*XC2H4 +          
     &          RB(111)*XH2O*XC2H2 +RB(173)*XHCO*XCH2O +RB(24)*XH*XCH2CO        
     &          +RF(129)*XCH*XCH3 +RB(73)*XH2*XC2H2 +RF(112)*XOH*XC2H4 +        
     &          RF(75)*XH*XC2H4 +RF(71)*XH*XC2H2                                
      DEN(15) = RB(164)*XCH4 +RB(141)*XCO +RF(72)*XH +RF(111)*XOH +             
     &          RF(173)*XO2 +RF(24)*XO +RB(129)*XH +RF(73)*XH +RB(112)          
     &          *XH2O +RB(75)*XH2 +RB(71)                                       
      IF( LITER .OR. .NOT. LBOUND(15) ) THEN                                    
      IF(DEN(15).LT.1.0)DEN(15)=MAX(ADJ*ABV(15),DEN(15),SMALL)                  
      VOLD = XC2H3                                                              
      XC2H3 = ABV(15)/DEN(15)                                                   
      DIFF = ABS( (XC2H3-VOLD)/MAX(XC2H3,VOLD,SMALL))                           
      CONMX = MAX( CONMX, DIFF )
      ENDIF

C     STEADY-STATE EXPRESSION FOR C2H2      

      ABV(16) = RF(177)*XHCCO*XHCCO +RF(134)*XCH*XHCCO +RF(70)*XH*XC2H +        
     &          RB(110)*XCH3*XCO +RB(22)*XOH*XC2H +RF(174)*XC2H4 +              
     &          RF(172)*XH2*XC2H +RB(108)*XH*XHCCOH +RB(107)*XH*XCH2CO +        
     &          RF(128)*XCH*XCH2 +RF(137)*XCH2*XCH2 +RB(109)*XH2O*XC2H +        
     &          RF(111)*XOH*XC2H3 +RF(124)*XC*XCH3 +RF(73)*XH*XC2H3 +           
     &          RB(21)*XH*XHCCO +RB(23)*XCH2*XCO +RB(71)*XC2H3                  
      DEN(16) = RB(177)*XCO*XCO +RB(134)*XCO +RB(70) +RF(110)*XOH +             
     &          RF(22)*XO +RB(174)*XH2 +RB(172)*XH +RF(108)*XOH +RF(107)        
     &          *XOH +RB(128)*XH +RB(137)*XH2 +RF(109)*XOH +RB(111)*XH2O        
     &          +RB(124)*XH +RB(73)*XH2 +RF(21)*XO +RF(23)*XO +RF(71)*XH        
      IF( LITER .OR. .NOT. LBOUND(16) ) THEN                                    
      IF(DEN(16).LT.1.0)DEN(16)=MAX(ADJ*ABV(16),DEN(16),SMALL)                  
      VOLD = XC2H2                                                              
      XC2H2 = ABV(16)/DEN(16)                                                   
      DIFF = ABS( (XC2H2-VOLD)/MAX(XC2H2,VOLD,SMALL))                           
      CONMX = MAX( CONMX, DIFF )
      ENDIF

C     STEADY-STATE EXPRESSION FOR HCCO      

      ABV(17) = RB(177)*XCO*XCO*XC2H2 +RB(134)*XCO*XC2H2 +RB(141)*XCO           
     &          *XC2H3 +RF(106)*XOH*XC2H +RF(29)*XO*XCH2CO +RB(176)*XOH         
     &          *XCO*XCO +RF(131)*XCH*XCO +RF(114)*XOH*XCH2CO +RF(80)*XH        
     &          *XCH2CO +RB(28)*XH*XCO*XCO +RF(21)*XO*XC2H2 +RB(79)             
     &          *XCH2S*XCO                                                      
      DEN(17) = RF(177)*XHCCO +RF(134)*XCH +RF(141)*XCH2 +RB(106)*XH +          
     &          RB(29)*XOH +RF(176)*XO2 +RB(131) +RB(114)*XH2O +RB(80)          
     &          *XH2 +RF(28)*XO +RB(21)*XH +RF(79)*XH                           
      IF( LITER .OR. .NOT. LBOUND(17) ) THEN                                    
      IF(DEN(17).LT.1.0)DEN(17)=MAX(ADJ*ABV(17),DEN(17),SMALL)                  
      VOLD = XHCCO                                                              
      AHCCO = RF(177)                                                           
      BHCCO = DEN(17) - AHCCO*XHCCO                                             
      ZHCCO = ABV(17)/BHCCO                                                     
      IF( AHCCO .GT. SMALL ) THEN                                               
        B2P4AC = BHCCO*BHCCO + 4.0*AHCCO*ABV(17)                                
        XHCCO = 0.5*( SQRT(B2P4AC) - BHCCO )/AHCCO                              
        IF( XHCCO .LT. SMALL ) XHCCO = ZHCCO                                    
      ELSE
        XHCCO = ZHCCO                                                           
      ENDIF
      DIFF = ABS( (XHCCO-VOLD)/MAX(XHCCO,VOLD,SMALL))                           
      CONMX = MAX( CONMX, DIFF )
      ENDIF

C     STEADY-STATE EXPRESSION FOR C2H       

      ABV(18) = RB(70)*XC2H2 +RF(22)*XO*XC2H2 +RB(172)*XH*XC2H2 +RB(20)         
     &          *XCH*XCO +RB(106)*XH*XHCCO +RF(123)*XC*XCH2 +RB(171)*XCO        
     &          *XHCO +RF(109)*XOH*XC2H2                                        
      DEN(18) = RF(70)*XH +RB(22)*XOH +RF(172)*XH2 +RF(20)*XO +RF(106)          
     &          *XOH +RB(123)*XH +RF(171)*XO2 +RB(109)*XH2O                     
      IF( LITER .OR. .NOT. LBOUND(18) ) THEN                                    
      IF(DEN(18).LT.1.0)DEN(18)=MAX(ADJ*ABV(18),DEN(18),SMALL)                  
      VOLD = XC2H                                                               
      XC2H = ABV(18)/DEN(18)                                                    
      DIFF = ABS( (XC2H-VOLD)/MAX(XC2H,VOLD,SMALL))                             
      CONMX = MAX( CONMX, DIFF )
      ENDIF

C     STEADY-STATE EXPRESSION FOR HCCOH     

      ABV(19) = RF(108)*XOH*XC2H2 +RB(82)*XCH2CO*XH                             
      DEN(19) = RB(108)*XH +RF(82)*XH                                           
      IF( LITER .OR. .NOT. LBOUND(19) ) THEN                                    
      IF(DEN(19).LT.1.0)DEN(19)=MAX(ADJ*ABV(19),DEN(19),SMALL)                  
      VOLD = XHCCOH                                                             
      XHCCOH = ABV(19)/DEN(19)                                                  
      DIFF = ABS( (XHCCOH-VOLD)/MAX(XHCCOH,VOLD,SMALL))                         
      CONMX = MAX( CONMX, DIFF )
      ENDIF

C     STEADY-STATE EXPRESSION FOR C         

      ABV(20) = RB(123)*XH*XC2H +RB(124)*XH*XC2H2 +RB(90)*XH*XCO +              
     &          RB(122)*XO*XCO +RF(49)*XH*XCH                                   
      DEN(20) = RF(123)*XCH2 +RF(124)*XCH3 +RF(90)*XOH +RF(122)*XO2 +           
     &          RB(49)*XH2                                                      
      IF( LITER .OR. .NOT. LBOUND(20) ) THEN                                    
      IF(DEN(20).LT.1.0)DEN(20)=MAX(ADJ*ABV(20),DEN(20),SMALL)                  
      VOLD = XC                                                                 
      XC = ABV(20)/DEN(20)                                                      
      DIFF = ABS( (XC-VOLD)/MAX(XC,VOLD,SMALL))                                 
      CONMX = MAX( CONMX, DIFF )
      ENDIF

C     STEADY-STATE EXPRESSION FOR C2H5      

      ABV(21) = RF(165)*XCH3*XC2H6 +RF(154)*XCH2S*XC2H6 +RB(175)*XHO2           
     &          *XC2H4 +RB(76)*XC2H6 +RB(77)*XH2*XC2H4 +RF(27)*XO*XC2H6         
     &          +RB(26)*XCH3*XCH2O +RF(78)*XH*XC2H6 +RF(113)*XOH*XC2H6 +        
     &          RF(159)*XCH3*XCH3 +RF(74)*XH*XC2H4                              
      DEN(21) = RB(165)*XCH4 +RB(154)*XCH3 +RF(175)*XO2 +RF(76)*XH +            
     &          RF(77)*XH +RB(27)*XOH +RF(26)*XO +RB(78)*XH2 +RB(113)           
     &          *XH2O +RB(159)*XH +RB(74)                                       
      IF( LITER .OR. .NOT. LBOUND(21) ) THEN                                    
      IF(DEN(21).LT.1.0)DEN(21)=MAX(ADJ*ABV(21),DEN(21),SMALL)                  
      VOLD = XC2H5                                                              
      XC2H5 = ABV(21)/DEN(21)                                                   
      DIFF = ABS( (XC2H5-VOLD)/MAX(XC2H5,VOLD,SMALL))                           
      CONMX = MAX( CONMX, DIFF )
      ENDIF

      RETURN
      END
C
C -------------------------------------------------------------------- C
C
      SUBROUTINE UPVALUE( IOPT, VNEW, BIG, LBOUND, CONMAX, XHO2, XH2O2,         
     &                    XCH3, XCH2O, XHCO, XCH2CO, XCH3OH, XCH3O,             
     &                    XCH2OH, XC2H6, XCH2S, XCH2, XCH, XC2H4, XC2H3,        
     &                    XC2H2, XHCCO, XC2H, XHCCOH, XC, XC2H5 )               

C*****double precision
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END double precision
C*****single precision
C       IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END single precision

      DIMENSION VNEW(*)
      LOGICAL LBOUND(*)

      IF( IOPT .EQ. 1 ) THEN
        XHO2 = VNEW(1)                                                          
        XH2O2 = VNEW(2)                                                         
        XCH3 = VNEW(3)                                                          
        XCH2O = VNEW(4)                                                         
        XHCO = VNEW(5)                                                          
        XCH2CO = VNEW(6)                                                        
        XCH3OH = VNEW(7)                                                        
        XCH3O = VNEW(8)                                                         
        XCH2OH = VNEW(9)                                                        
        XC2H6 = VNEW(10)                                                        
        XCH2S = VNEW(11)                                                        
        XCH2 = VNEW(12)                                                         
        XCH = VNEW(13)                                                          
        XC2H4 = VNEW(14)                                                        
        XC2H3 = VNEW(15)                                                        
        XC2H2 = VNEW(16)                                                        
        XHCCO = VNEW(17)                                                        
        XC2H = VNEW(18)                                                         
        XHCCOH = VNEW(19)                                                       
        XC = VNEW(20)                                                           
        XC2H5 = VNEW(21)                                                        
      ELSEIF( IOPT .EQ. 2 ) THEN
        CALL EXBOUND( XHO2, VNEW(1), LBOUND(1),CONMAX,BIG)                      
        CALL EXBOUND( XH2O2, VNEW(2), LBOUND(2),CONMAX,BIG)                     
        CALL EXBOUND( XCH3, VNEW(3), LBOUND(3),CONMAX,BIG)                      
        CALL EXBOUND( XCH2O, VNEW(4), LBOUND(4),CONMAX,BIG)                     
        CALL EXBOUND( XHCO, VNEW(5), LBOUND(5),CONMAX,BIG)                      
        CALL EXBOUND( XCH2CO, VNEW(6), LBOUND(6),CONMAX,BIG)                    
        CALL EXBOUND( XCH3OH, VNEW(7), LBOUND(7),CONMAX,BIG)                    
        CALL EXBOUND( XCH3O, VNEW(8), LBOUND(8),CONMAX,BIG)                     
        CALL EXBOUND( XCH2OH, VNEW(9), LBOUND(9),CONMAX,BIG)                    
        CALL EXBOUND( XC2H6, VNEW(10), LBOUND(10),CONMAX,BIG)                   
        CALL EXBOUND( XCH2S, VNEW(11), LBOUND(11),CONMAX,BIG)                   
        CALL EXBOUND( XCH2, VNEW(12), LBOUND(12),CONMAX,BIG)                    
        CALL EXBOUND( XCH, VNEW(13), LBOUND(13),CONMAX,BIG)                     
        CALL EXBOUND( XC2H4, VNEW(14), LBOUND(14),CONMAX,BIG)                   
        CALL EXBOUND( XC2H3, VNEW(15), LBOUND(15),CONMAX,BIG)                   
        CALL EXBOUND( XC2H2, VNEW(16), LBOUND(16),CONMAX,BIG)                   
        CALL EXBOUND( XHCCO, VNEW(17), LBOUND(17),CONMAX,BIG)                   
        CALL EXBOUND( XC2H, VNEW(18), LBOUND(18),CONMAX,BIG)                    
        CALL EXBOUND( XHCCOH, VNEW(19), LBOUND(19),CONMAX,BIG)                  
        CALL EXBOUND( XC, VNEW(20), LBOUND(20),CONMAX,BIG)                      
        CALL EXBOUND( XC2H5, VNEW(21), LBOUND(21),CONMAX,BIG)                   
      ELSEIF( IOPT .EQ. 3 ) THEN
        VNEW(1) = XHO2                                                          
        VNEW(2) = XH2O2                                                         
        VNEW(3) = XCH3                                                          
        VNEW(4) = XCH2O                                                         
        VNEW(5) = XHCO                                                          
        VNEW(6) = XCH2CO                                                        
        VNEW(7) = XCH3OH                                                        
        VNEW(8) = XCH3O                                                         
        VNEW(9) = XCH2OH                                                        
        VNEW(10) = XC2H6                                                        
        VNEW(11) = XCH2S                                                        
        VNEW(12) = XCH2                                                         
        VNEW(13) = XCH                                                          
        VNEW(14) = XC2H4                                                        
        VNEW(15) = XC2H3                                                        
        VNEW(16) = XC2H2                                                        
        VNEW(17) = XHCCO                                                        
        VNEW(18) = XC2H                                                         
        VNEW(19) = XHCCOH                                                       
        VNEW(20) = XC                                                           
        VNEW(21) = XC2H5                                                        
      ENDIF

      RETURN
      END
C
C -------------------------------------------------------------------- C
C
      SUBROUTINE EXBOUND( VOLD, B, LBOUND, CONMAX, BIG )

C*****double precision
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END double precision
C*****single precision
C       IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END single precision

      PARAMETER (VBOUND = 100.0)

      LOGICAL LBOUND
      AA = LOG(VOLD) - B
      IF( ABS(AA) .LT. VBOUND ) THEN

        IF( EXP(AA) .LT. BIG ) THEN
          DIFF = ABS( (EXP(AA)-VOLD)/MAX(EXP(AA),VOLD) )
          CONMAX = MAX( CONMAX, DIFF )
          VOLD = EXP(AA)
          LBOUND = .TRUE.
        ELSE
          LBOUND = .FALSE.
        ENDIF

      ELSE
        LBOUND = .FALSE.
      ENDIF

      RETURN
      END
C
C -------------------------------------------------------------------- C
C
      SUBROUTINE NETRATE( W, RF, RB, XM, XH2, XH, XO, XO2, XOH, XH2O,           
     &  XHO2, XH2O2, XC, XCH, XCH2, XCH2S, XCH3, XCH4, XCO, XCO2, XHCO,         
     &  XCH2O, XCH2OH, XCH3O, XCH3OH, XC2H, XC2H2, XC2H3, XC2H4, XC2H5,         
     &  XC2H6, XHCCO, XCH2CO, XHCCOH, XN2 )                                     

C*****double precision
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END double precision
C*****single precision
C       IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END single precision

      DIMENSION W(*), RF(*), RB(*), XM(*)

C   NET PRODUCTION RATES FOR SKELETAL MECHANSIM

      W(1)=RF(1)*XO*XO*XM(1)-RB(1)*XO2*XM(1)                            
      W(2)=RF(2)*XH*XO*XM(2)-RB(2)*XOH*XM(2)                            
      W(3)=RF(3)*XH2*XO-RB(3)*XH*XOH                                    
      W(4)=RF(4)*XO*XHO2-RB(4)*XO2*XOH                                  
      W(5)=RF(5)*XO*XH2O2-RB(5)*XOH*XHO2                                
      W(6)=RF(6)*XO*XCH-RB(6)*XH*XCO                                    
      W(7)=RF(7)*XO*XCH2-RB(7)*XH*XHCO                                  
      W(8)=RF(8)*XO*XCH2S-RB(8)*XH2*XCO                                 
      W(9)=RF(9)*XO*XCH2S-RB(9)*XH*XHCO                                 
      W(10)=RF(10)*XO*XCH3-RB(10)*XH*XCH2O                              
      W(11)=RF(11)*XO*XCH4-RB(11)*XOH*XCH3                              
      W(12)=RF(12)*XO*XCO*XM(12)-RB(12)*XCO2*XM(12)                     
      W(13)=RF(13)*XO*XHCO-RB(13)*XOH*XCO                               
      W(14)=RF(14)*XO*XHCO-RB(14)*XH*XCO2                               
      W(15)=RF(15)*XO*XCH2O-RB(15)*XOH*XHCO                             
      W(16)=RF(16)*XO*XCH2OH-RB(16)*XOH*XCH2O                           
      W(17)=RF(17)*XO*XCH3O-RB(17)*XOH*XCH2O                            
      W(18)=RF(18)*XO*XCH3OH-RB(18)*XOH*XCH2OH                          
      W(19)=RF(19)*XO*XCH3OH-RB(19)*XOH*XCH3O                           
      W(20)=RF(20)*XO*XC2H-RB(20)*XCH*XCO                               
      W(21)=RF(21)*XO*XC2H2-RB(21)*XH*XHCCO                             
      W(22)=RF(22)*XO*XC2H2-RB(22)*XOH*XC2H                             
      W(23)=RF(23)*XO*XC2H2-RB(23)*XCH2*XCO                             
      W(24)=RF(24)*XO*XC2H3-RB(24)*XH*XCH2CO                            
      W(25)=RF(25)*XO*XC2H4-RB(25)*XCH3*XHCO                            
      W(26)=RF(26)*XO*XC2H5-RB(26)*XCH3*XCH2O                           
      W(27)=RF(27)*XO*XC2H6-RB(27)*XOH*XC2H5                            
      W(28)=RF(28)*XO*XHCCO-RB(28)*XH*XCO*XCO                           
      W(29)=RF(29)*XO*XCH2CO-RB(29)*XOH*XHCCO                           
      W(30)=RF(30)*XO*XCH2CO-RB(30)*XCH2*XCO2                           
      W(31)=RF(31)*XO2*XCO-RB(31)*XO*XCO2                               
      W(32)=RF(32)*XO2*XCH2O-RB(32)*XHO2*XHCO                           
      W(33)=RF(33)*XH*XO2*XM(33)-RB(33)*XHO2*XM(33)                     
      W(34)=RF(34)*XH*XO2*XO2-RB(34)*XHO2*XO2                           
      W(35)=RF(35)*XH*XO2*XH2O-RB(35)*XHO2*XH2O                         
      W(36)=RF(36)*XH*XO2*XN2-RB(36)*XHO2*XN2                           
      W(38)=RF(38)*XH*XO2-RB(38)*XO*XOH                                 
      W(39)=RF(39)*XH*XH*XM(39)-RB(39)*XH2*XM(39)                       
      W(40)=RF(40)*XH*XH*XH2-RB(40)*XH2*XH2                             
      W(41)=RF(41)*XH*XH*XH2O-RB(41)*XH2*XH2O                           
      W(42)=RF(42)*XH*XH*XCO2-RB(42)*XH2*XCO2                           
      W(43)=RF(43)*XH*XOH*XM(43)-RB(43)*XH2O*XM(43)                     
      W(44)=RF(44)*XH*XHO2-RB(44)*XO*XH2O                               
      W(45)=RF(45)*XH*XHO2-RB(45)*XH2*XO2                               
      W(46)=RF(46)*XH*XHO2-RB(46)*XOH*XOH                               
      W(47)=RF(47)*XH*XH2O2-RB(47)*XH2*XHO2                             
      W(48)=RF(48)*XH*XH2O2-RB(48)*XOH*XH2O                             
      W(49)=RF(49)*XH*XCH-RB(49)*XH2*XC                                 
      W(50)=RF(50)*XH*XCH2-RB(50)*XCH3                                  
      W(51)=RF(51)*XH*XCH2S-RB(51)*XH2*XCH                              
      W(52)=RF(52)*XH*XCH3-RB(52)*XCH4                                  
      W(53)=RF(53)*XH*XCH4-RB(53)*XH2*XCH3                              
      W(54)=RF(54)*XH*XHCO-RB(54)*XCH2O                                 
      W(55)=RF(55)*XH*XHCO-RB(55)*XH2*XCO                               
      W(56)=RF(56)*XH*XCH2O-RB(56)*XCH2OH                               
      W(57)=RF(57)*XH*XCH2O-RB(57)*XCH3O                                
      W(58)=RF(58)*XH*XCH2O-RB(58)*XH2*XHCO                             
      W(59)=RF(59)*XH*XCH2OH-RB(59)*XCH3OH                              
      W(60)=RF(60)*XH*XCH2OH-RB(60)*XH2*XCH2O                           
      W(61)=RF(61)*XH*XCH2OH-RB(61)*XOH*XCH3                            
      W(62)=RF(62)*XH*XCH2OH-RB(62)*XH2O*XCH2S                          
      W(63)=RF(63)*XH*XCH3O-RB(63)*XCH3OH                               
      W(64)=RF(64)*XCH3O*XH-RB(64)*XCH2OH*XH                            
      W(65)=RF(65)*XH*XCH3O-RB(65)*XH2*XCH2O                            
      W(66)=RF(66)*XH*XCH3O-RB(66)*XOH*XCH3                             
      W(67)=RF(67)*XH*XCH3O-RB(67)*XH2O*XCH2S                           
      W(68)=RF(68)*XH*XCH3OH-RB(68)*XH2*XCH2OH                          
      W(69)=RF(69)*XH*XCH3OH-RB(69)*XH2*XCH3O                           
      W(70)=RF(70)*XH*XC2H-RB(70)*XC2H2                                 
      W(71)=RF(71)*XH*XC2H2-RB(71)*XC2H3                                
      W(72)=RF(72)*XH*XC2H3-RB(72)*XC2H4                                
      W(73)=RF(73)*XH*XC2H3-RB(73)*XH2*XC2H2                            
      W(74)=RF(74)*XH*XC2H4-RB(74)*XC2H5                                
      W(75)=RF(75)*XH*XC2H4-RB(75)*XH2*XC2H3                            
      W(76)=RF(76)*XH*XC2H5-RB(76)*XC2H6                                
      W(77)=RF(77)*XH*XC2H5-RB(77)*XH2*XC2H4                            
      W(78)=RF(78)*XH*XC2H6-RB(78)*XH2*XC2H5                            
      W(79)=RF(79)*XH*XHCCO-RB(79)*XCH2S*XCO                            
      W(80)=RF(80)*XH*XCH2CO-RB(80)*XH2*XHCCO                           
      W(81)=RF(81)*XH*XCH2CO-RB(81)*XCH3*XCO                            
      W(82)=RF(82)*XHCCOH*XH-RB(82)*XCH2CO*XH                           
      W(83)=RF(83)*XH2*XCO-RB(83)*XCH2O                                 
      W(84)=RF(84)*XH2*XOH-RB(84)*XH*XH2O                               
      W(85)=RF(85)*XOH*XOH-RB(85)*XH2O2                                 
      W(86)=RF(86)*XOH*XOH-RB(86)*XO*XH2O                               
      W(87)=RF(87)*XOH*XHO2-RB(87)*XO2*XH2O                             
      W(88)=RF(88)*XOH*XH2O2-RB(88)*XH2O*XHO2                           
      W(89)=RF(89)*XOH*XH2O2-RB(89)*XH2O*XHO2                           
      W(90)=RF(90)*XOH*XC-RB(90)*XH*XCO                                 
      W(91)=RF(91)*XOH*XCH-RB(91)*XH*XHCO                               
      W(92)=RF(92)*XOH*XCH2-RB(92)*XH*XCH2O                             
      W(93)=RF(93)*XOH*XCH2-RB(93)*XH2O*XCH                             
      W(94)=RF(94)*XOH*XCH2S-RB(94)*XH*XCH2O                            
      W(95)=RF(95)*XOH*XCH3-RB(95)*XCH3OH                               
      W(96)=RF(96)*XOH*XCH3-RB(96)*XH2O*XCH2                            
      W(97)=RF(97)*XOH*XCH3-RB(97)*XH2O*XCH2S                           
      W(98)=RF(98)*XOH*XCH4-RB(98)*XH2O*XCH3                            
      W(99)=RF(99)*XOH*XCO-RB(99)*XH*XCO2                               
      W(100)=RF(100)*XOH*XHCO-RB(100)*XH2O*XCO                          
      W(101)=RF(101)*XOH*XCH2O-RB(101)*XH2O*XHCO                        
      W(102)=RF(102)*XOH*XCH2OH-RB(102)*XH2O*XCH2O                      
      W(103)=RF(103)*XOH*XCH3O-RB(103)*XH2O*XCH2O                       
      W(104)=RF(104)*XOH*XCH3OH-RB(104)*XH2O*XCH2OH                     
      W(105)=RF(105)*XOH*XCH3OH-RB(105)*XH2O*XCH3O                      
      W(106)=RF(106)*XOH*XC2H-RB(106)*XH*XHCCO                          
      W(107)=RF(107)*XOH*XC2H2-RB(107)*XH*XCH2CO                        
      W(108)=RF(108)*XOH*XC2H2-RB(108)*XH*XHCCOH                        
      W(109)=RF(109)*XOH*XC2H2-RB(109)*XH2O*XC2H                        
      W(110)=RF(110)*XOH*XC2H2-RB(110)*XCH3*XCO                         
      W(111)=RF(111)*XOH*XC2H3-RB(111)*XH2O*XC2H2                       
      W(112)=RF(112)*XOH*XC2H4-RB(112)*XH2O*XC2H3                       
      W(113)=RF(113)*XOH*XC2H6-RB(113)*XH2O*XC2H5                       
      W(114)=RF(114)*XOH*XCH2CO-RB(114)*XH2O*XHCCO                      
      W(115)=RF(115)*XHO2*XHO2-RB(115)*XO2*XH2O2                        
      W(116)=RF(116)*XHO2*XHO2-RB(116)*XO2*XH2O2                        
      W(117)=RF(117)*XHO2*XCH2-RB(117)*XOH*XCH2O                        
      W(118)=RF(118)*XHO2*XCH3-RB(118)*XO2*XCH4                         
      W(119)=RF(119)*XHO2*XCH3-RB(119)*XOH*XCH3O                        
      W(120)=RF(120)*XHO2*XCO-RB(120)*XOH*XCO2                          
      W(121)=RF(121)*XHO2*XCH2O-RB(121)*XH2O2*XHCO                      
      W(122)=RF(122)*XO2*XC-RB(122)*XO*XCO                              
      W(123)=RF(123)*XC*XCH2-RB(123)*XH*XC2H                            
      W(124)=RF(124)*XC*XCH3-RB(124)*XH*XC2H2                           
      W(125)=RF(125)*XO2*XCH-RB(125)*XO*XHCO                            
      W(126)=RF(126)*XH2*XCH-RB(126)*XH*XCH2                            
      W(127)=RF(127)*XH2O*XCH-RB(127)*XH*XCH2O                          
      W(128)=RF(128)*XCH*XCH2-RB(128)*XH*XC2H2                          
      W(129)=RF(129)*XCH*XCH3-RB(129)*XH*XC2H3                          
      W(130)=RF(130)*XCH*XCH4-RB(130)*XH*XC2H4                          
      W(131)=RF(131)*XCH*XCO-RB(131)*XHCCO                              
      W(132)=RF(132)*XCH*XCO2-RB(132)*XCO*XHCO                          
      W(133)=RF(133)*XCH*XCH2O-RB(133)*XH*XCH2CO                        
      W(134)=RF(134)*XCH*XHCCO-RB(134)*XCO*XC2H2                        
      W(135)=RF(135)*XO2*XCH2-RB(135)*XOH*XHCO                          
      W(136)=RF(136)*XH2*XCH2-RB(136)*XH*XCH3                           
      W(137)=RF(137)*XCH2*XCH2-RB(137)*XH2*XC2H2                        
      W(138)=RF(138)*XCH2*XCH3-RB(138)*XH*XC2H4                         
      W(139)=RF(139)*XCH2*XCH4-RB(139)*XCH3*XCH3                        
      W(140)=RF(140)*XCH2*XCO-RB(140)*XCH2CO                            
      W(141)=RF(141)*XCH2*XHCCO-RB(141)*XCO*XC2H3                       
      W(142)=RF(142)*XCH2S*XN2-RB(142)*XCH2*XN2                         
      W(144)=RF(144)*XO2*XCH2S-RB(144)*XH*XOH*XCO                       
      W(145)=RF(145)*XO2*XCH2S-RB(145)*XH2O*XCO                         
      W(146)=RF(146)*XH2*XCH2S-RB(146)*XH*XCH3                          
      W(147)=RF(147)*XH2O*XCH2S-RB(147)*XCH3OH                          
      W(148)=RF(148)*XCH2S*XH2O-RB(148)*XCH2*XH2O                       
      W(149)=RF(149)*XCH2S*XCH3-RB(149)*XH*XC2H4                        
      W(150)=RF(150)*XCH2S*XCH4-RB(150)*XCH3*XCH3                       
      W(151)=RF(151)*XCH2S*XCO-RB(151)*XCH2*XCO                         
      W(152)=RF(152)*XCH2S*XCO2-RB(152)*XCH2*XCO2                       
      W(153)=RF(153)*XCH2S*XCO2-RB(153)*XCO*XCH2O                       
      W(154)=RF(154)*XCH2S*XC2H6-RB(154)*XCH3*XC2H5                     
      W(155)=RF(155)*XO2*XCH3-RB(155)*XO*XCH3O                          
      W(156)=RF(156)*XO2*XCH3-RB(156)*XOH*XCH2O                         
      W(157)=RF(157)*XH2O2*XCH3-RB(157)*XHO2*XCH4                       
      W(158)=RF(158)*XCH3*XCH3-RB(158)*XC2H6                            
      W(159)=RF(159)*XCH3*XCH3-RB(159)*XH*XC2H5                         
      W(160)=RF(160)*XCH3*XHCO-RB(160)*XCH4*XCO                         
      W(161)=RF(161)*XCH3*XCH2O-RB(161)*XCH4*XHCO                       
      W(162)=RF(162)*XCH3*XCH3OH-RB(162)*XCH4*XCH2OH                    
      W(163)=RF(163)*XCH3*XCH3OH-RB(163)*XCH4*XCH3O                     
      W(164)=RF(164)*XCH3*XC2H4-RB(164)*XCH4*XC2H3                      
      W(165)=RF(165)*XCH3*XC2H6-RB(165)*XCH4*XC2H5                      
      W(166)=RF(166)*XHCO*XH2O-RB(166)*XH*XCO*XH2O                      
      W(167)=RF(167)*XHCO*XM(167)-RB(167)*XH*XCO*XM(167)                
      W(168)=RF(168)*XO2*XHCO-RB(168)*XHO2*XCO                          
      W(169)=RF(169)*XO2*XCH2OH-RB(169)*XHO2*XCH2O                      
      W(170)=RF(170)*XO2*XCH3O-RB(170)*XHO2*XCH2O                       
      W(171)=RF(171)*XO2*XC2H-RB(171)*XCO*XHCO                          
      W(172)=RF(172)*XH2*XC2H-RB(172)*XH*XC2H2                          
      W(173)=RF(173)*XO2*XC2H3-RB(173)*XHCO*XCH2O                       
      W(174)=RF(174)*XC2H4-RB(174)*XH2*XC2H2                            
      W(175)=RF(175)*XO2*XC2H5-RB(175)*XHO2*XC2H4                       
      W(176)=RF(176)*XO2*XHCCO-RB(176)*XOH*XCO*XCO                      
      W(177)=RF(177)*XHCCO*XHCCO-RB(177)*XCO*XCO*XC2H2                  

      RETURN
      END
C -------------------------------------------------------------_C

