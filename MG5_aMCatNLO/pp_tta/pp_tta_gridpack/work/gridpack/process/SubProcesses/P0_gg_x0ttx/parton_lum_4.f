      SUBROUTINE DLUM_4(LUM)
C     ****************************************************            
C         
C     Generated by MadGraph5_aMC@NLO v. 2.5.5, 2017-05-26
C     By the MadGraph5_aMC@NLO Development Team
C     Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
C     RETURNS PARTON LUMINOSITIES FOR MADFKS                          
C        
C     
C     Process: g d~ > x0 t t~ d~ QNP<=1 WEIGHTED<=5 [ all = QCD ]
C     Process: g u~ > x0 t t~ u~ QNP<=1 WEIGHTED<=5 [ all = QCD ]
C     Process: g s~ > x0 t t~ s~ QNP<=1 WEIGHTED<=5 [ all = QCD ]
C     Process: g c~ > x0 t t~ c~ QNP<=1 WEIGHTED<=5 [ all = QCD ]
C     
C     ****************************************************            
C         
      IMPLICIT NONE
C     
C     CONSTANTS                                                       
C         
C     
      INCLUDE 'genps.inc'
      INCLUDE 'nexternal.inc'
      DOUBLE PRECISION       CONV
      PARAMETER (CONV=389379660D0)  !CONV TO PICOBARNS             
C     
C     ARGUMENTS                                                       
C         
C     
      DOUBLE PRECISION PP(0:3,NEXTERNAL), LUM
C     
C     LOCAL VARIABLES                                                 
C         
C     
      INTEGER I, ICROSS,ITYPE,LP
      DOUBLE PRECISION P1(0:3,NEXTERNAL)
      DOUBLE PRECISION G1
      DOUBLE PRECISION CX2,SX2,UX2,DX2
      DOUBLE PRECISION XPQ(-7:7)
C     
C     EXTERNAL FUNCTIONS                                              
C         
C     
      DOUBLE PRECISION ALPHAS2,REWGT,PDG2PDF
C     
C     GLOBAL VARIABLES                                                
C         
C     
      INTEGER              IPROC
      DOUBLE PRECISION PD(0:MAXPROC)
      COMMON /SUBPROC/ PD, IPROC
      INCLUDE 'coupl.inc'
      INCLUDE 'run.inc'
      INTEGER IMIRROR
      COMMON/CMIRROR/IMIRROR
C     
C     DATA                                                            
C         
C     
      DATA G1/1*1D0/
      DATA CX2,SX2,UX2,DX2/4*1D0/
      DATA ICROSS/1/
C     ----------                                                      
C         
C     BEGIN CODE                                                      
C         
C     ----------                                                      
C         
      LUM = 0D0
      IF (IMIRROR.EQ.2) THEN
        IF (ABS(LPP(2)) .GE. 1) THEN
          LP=SIGN(1,LPP(2))
          G1=PDG2PDF(ABS(LPP(2)),0*LP,XBK(2),DSQRT(Q2FACT(2)))
        ENDIF
        IF (ABS(LPP(1)) .GE. 1) THEN
          LP=SIGN(1,LPP(1))
          CX2=PDG2PDF(ABS(LPP(1)),-4*LP,XBK(1),DSQRT(Q2FACT(1)))
          SX2=PDG2PDF(ABS(LPP(1)),-3*LP,XBK(1),DSQRT(Q2FACT(1)))
          UX2=PDG2PDF(ABS(LPP(1)),-2*LP,XBK(1),DSQRT(Q2FACT(1)))
          DX2=PDG2PDF(ABS(LPP(1)),-1*LP,XBK(1),DSQRT(Q2FACT(1)))
        ENDIF
        PD(0) = 0D0
        IPROC = 0
        IPROC=IPROC+1  ! g d~ > x0 t t~ d~
        PD(IPROC) = G1*DX2
        IPROC=IPROC+1  ! g u~ > x0 t t~ u~
        PD(IPROC) = G1*UX2
        IPROC=IPROC+1  ! g s~ > x0 t t~ s~
        PD(IPROC) = G1*SX2
        IPROC=IPROC+1  ! g c~ > x0 t t~ c~
        PD(IPROC) = G1*CX2
      ELSE
        IF (ABS(LPP(1)) .GE. 1) THEN
          LP=SIGN(1,LPP(1))
          G1=PDG2PDF(ABS(LPP(1)),0*LP,XBK(1),DSQRT(Q2FACT(1)))
        ENDIF
        IF (ABS(LPP(2)) .GE. 1) THEN
          LP=SIGN(1,LPP(2))
          CX2=PDG2PDF(ABS(LPP(2)),-4*LP,XBK(2),DSQRT(Q2FACT(2)))
          SX2=PDG2PDF(ABS(LPP(2)),-3*LP,XBK(2),DSQRT(Q2FACT(2)))
          UX2=PDG2PDF(ABS(LPP(2)),-2*LP,XBK(2),DSQRT(Q2FACT(2)))
          DX2=PDG2PDF(ABS(LPP(2)),-1*LP,XBK(2),DSQRT(Q2FACT(2)))
        ENDIF
        PD(0) = 0D0
        IPROC = 0
        IPROC=IPROC+1  ! g d~ > x0 t t~ d~
        PD(IPROC) = G1*DX2
        IPROC=IPROC+1  ! g u~ > x0 t t~ u~
        PD(IPROC) = G1*UX2
        IPROC=IPROC+1  ! g s~ > x0 t t~ s~
        PD(IPROC) = G1*SX2
        IPROC=IPROC+1  ! g c~ > x0 t t~ c~
        PD(IPROC) = G1*CX2
      ENDIF
      DO I=1,IPROC
        IF (NINCOMING.EQ.2) THEN
          LUM = LUM + PD(I) * CONV
        ELSE
          LUM = LUM + PD(I)
        ENDIF
      ENDDO
      RETURN
      END

