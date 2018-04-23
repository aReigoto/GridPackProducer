      SUBROUTINE COEF_CONSTRUCTION_1(P,NHEL,H,IC)
C     
C     Modules
C     
      USE POLYNOMIAL_CONSTANTS
C     
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=4)
      INTEGER    NCOMB
      PARAMETER (NCOMB=16)
      INTEGER NBORNAMPS
      PARAMETER (NBORNAMPS=1)
      INTEGER    NLOOPS, NLOOPGROUPS, NCTAMPS
      PARAMETER (NLOOPS=11, NLOOPGROUPS=8, NCTAMPS=29)
      INTEGER    NLOOPAMPS
      PARAMETER (NLOOPAMPS=40)
      INTEGER    NWAVEFUNCS,NLOOPWAVEFUNCS
      PARAMETER (NWAVEFUNCS=6,NLOOPWAVEFUNCS=26)
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
      REAL*16     MP__ZERO
      PARAMETER (MP__ZERO=0.0E0_16)
C     These are constants related to the split orders
      INTEGER    NSO, NSQUAREDSO, NAMPSO
      PARAMETER (NSO=1, NSQUAREDSO=1, NAMPSO=2)
C     
C     ARGUMENTS
C     
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
      INTEGER H
C     
C     LOCAL VARIABLES
C     
      INTEGER I,J,K
      COMPLEX*16 COEFS(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)

      LOGICAL DUMMYFALSE
      DATA DUMMYFALSE/.FALSE./
C     
C     GLOBAL VARIABLES
C     
      INCLUDE 'coupl.inc'
      INCLUDE 'mp_coupl.inc'

      INTEGER HELOFFSET
      INTEGER GOODHEL(NCOMB)
      LOGICAL GOODAMP(NSQUAREDSO,NLOOPGROUPS)
      COMMON/FILTERS/GOODAMP,GOODHEL,HELOFFSET

      LOGICAL CHECKPHASE
      LOGICAL HELDOUBLECHECKED
      COMMON/INIT/CHECKPHASE, HELDOUBLECHECKED

      INTEGER SQSO_TARGET
      COMMON/SOCHOICE/SQSO_TARGET

      LOGICAL UVCT_REQ_SO_DONE,MP_UVCT_REQ_SO_DONE,CT_REQ_SO_DONE
     $ ,MP_CT_REQ_SO_DONE,LOOP_REQ_SO_DONE,MP_LOOP_REQ_SO_DONE
     $ ,CTCALL_REQ_SO_DONE,FILTER_SO
      COMMON/SO_REQS/UVCT_REQ_SO_DONE,MP_UVCT_REQ_SO_DONE
     $ ,CT_REQ_SO_DONE,MP_CT_REQ_SO_DONE,LOOP_REQ_SO_DONE
     $ ,MP_LOOP_REQ_SO_DONE,CTCALL_REQ_SO_DONE,FILTER_SO

      INTEGER I_SO
      COMMON/I_SO/I_SO
      INTEGER I_LIB
      COMMON/I_LIB/I_LIB

      COMPLEX*16 AMP(NBORNAMPS)
      COMMON/AMPS/AMP
      COMPLEX*16 W(20,NWAVEFUNCS)
      COMMON/W/W

      COMPLEX*16 WL(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE
     $ ,0:NLOOPWAVEFUNCS)
      COMPLEX*16 PL(0:3,0:NLOOPWAVEFUNCS)
      COMMON/WL/WL,PL

      COMPLEX*16 AMPL(3,NCTAMPS)
      COMMON/AMPL/AMPL

C     
C     ----------
C     BEGIN CODE
C     ----------

C     The target squared split order contribution is already reached
C      if true.
      IF (FILTER_SO.AND.LOOP_REQ_SO_DONE) THEN
        GOTO 1001
      ENDIF

C     Coefficient construction for loop diagram with ID 2
      CALL FFV1L2_1(PL(0,0),W(1,5),GC_5,ZERO,ZERO,PL(0,1),COEFS)
      CALL UPDATE_WL_0_1(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,1))
      CALL FFV1L2_1(PL(0,1),W(1,6),GC_5,ZERO,ZERO,PL(0,2),COEFS)
      CALL UPDATE_WL_1_1(WL(1,0,1,1),4,COEFS,4,4,WL(1,0,1,2))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,2),2,4,1,1,4,30,H)
C     Coefficient construction for loop diagram with ID 3
      CALL FFV1L3_1(PL(0,0),W(1,1),GC_5,ZERO,ZERO,PL(0,3),COEFS)
      CALL UPDATE_WL_0_1(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,3))
      CALL FFV1L2P0_3(PL(0,3),W(1,2),GC_5,ZERO,ZERO,PL(0,4),COEFS)
      CALL UPDATE_WL_1_0(WL(1,0,1,3),4,COEFS,4,4,WL(1,0,1,4))
      CALL FFV1L3_2(PL(0,4),W(1,4),GC_5,MDL_MT,MDL_WT,PL(0,5),COEFS)
      CALL UPDATE_WL_1_1(WL(1,0,1,4),4,COEFS,4,4,WL(1,0,1,5))
      CALL FFV1L1P0_3(PL(0,5),W(1,3),GC_5,ZERO,ZERO,PL(0,6),COEFS)
      CALL UPDATE_WL_2_0(WL(1,0,1,5),4,COEFS,4,4,WL(1,0,1,6))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,6),2,4,2,1,1,31,H)
C     Coefficient construction for loop diagram with ID 4
      CALL FFV1L3_1(PL(0,4),W(1,3),GC_5,MDL_MT,MDL_WT,PL(0,7),COEFS)
      CALL UPDATE_WL_1_1(WL(1,0,1,4),4,COEFS,4,4,WL(1,0,1,7))
      CALL FFV1L2P0_3(PL(0,7),W(1,4),GC_5,ZERO,ZERO,PL(0,8),COEFS)
      CALL UPDATE_WL_2_0(WL(1,0,1,7),4,COEFS,4,4,WL(1,0,1,8))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,8),2,4,3,1,1,32,H)
C     Coefficient construction for loop diagram with ID 5
      CALL VVV1L2P0_1(PL(0,4),W(1,6),GC_4,ZERO,ZERO,PL(0,9),COEFS)
      CALL UPDATE_WL_1_1(WL(1,0,1,4),4,COEFS,4,4,WL(1,0,1,9))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,9),2,4,4,1,1,33,H)
C     Coefficient construction for loop diagram with ID 6
      CALL FFV1L1P0_3(PL(0,0),W(1,1),GC_5,ZERO,ZERO,PL(0,10),COEFS)
      CALL UPDATE_WL_0_0(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,10))
      CALL FFV1L3_2(PL(0,10),W(1,2),GC_5,ZERO,ZERO,PL(0,11),COEFS)
      CALL UPDATE_WL_0_1(WL(1,0,1,10),4,COEFS,4,4,WL(1,0,1,11))
      CALL FFV1L1_2(PL(0,11),W(1,6),GC_5,ZERO,ZERO,PL(0,12),COEFS)
      CALL UPDATE_WL_1_1(WL(1,0,1,11),4,COEFS,4,4,WL(1,0,1,12))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,12),2,4,4,1,1,34,H)
C     Coefficient construction for loop diagram with ID 7
      CALL FFV1L2_1(PL(0,0),W(1,5),GC_5,MDL_MB,ZERO,PL(0,13),COEFS)
      CALL UPDATE_WL_0_1(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,13))
      CALL FFV1L2_1(PL(0,13),W(1,6),GC_5,MDL_MB,ZERO,PL(0,14),COEFS)
      CALL UPDATE_WL_1_1(WL(1,0,1,13),4,COEFS,4,4,WL(1,0,1,14))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,14),2,4,5,1,1,35,H)
C     Coefficient construction for loop diagram with ID 8
      CALL FFV1L2_1(PL(0,0),W(1,5),GC_5,MDL_MT,MDL_WT,PL(0,15),COEFS)
      CALL UPDATE_WL_0_1(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,15))
      CALL FFV1L2_1(PL(0,15),W(1,6),GC_5,MDL_MT,MDL_WT,PL(0,16),COEFS)
      CALL UPDATE_WL_1_1(WL(1,0,1,15),4,COEFS,4,4,WL(1,0,1,16))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,16),2,4,6,1,1,36,H)
C     Coefficient construction for loop diagram with ID 9
      CALL FFV1L1P0_3(PL(0,0),W(1,3),GC_5,ZERO,ZERO,PL(0,17),COEFS)
      CALL UPDATE_WL_0_0(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,17))
      CALL FFV1L3_2(PL(0,17),W(1,4),GC_5,MDL_MT,MDL_WT,PL(0,18),COEFS)
      CALL UPDATE_WL_0_1(WL(1,0,1,17),4,COEFS,4,4,WL(1,0,1,18))
      CALL FFV1L1_2(PL(0,18),W(1,5),GC_5,MDL_MT,MDL_WT,PL(0,19),COEFS)
      CALL UPDATE_WL_1_1(WL(1,0,1,18),4,COEFS,4,4,WL(1,0,1,19))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,19),2,4,7,1,1,37,H)
C     Coefficient construction for loop diagram with ID 10
      CALL FFV1L3_1(PL(0,0),W(1,3),GC_5,MDL_MT,MDL_WT,PL(0,20),COEFS)
      CALL UPDATE_WL_0_1(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,20))
      CALL FFV1L2P0_3(PL(0,20),W(1,4),GC_5,ZERO,ZERO,PL(0,21),COEFS)
      CALL UPDATE_WL_1_0(WL(1,0,1,20),4,COEFS,4,4,WL(1,0,1,21))
      CALL VVV1L2P0_1(PL(0,21),W(1,5),GC_4,ZERO,ZERO,PL(0,22),COEFS)
      CALL UPDATE_WL_1_1(WL(1,0,1,21),4,COEFS,4,4,WL(1,0,1,22))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,22),2,4,8,1,1,38,H)
C     Coefficient construction for loop diagram with ID 11
      CALL VVV1L2P0_1(PL(0,0),W(1,5),GC_4,ZERO,ZERO,PL(0,23),COEFS)
      CALL UPDATE_WL_0_1(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,23))
      CALL VVV1L2P0_1(PL(0,23),W(1,6),GC_4,ZERO,ZERO,PL(0,24),COEFS)
      CALL UPDATE_WL_1_1(WL(1,0,1,23),4,COEFS,4,4,WL(1,0,1,24))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,24),2,4,1,2,1,39,H)
C     Coefficient construction for loop diagram with ID 12
      CALL GHGHGL2_1(PL(0,0),W(1,5),GC_4,ZERO,ZERO,PL(0,25),COEFS)
      CALL UPDATE_WL_0_1(WL(1,0,1,0),1,COEFS,1,1,WL(1,0,1,25))
      CALL GHGHGL2_1(PL(0,25),W(1,6),GC_4,ZERO,ZERO,PL(0,26),COEFS)
      CALL UPDATE_WL_1_1(WL(1,0,1,25),1,COEFS,1,1,WL(1,0,1,26))
      CALL CREATE_LOOP_COEFS(WL(1,0,1,26),2,1,1,1,1,40,H)
C     At this point, all loop coefficients needed for (QCD=6), i.e. of
C      split order ID=1, are computed.
      IF(FILTER_SO.AND.SQSO_TARGET.EQ.1) GOTO 4000

      GOTO 1001
 4000 CONTINUE
      LOOP_REQ_SO_DONE=.TRUE.
 1001 CONTINUE
      END

