      SUBROUTINE MP_COEF_CONSTRUCTION_1(P,NHEL,H,IC)
C     
      USE POLYNOMIAL_CONSTANTS
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=5)
      INTEGER    NCOMB
      PARAMETER (NCOMB=16)

      INTEGER NBORNAMPS
      PARAMETER (NBORNAMPS=2)
      INTEGER    NLOOPS, NLOOPGROUPS, NCTAMPS
      PARAMETER (NLOOPS=35, NLOOPGROUPS=23, NCTAMPS=72)
      INTEGER    NLOOPAMPS
      PARAMETER (NLOOPAMPS=107)
      INTEGER    NWAVEFUNCS,NLOOPWAVEFUNCS
      PARAMETER (NWAVEFUNCS=13,NLOOPWAVEFUNCS=73)
      REAL*16     ZERO
      PARAMETER (ZERO=0.0E0_16)
      COMPLEX*32     IZERO
      PARAMETER (IZERO=CMPLX(0.0E0_16,0.0E0_16,KIND=16))
C     These are constants related to the split orders
      INTEGER    NSO, NSQUAREDSO, NAMPSO
      PARAMETER (NSO=1, NSQUAREDSO=1, NAMPSO=2)
C     
C     ARGUMENTS
C     
      REAL*16 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
      INTEGER H
C     
C     LOCAL VARIABLES
C     
      INTEGER I,J,K
      COMPLEX*32 COEFS(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)
C     
C     GLOBAL VARIABLES
C     
      INCLUDE 'mp_coupl_same_name.inc'

      INTEGER GOODHEL(NCOMB)
      LOGICAL GOODAMP(NSQUAREDSO,NLOOPGROUPS)
      COMMON/FILTERS/GOODAMP,GOODHEL

      INTEGER SQSO_TARGET
      COMMON/SOCHOICE/SQSO_TARGET

      LOGICAL UVCT_REQ_SO_DONE,MP_UVCT_REQ_SO_DONE,CT_REQ_SO_DONE
     $ ,MP_CT_REQ_SO_DONE,LOOP_REQ_SO_DONE,MP_LOOP_REQ_SO_DONE
     $ ,CTCALL_REQ_SO_DONE,FILTER_SO
      COMMON/SO_REQS/UVCT_REQ_SO_DONE,MP_UVCT_REQ_SO_DONE
     $ ,CT_REQ_SO_DONE,MP_CT_REQ_SO_DONE,LOOP_REQ_SO_DONE
     $ ,MP_LOOP_REQ_SO_DONE,CTCALL_REQ_SO_DONE,FILTER_SO

      COMPLEX*32 AMP(NBORNAMPS)
      COMMON/MP_AMPS/AMP
      COMPLEX*32 W(20,NWAVEFUNCS)
      COMMON/MP_W/W

      COMPLEX*32 WL(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE
     $ ,0:NLOOPWAVEFUNCS)
      COMPLEX*32 PL(0:3,0:NLOOPWAVEFUNCS)
      COMMON/MP_WL/WL,PL

      COMPLEX*32 AMPL(3,NCTAMPS)
      COMMON/MP_AMPL/AMPL

C     
C     ----------
C     BEGIN CODE
C     ----------

C     The target squared split order contribution is already reached
C      if true.
      IF (FILTER_SO.AND.MP_LOOP_REQ_SO_DONE) THEN
        GOTO 1001
      ENDIF

C     Coefficient construction for loop diagram with ID 3
      CALL MP_FFV1L2_1(PL(0,0),W(1,6),GC_5,ZERO,ZERO,PL(0,1),COEFS)
      CALL MP_UPDATE_WL_0_1(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,1))
      CALL MP_FFV1L2_1(PL(0,1),W(1,9),GC_5,ZERO,ZERO,PL(0,2),COEFS)
      CALL MP_UPDATE_WL_1_1(WL(1,0,1,1),4,COEFS,4,4,WL(1,0,1,2))
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,2),2,4,1,1,4,73,H)
C     Coefficient construction for loop diagram with ID 4
      CALL MP_FFV1L2_1(PL(0,1),W(1,10),GC_5,ZERO,ZERO,PL(0,3),COEFS)
      CALL MP_UPDATE_WL_1_1(WL(1,0,1,1),4,COEFS,4,4,WL(1,0,1,3))
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,3),2,4,1,1,4,74,H)
C     Coefficient construction for loop diagram with ID 5
      CALL MP_FFV1L3_2(PL(0,0),W(1,1),GC_5,ZERO,ZERO,PL(0,4),COEFS)
      CALL MP_UPDATE_WL_0_1(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,4))
      CALL MP_FFV1L1P0_3(PL(0,4),W(1,2),GC_5,ZERO,ZERO,PL(0,5),COEFS)
      CALL MP_UPDATE_WL_1_0(WL(1,0,1,4),4,COEFS,4,4,WL(1,0,1,5))
      CALL MP_FFV1L3_2(PL(0,5),W(1,5),GC_5,MDL_MT,MDL_WT,PL(0,6),COEFS)
      CALL MP_UPDATE_WL_1_1(WL(1,0,1,5),4,COEFS,4,4,WL(1,0,1,6))
      CALL MP_FFS1_5L1_2(PL(0,6),W(1,3),GC_3010H,GC_3010A,MDL_MT
     $ ,MDL_WT,PL(0,7),COEFS)
      CALL MP_UPDATE_WL_2_1(WL(1,0,1,6),4,COEFS,4,4,WL(1,0,1,7))
      CALL MP_FFV1L1P0_3(PL(0,7),W(1,4),GC_5,ZERO,ZERO,PL(0,8),COEFS)
      CALL MP_UPDATE_WL_3_0(WL(1,0,1,7),4,COEFS,4,4,WL(1,0,1,8))
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,8),3,4,2,1,1,75,H)
C     Coefficient construction for loop diagram with ID 6
      CALL MP_FFV1L3_1(PL(0,5),W(1,4),GC_5,MDL_MT,MDL_WT,PL(0,9),COEFS)
      CALL MP_UPDATE_WL_1_1(WL(1,0,1,5),4,COEFS,4,4,WL(1,0,1,9))
      CALL MP_FFS1_5L2_1(PL(0,9),W(1,3),GC_3010H,GC_3010A,MDL_MT
     $ ,MDL_WT,PL(0,10),COEFS)
      CALL MP_UPDATE_WL_2_1(WL(1,0,1,9),4,COEFS,4,4,WL(1,0,1,10))
      CALL MP_FFV1L2P0_3(PL(0,10),W(1,5),GC_5,ZERO,ZERO,PL(0,11),COEFS)
      CALL MP_UPDATE_WL_3_0(WL(1,0,1,10),4,COEFS,4,4,WL(1,0,1,11))
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,11),3,4,3,1,1,76,H)
C     Coefficient construction for loop diagram with ID 7
      CALL MP_VVV1L2P0_1(PL(0,5),W(1,9),GC_4,ZERO,ZERO,PL(0,12),COEFS)
      CALL MP_UPDATE_WL_1_1(WL(1,0,1,5),4,COEFS,4,4,WL(1,0,1,12))
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,12),2,4,4,1,1,77,H)
C     Coefficient construction for loop diagram with ID 8
      CALL MP_FFV1L1P0_3(PL(0,6),W(1,7),GC_5,ZERO,ZERO,PL(0,13),COEFS)
      CALL MP_UPDATE_WL_2_0(WL(1,0,1,6),4,COEFS,4,4,WL(1,0,1,13))
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,13),2,4,5,1,1,78,H)
C     Coefficient construction for loop diagram with ID 9
      CALL MP_FFV1L3_1(PL(0,5),W(1,7),GC_5,MDL_MT,MDL_WT,PL(0,14)
     $ ,COEFS)
      CALL MP_UPDATE_WL_1_1(WL(1,0,1,5),4,COEFS,4,4,WL(1,0,1,14))
      CALL MP_FFV1L2P0_3(PL(0,14),W(1,5),GC_5,ZERO,ZERO,PL(0,15),COEFS)
      CALL MP_UPDATE_WL_2_0(WL(1,0,1,14),4,COEFS,4,4,WL(1,0,1,15))
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,15),2,4,6,1,1,79,H)
C     Coefficient construction for loop diagram with ID 10
      CALL MP_VVV1L2P0_1(PL(0,5),W(1,10),GC_4,ZERO,ZERO,PL(0,16),COEFS)
      CALL MP_UPDATE_WL_1_1(WL(1,0,1,5),4,COEFS,4,4,WL(1,0,1,16))
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,16),2,4,4,1,1,80,H)
C     Coefficient construction for loop diagram with ID 11
      CALL MP_FFV1L2P0_3(PL(0,9),W(1,8),GC_5,ZERO,ZERO,PL(0,17),COEFS)
      CALL MP_UPDATE_WL_2_0(WL(1,0,1,9),4,COEFS,4,4,WL(1,0,1,17))
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,17),2,4,7,1,1,81,H)
C     Coefficient construction for loop diagram with ID 12
      CALL MP_FFV1L3_2(PL(0,5),W(1,8),GC_5,MDL_MT,MDL_WT,PL(0,18)
     $ ,COEFS)
      CALL MP_UPDATE_WL_1_1(WL(1,0,1,5),4,COEFS,4,4,WL(1,0,1,18))
      CALL MP_FFV1L1P0_3(PL(0,18),W(1,4),GC_5,ZERO,ZERO,PL(0,19),COEFS)
      CALL MP_UPDATE_WL_2_0(WL(1,0,1,18),4,COEFS,4,4,WL(1,0,1,19))
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,19),2,4,8,1,1,82,H)
C     Coefficient construction for loop diagram with ID 13
      CALL MP_FFV1L2P0_3(PL(0,0),W(1,1),GC_5,ZERO,ZERO,PL(0,20),COEFS)
      CALL MP_UPDATE_WL_0_0(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,20))
      CALL MP_FFV1L3_1(PL(0,20),W(1,2),GC_5,ZERO,ZERO,PL(0,21),COEFS)
      CALL MP_UPDATE_WL_0_1(WL(1,0,1,20),4,COEFS,4,4,WL(1,0,1,21))
      CALL MP_FFV1L2_1(PL(0,21),W(1,9),GC_5,ZERO,ZERO,PL(0,22),COEFS)
      CALL MP_UPDATE_WL_1_1(WL(1,0,1,21),4,COEFS,4,4,WL(1,0,1,22))
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,22),2,4,4,1,1,83,H)
C     Coefficient construction for loop diagram with ID 14
      CALL MP_FFV1L2_1(PL(0,21),W(1,10),GC_5,ZERO,ZERO,PL(0,23),COEFS)
      CALL MP_UPDATE_WL_1_1(WL(1,0,1,21),4,COEFS,4,4,WL(1,0,1,23))
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,23),2,4,4,1,1,84,H)
C     Coefficient construction for loop diagram with ID 15
      CALL MP_FFV1L2_1(PL(0,0),W(1,6),GC_5,MDL_MB,ZERO,PL(0,24),COEFS)
      CALL MP_UPDATE_WL_0_1(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,24))
      CALL MP_FFV1L2_1(PL(0,24),W(1,9),GC_5,MDL_MB,ZERO,PL(0,25),COEFS)
      CALL MP_UPDATE_WL_1_1(WL(1,0,1,24),4,COEFS,4,4,WL(1,0,1,25))
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,25),2,4,9,1,1,85,H)
C     Coefficient construction for loop diagram with ID 16
      CALL MP_FFV1L2_1(PL(0,24),W(1,10),GC_5,MDL_MB,ZERO,PL(0,26)
     $ ,COEFS)
      CALL MP_UPDATE_WL_1_1(WL(1,0,1,24),4,COEFS,4,4,WL(1,0,1,26))
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,26),2,4,9,1,1,86,H)
C     Coefficient construction for loop diagram with ID 17
      CALL MP_FFS1_5L2_1(PL(0,0),W(1,3),GC_3009H,GC_3009A,MDL_MB,ZERO
     $ ,PL(0,27),COEFS)
      CALL MP_UPDATE_WL_0_1(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,27))
      CALL MP_FFV1L2_1(PL(0,27),W(1,6),GC_5,MDL_MB,ZERO,PL(0,28),COEFS)
      CALL MP_UPDATE_WL_1_1(WL(1,0,1,27),4,COEFS,4,4,WL(1,0,1,28))
      CALL MP_FFV1L2_1(PL(0,28),W(1,11),GC_5,MDL_MB,ZERO,PL(0,29)
     $ ,COEFS)
      CALL MP_UPDATE_WL_2_1(WL(1,0,1,28),4,COEFS,4,4,WL(1,0,1,29))
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,29),3,4,10,1,1,87,H)
C     Coefficient construction for loop diagram with ID 18
      CALL MP_FFS1_5L1_2(PL(0,0),W(1,3),GC_3009H,GC_3009A,MDL_MB,ZERO
     $ ,PL(0,30),COEFS)
      CALL MP_UPDATE_WL_0_1(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,30))
      CALL MP_FFV1L1_2(PL(0,30),W(1,6),GC_5,MDL_MB,ZERO,PL(0,31),COEFS)
      CALL MP_UPDATE_WL_1_1(WL(1,0,1,30),4,COEFS,4,4,WL(1,0,1,31))
      CALL MP_FFV1L1_2(PL(0,31),W(1,11),GC_5,MDL_MB,ZERO,PL(0,32)
     $ ,COEFS)
      CALL MP_UPDATE_WL_2_1(WL(1,0,1,31),4,COEFS,4,4,WL(1,0,1,32))
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,32),3,4,10,1,1,88,H)
C     Coefficient construction for loop diagram with ID 19
      CALL MP_FFV1L1P0_3(PL(0,0),W(1,7),GC_5,ZERO,ZERO,PL(0,33),COEFS)
      CALL MP_UPDATE_WL_0_0(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,33))
      CALL MP_FFV1L3_2(PL(0,33),W(1,12),GC_5,MDL_MT,MDL_WT,PL(0,34)
     $ ,COEFS)
      CALL MP_UPDATE_WL_0_1(WL(1,0,1,33),4,COEFS,4,4,WL(1,0,1,34))
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,34),1,4,11,1,1,89,H)
C     Coefficient construction for loop diagram with ID 20
      CALL MP_FFV1L2_1(PL(0,0),W(1,6),GC_5,MDL_MT,MDL_WT,PL(0,35)
     $ ,COEFS)
      CALL MP_UPDATE_WL_0_1(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,35))
      CALL MP_FFV1L2_1(PL(0,35),W(1,9),GC_5,MDL_MT,MDL_WT,PL(0,36)
     $ ,COEFS)
      CALL MP_UPDATE_WL_1_1(WL(1,0,1,35),4,COEFS,4,4,WL(1,0,1,36))
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,36),2,4,12,1,1,90,H)
C     Coefficient construction for loop diagram with ID 21
      CALL MP_FFV1L3_2(PL(0,0),W(1,5),GC_5,MDL_MT,MDL_WT,PL(0,37)
     $ ,COEFS)
      CALL MP_UPDATE_WL_0_1(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,37))
      CALL MP_FFV1L1_2(PL(0,37),W(1,6),GC_5,MDL_MT,MDL_WT,PL(0,38)
     $ ,COEFS)
      CALL MP_UPDATE_WL_1_1(WL(1,0,1,37),4,COEFS,4,4,WL(1,0,1,38))
      CALL MP_FFV1L1P0_3(PL(0,38),W(1,7),GC_5,ZERO,ZERO,PL(0,39),COEFS)
      CALL MP_UPDATE_WL_2_0(WL(1,0,1,38),4,COEFS,4,4,WL(1,0,1,39))
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,39),2,4,13,1,1,91,H)
C     Coefficient construction for loop diagram with ID 22
      CALL MP_FFV1L2P0_3(PL(0,0),W(1,5),GC_5,ZERO,ZERO,PL(0,40),COEFS)
      CALL MP_UPDATE_WL_0_0(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,40))
      CALL MP_VVV1L2P0_1(PL(0,40),W(1,6),GC_4,ZERO,ZERO,PL(0,41),COEFS)
      CALL MP_UPDATE_WL_0_1(WL(1,0,1,40),4,COEFS,4,4,WL(1,0,1,41))
      CALL MP_FFV1L3_1(PL(0,41),W(1,7),GC_5,MDL_MT,MDL_WT,PL(0,42)
     $ ,COEFS)
      CALL MP_UPDATE_WL_1_1(WL(1,0,1,41),4,COEFS,4,4,WL(1,0,1,42))
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,42),2,4,14,1,1,92,H)
C     Coefficient construction for loop diagram with ID 23
      CALL MP_FFV1L2P0_3(PL(0,0),W(1,8),GC_5,ZERO,ZERO,PL(0,43),COEFS)
      CALL MP_UPDATE_WL_0_0(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,43))
      CALL MP_FFV1L3_1(PL(0,43),W(1,13),GC_5,MDL_MT,MDL_WT,PL(0,44)
     $ ,COEFS)
      CALL MP_UPDATE_WL_0_1(WL(1,0,1,43),4,COEFS,4,4,WL(1,0,1,44))
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,44),1,4,15,1,1,93,H)
C     Coefficient construction for loop diagram with ID 24
      CALL MP_FFV1L2_1(PL(0,35),W(1,10),GC_5,MDL_MT,MDL_WT,PL(0,45)
     $ ,COEFS)
      CALL MP_UPDATE_WL_1_1(WL(1,0,1,35),4,COEFS,4,4,WL(1,0,1,45))
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,45),2,4,12,1,1,94,H)
C     Coefficient construction for loop diagram with ID 25
      CALL MP_FFV1L3_1(PL(0,0),W(1,4),GC_5,MDL_MT,MDL_WT,PL(0,46)
     $ ,COEFS)
      CALL MP_UPDATE_WL_0_1(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,46))
      CALL MP_FFV1L2_1(PL(0,46),W(1,6),GC_5,MDL_MT,MDL_WT,PL(0,47)
     $ ,COEFS)
      CALL MP_UPDATE_WL_1_1(WL(1,0,1,46),4,COEFS,4,4,WL(1,0,1,47))
      CALL MP_FFV1L2P0_3(PL(0,47),W(1,8),GC_5,ZERO,ZERO,PL(0,48),COEFS)
      CALL MP_UPDATE_WL_2_0(WL(1,0,1,47),4,COEFS,4,4,WL(1,0,1,48))
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,48),2,4,16,1,1,95,H)
C     Coefficient construction for loop diagram with ID 26
      CALL MP_FFV1L1P0_3(PL(0,0),W(1,4),GC_5,ZERO,ZERO,PL(0,49),COEFS)
      CALL MP_UPDATE_WL_0_0(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,49))
      CALL MP_VVV1L2P0_1(PL(0,49),W(1,6),GC_4,ZERO,ZERO,PL(0,50),COEFS)
      CALL MP_UPDATE_WL_0_1(WL(1,0,1,49),4,COEFS,4,4,WL(1,0,1,50))
      CALL MP_FFV1L3_2(PL(0,50),W(1,8),GC_5,MDL_MT,MDL_WT,PL(0,51)
     $ ,COEFS)
      CALL MP_UPDATE_WL_1_1(WL(1,0,1,50),4,COEFS,4,4,WL(1,0,1,51))
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,51),2,4,17,1,1,96,H)
C     Coefficient construction for loop diagram with ID 27
      CALL MP_FFS1_5L2_1(PL(0,0),W(1,3),GC_3010H,GC_3010A,MDL_MT
     $ ,MDL_WT,PL(0,52),COEFS)
      CALL MP_UPDATE_WL_0_1(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,52))
      CALL MP_FFV1L2P0_3(PL(0,52),W(1,5),GC_5,ZERO,ZERO,PL(0,53),COEFS)
      CALL MP_UPDATE_WL_1_0(WL(1,0,1,52),4,COEFS,4,4,WL(1,0,1,53))
      CALL MP_FFV1L3_1(PL(0,53),W(1,13),GC_5,MDL_MT,MDL_WT,PL(0,54)
     $ ,COEFS)
      CALL MP_UPDATE_WL_1_1(WL(1,0,1,53),4,COEFS,4,4,WL(1,0,1,54))
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,54),2,4,18,1,1,97,H)
C     Coefficient construction for loop diagram with ID 28
      CALL MP_FFV1L3_1(PL(0,53),W(1,4),GC_5,MDL_MT,MDL_WT,PL(0,55)
     $ ,COEFS)
      CALL MP_UPDATE_WL_1_1(WL(1,0,1,53),4,COEFS,4,4,WL(1,0,1,55))
      CALL MP_FFV1L2_1(PL(0,55),W(1,6),GC_5,MDL_MT,MDL_WT,PL(0,56)
     $ ,COEFS)
      CALL MP_UPDATE_WL_2_1(WL(1,0,1,55),4,COEFS,4,4,WL(1,0,1,56))
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,56),3,4,19,1,1,98,H)
C     Coefficient construction for loop diagram with ID 29
      CALL MP_FFV1L2_1(PL(0,52),W(1,6),GC_5,MDL_MT,MDL_WT,PL(0,57)
     $ ,COEFS)
      CALL MP_UPDATE_WL_1_1(WL(1,0,1,52),4,COEFS,4,4,WL(1,0,1,57))
      CALL MP_FFV1L2_1(PL(0,57),W(1,11),GC_5,MDL_MT,MDL_WT,PL(0,58)
     $ ,COEFS)
      CALL MP_UPDATE_WL_2_1(WL(1,0,1,57),4,COEFS,4,4,WL(1,0,1,58))
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,58),3,4,20,1,1,99,H)
C     Coefficient construction for loop diagram with ID 30
      CALL MP_FFS1_5L1_2(PL(0,0),W(1,3),GC_3010H,GC_3010A,MDL_MT
     $ ,MDL_WT,PL(0,59),COEFS)
      CALL MP_UPDATE_WL_0_1(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,59))
      CALL MP_FFV1L1_2(PL(0,59),W(1,6),GC_5,MDL_MT,MDL_WT,PL(0,60)
     $ ,COEFS)
      CALL MP_UPDATE_WL_1_1(WL(1,0,1,59),4,COEFS,4,4,WL(1,0,1,60))
      CALL MP_FFV1L1_2(PL(0,60),W(1,11),GC_5,MDL_MT,MDL_WT,PL(0,61)
     $ ,COEFS)
      CALL MP_UPDATE_WL_2_1(WL(1,0,1,60),4,COEFS,4,4,WL(1,0,1,61))
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,61),3,4,20,1,1,100,H)
C     Coefficient construction for loop diagram with ID 31
      CALL MP_FFV1L1P0_3(PL(0,59),W(1,4),GC_5,ZERO,ZERO,PL(0,62),COEFS)
      CALL MP_UPDATE_WL_1_0(WL(1,0,1,59),4,COEFS,4,4,WL(1,0,1,62))
      CALL MP_FFV1L3_2(PL(0,62),W(1,5),GC_5,MDL_MT,MDL_WT,PL(0,63)
     $ ,COEFS)
      CALL MP_UPDATE_WL_1_1(WL(1,0,1,62),4,COEFS,4,4,WL(1,0,1,63))
      CALL MP_FFV1L1_2(PL(0,63),W(1,6),GC_5,MDL_MT,MDL_WT,PL(0,64)
     $ ,COEFS)
      CALL MP_UPDATE_WL_2_1(WL(1,0,1,63),4,COEFS,4,4,WL(1,0,1,64))
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,64),3,4,21,1,1,101,H)
C     Coefficient construction for loop diagram with ID 32
      CALL MP_VVV1L2P0_1(PL(0,62),W(1,6),GC_4,ZERO,ZERO,PL(0,65),COEFS)
      CALL MP_UPDATE_WL_1_1(WL(1,0,1,62),4,COEFS,4,4,WL(1,0,1,65))
      CALL MP_FFV1L3_2(PL(0,65),W(1,5),GC_5,MDL_MT,MDL_WT,PL(0,66)
     $ ,COEFS)
      CALL MP_UPDATE_WL_2_1(WL(1,0,1,65),4,COEFS,4,4,WL(1,0,1,66))
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,66),3,4,22,1,1,102,H)
C     Coefficient construction for loop diagram with ID 33
      CALL MP_FFV1L3_2(PL(0,62),W(1,12),GC_5,MDL_MT,MDL_WT,PL(0,67)
     $ ,COEFS)
      CALL MP_UPDATE_WL_1_1(WL(1,0,1,62),4,COEFS,4,4,WL(1,0,1,67))
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,67),2,4,23,1,1,103,H)
C     Coefficient construction for loop diagram with ID 34
      CALL MP_VVV1L2P0_1(PL(0,0),W(1,6),GC_4,ZERO,ZERO,PL(0,68),COEFS)
      CALL MP_UPDATE_WL_0_1(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,68))
      CALL MP_VVV1L2P0_1(PL(0,68),W(1,9),GC_4,ZERO,ZERO,PL(0,69),COEFS)
      CALL MP_UPDATE_WL_1_1(WL(1,0,1,68),4,COEFS,4,4,WL(1,0,1,69))
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,69),2,4,1,2,1,104,H)
C     Coefficient construction for loop diagram with ID 35
      CALL MP_VVV1L2P0_1(PL(0,68),W(1,10),GC_4,ZERO,ZERO,PL(0,70)
     $ ,COEFS)
      CALL MP_UPDATE_WL_1_1(WL(1,0,1,68),4,COEFS,4,4,WL(1,0,1,70))
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,70),2,4,1,2,1,105,H)
C     Coefficient construction for loop diagram with ID 36
      CALL MP_GHGHGL2_1(PL(0,0),W(1,6),GC_4,ZERO,ZERO,PL(0,71),COEFS)
      CALL MP_UPDATE_WL_0_1(WL(1,0,1,0),1,COEFS,1,1,WL(1,0,1,71))
      CALL MP_GHGHGL2_1(PL(0,71),W(1,9),GC_4,ZERO,ZERO,PL(0,72),COEFS)
      CALL MP_UPDATE_WL_1_1(WL(1,0,1,71),1,COEFS,1,1,WL(1,0,1,72))
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,72),2,1,1,1,1,106,H)
C     Coefficient construction for loop diagram with ID 37
      CALL MP_GHGHGL2_1(PL(0,71),W(1,10),GC_4,ZERO,ZERO,PL(0,73),COEFS)
      CALL MP_UPDATE_WL_1_1(WL(1,0,1,71),1,COEFS,1,1,WL(1,0,1,73))
      CALL MP_CREATE_LOOP_COEFS(WL(1,0,1,73),2,1,1,1,1,107,H)
C     At this point, all loop coefficients needed for (QCD=6), i.e. of
C      split order ID=1, are computed.
      IF(FILTER_SO.AND.SQSO_TARGET.EQ.1) GOTO 4000

      GOTO 1001
 4000 CONTINUE
      MP_LOOP_REQ_SO_DONE=.TRUE.
 1001 CONTINUE
      END

