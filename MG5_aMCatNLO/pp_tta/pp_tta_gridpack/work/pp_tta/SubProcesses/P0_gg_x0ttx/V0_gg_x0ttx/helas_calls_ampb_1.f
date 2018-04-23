      SUBROUTINE HELAS_CALLS_AMPB_1(P,NHEL,H,IC)
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
      PARAMETER (NEXTERNAL=5)
      INTEGER    NCOMB
      PARAMETER (NCOMB=16)
      INTEGER NBORNAMPS
      PARAMETER (NBORNAMPS=8)
      INTEGER    NLOOPS, NLOOPGROUPS, NCTAMPS
      PARAMETER (NLOOPS=184, NLOOPGROUPS=89, NCTAMPS=282)
      INTEGER    NLOOPAMPS
      PARAMETER (NLOOPAMPS=466)
      INTEGER    NWAVEFUNCS,NLOOPWAVEFUNCS
      PARAMETER (NWAVEFUNCS=35,NLOOPWAVEFUNCS=325)
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
      IF (FILTER_SO.AND.CT_REQ_SO_DONE) THEN
        GOTO 1001
      ENDIF

      CALL VXXXXX(P(0,1),ZERO,NHEL(1),-1*IC(1),W(1,1))
      CALL VXXXXX(P(0,2),ZERO,NHEL(2),-1*IC(2),W(1,2))
      CALL SXXXXX(P(0,3),+1*IC(3),W(1,3))
      CALL OXXXXX(P(0,4),MDL_MT,NHEL(4),+1*IC(4),W(1,4))
      CALL IXXXXX(P(0,5),MDL_MT,NHEL(5),-1*IC(5),W(1,5))
      CALL VVV1P0_1(W(1,1),W(1,2),GC_4,ZERO,ZERO,W(1,6))
      CALL FFS1_5_1(W(1,4),W(1,3),GC_3010H,GC_3010A,MDL_MT,MDL_WT,W(1
     $ ,7))
C     Amplitude(s) for born diagram with ID 1
      CALL FFV1_0(W(1,5),W(1,7),W(1,6),GC_5,AMP(1))
      CALL FFS1_5_2(W(1,5),W(1,3),GC_3010H,GC_3010A,MDL_MT,MDL_WT,W(1
     $ ,8))
C     Amplitude(s) for born diagram with ID 2
      CALL FFV1_0(W(1,8),W(1,4),W(1,6),GC_5,AMP(2))
      CALL FFV1_1(W(1,4),W(1,1),GC_5,MDL_MT,MDL_WT,W(1,9))
      CALL FFV1_2(W(1,5),W(1,2),GC_5,MDL_MT,MDL_WT,W(1,10))
C     Amplitude(s) for born diagram with ID 3
      CALL FFS1_5_0(W(1,10),W(1,9),W(1,3),GC_3010H,GC_3010A,AMP(3))
C     Amplitude(s) for born diagram with ID 4
      CALL FFV1_0(W(1,8),W(1,9),W(1,2),GC_5,AMP(4))
      CALL FFV1_2(W(1,5),W(1,1),GC_5,MDL_MT,MDL_WT,W(1,11))
      CALL FFV1_1(W(1,4),W(1,2),GC_5,MDL_MT,MDL_WT,W(1,12))
C     Amplitude(s) for born diagram with ID 5
      CALL FFS1_5_0(W(1,11),W(1,12),W(1,3),GC_3010H,GC_3010A,AMP(5))
C     Amplitude(s) for born diagram with ID 6
      CALL FFV1_0(W(1,11),W(1,7),W(1,2),GC_5,AMP(6))
C     Amplitude(s) for born diagram with ID 7
      CALL FFV1_0(W(1,8),W(1,12),W(1,1),GC_5,AMP(7))
C     Amplitude(s) for born diagram with ID 8
      CALL FFV1_0(W(1,10),W(1,7),W(1,1),GC_5,AMP(8))
      CALL FFV1P0_3(W(1,5),W(1,7),GC_5,ZERO,ZERO,W(1,13))
C     Counter-term amplitude(s) for loop diagram number 9
      CALL R2_GG_1_0(W(1,6),W(1,13),R2_GGQ,AMPL(1,1))
      CALL R2_GG_1_0(W(1,6),W(1,13),R2_GGQ,AMPL(1,2))
      CALL R2_GG_1_0(W(1,6),W(1,13),R2_GGQ,AMPL(1,3))
      CALL R2_GG_1_0(W(1,6),W(1,13),R2_GGQ,AMPL(1,4))
      CALL FFV1P0_3(W(1,8),W(1,4),GC_5,ZERO,ZERO,W(1,14))
C     Counter-term amplitude(s) for loop diagram number 10
      CALL R2_GG_1_0(W(1,6),W(1,14),R2_GGQ,AMPL(1,5))
      CALL R2_GG_1_0(W(1,6),W(1,14),R2_GGQ,AMPL(1,6))
      CALL R2_GG_1_0(W(1,6),W(1,14),R2_GGQ,AMPL(1,7))
      CALL R2_GG_1_0(W(1,6),W(1,14),R2_GGQ,AMPL(1,8))
C     Counter-term amplitude(s) for loop diagram number 11
      CALL VVV1_0(W(1,1),W(1,2),W(1,13),UV_3GB_1EPS,AMPL(2,9))
      CALL VVV1_0(W(1,1),W(1,2),W(1,13),UV_3GB_1EPS,AMPL(2,10))
      CALL VVV1_0(W(1,1),W(1,2),W(1,13),UV_3GB_1EPS,AMPL(2,11))
      CALL VVV1_0(W(1,1),W(1,2),W(1,13),UV_3GB_1EPS,AMPL(2,12))
      CALL VVV1_0(W(1,1),W(1,2),W(1,13),UV_3GB,AMPL(1,13))
      CALL VVV1_0(W(1,1),W(1,2),W(1,13),UV_3GB_1EPS,AMPL(2,14))
      CALL VVV1_0(W(1,1),W(1,2),W(1,13),UV_3GT,AMPL(1,15))
      CALL VVV1_0(W(1,1),W(1,2),W(1,13),UV_3GB_1EPS,AMPL(2,16))
      CALL VVV1_0(W(1,1),W(1,2),W(1,13),UV_3GG_1EPS,AMPL(2,17))
      CALL VVV1_0(W(1,1),W(1,2),W(1,13),R2_3GQ,AMPL(1,18))
      CALL VVV1_0(W(1,1),W(1,2),W(1,13),R2_3GQ,AMPL(1,19))
      CALL VVV1_0(W(1,1),W(1,2),W(1,13),R2_3GQ,AMPL(1,20))
      CALL VVV1_0(W(1,1),W(1,2),W(1,13),R2_3GQ,AMPL(1,21))
C     Counter-term amplitude(s) for loop diagram number 12
      CALL VVV1_0(W(1,1),W(1,2),W(1,14),UV_3GB_1EPS,AMPL(2,22))
      CALL VVV1_0(W(1,1),W(1,2),W(1,14),UV_3GB_1EPS,AMPL(2,23))
      CALL VVV1_0(W(1,1),W(1,2),W(1,14),UV_3GB_1EPS,AMPL(2,24))
      CALL VVV1_0(W(1,1),W(1,2),W(1,14),UV_3GB_1EPS,AMPL(2,25))
      CALL VVV1_0(W(1,1),W(1,2),W(1,14),UV_3GB,AMPL(1,26))
      CALL VVV1_0(W(1,1),W(1,2),W(1,14),UV_3GB_1EPS,AMPL(2,27))
      CALL VVV1_0(W(1,1),W(1,2),W(1,14),UV_3GT,AMPL(1,28))
      CALL VVV1_0(W(1,1),W(1,2),W(1,14),UV_3GB_1EPS,AMPL(2,29))
      CALL VVV1_0(W(1,1),W(1,2),W(1,14),UV_3GG_1EPS,AMPL(2,30))
      CALL VVV1_0(W(1,1),W(1,2),W(1,14),R2_3GQ,AMPL(1,31))
      CALL VVV1_0(W(1,1),W(1,2),W(1,14),R2_3GQ,AMPL(1,32))
      CALL VVV1_0(W(1,1),W(1,2),W(1,14),R2_3GQ,AMPL(1,33))
      CALL VVV1_0(W(1,1),W(1,2),W(1,14),R2_3GQ,AMPL(1,34))
C     Counter-term amplitude(s) for loop diagram number 15
      CALL R2_GG_1_R2_GG_3_0(W(1,6),W(1,13),R2_GGQ,R2_GGB,AMPL(1,35))
C     Counter-term amplitude(s) for loop diagram number 16
      CALL R2_GG_1_R2_GG_3_0(W(1,6),W(1,14),R2_GGQ,R2_GGB,AMPL(1,36))
      CALL FFV1P0_3(W(1,5),W(1,4),GC_5,ZERO,ZERO,W(1,15))
C     Counter-term amplitude(s) for loop diagram number 17
      CALL VVS1_0(W(1,6),W(1,15),W(1,3),R2_GGX0_HB,AMPL(1,37))
      CALL FFV1P0_3(W(1,5),W(1,9),GC_5,ZERO,ZERO,W(1,16))
C     Counter-term amplitude(s) for loop diagram number 19
      CALL VVS1_0(W(1,2),W(1,16),W(1,3),R2_GGX0_HB,AMPL(1,38))
      CALL FFV1P0_3(W(1,11),W(1,4),GC_5,ZERO,ZERO,W(1,17))
C     Counter-term amplitude(s) for loop diagram number 21
      CALL VVS1_0(W(1,2),W(1,17),W(1,3),R2_GGX0_HB,AMPL(1,39))
      CALL FFV1P0_3(W(1,5),W(1,12),GC_5,ZERO,ZERO,W(1,18))
C     Counter-term amplitude(s) for loop diagram number 23
      CALL VVS1_0(W(1,1),W(1,18),W(1,3),R2_GGX0_HB,AMPL(1,40))
      CALL FFV1P0_3(W(1,10),W(1,4),GC_5,ZERO,ZERO,W(1,19))
C     Counter-term amplitude(s) for loop diagram number 25
      CALL VVS1_0(W(1,1),W(1,19),W(1,3),R2_GGX0_HB,AMPL(1,41))
C     Counter-term amplitude(s) for loop diagram number 27
      CALL VVV1_0(W(1,1),W(1,2),W(1,13),R2_3GQ,AMPL(1,42))
C     Counter-term amplitude(s) for loop diagram number 28
      CALL VVV1_0(W(1,1),W(1,2),W(1,14),R2_3GQ,AMPL(1,43))
      CALL VVV1P0_1(W(1,2),W(1,15),GC_4,ZERO,ZERO,W(1,20))
C     Counter-term amplitude(s) for loop diagram number 34
      CALL VVS1_0(W(1,1),W(1,20),W(1,3),R2_GGX0_HB,AMPL(1,44))
      CALL VVV1P0_1(W(1,1),W(1,15),GC_4,ZERO,ZERO,W(1,21))
C     Counter-term amplitude(s) for loop diagram number 39
      CALL VVS1_0(W(1,2),W(1,21),W(1,3),R2_GGX0_HB,AMPL(1,45))
      CALL FFV1_2(W(1,5),W(1,6),GC_5,MDL_MT,MDL_WT,W(1,22))
C     Counter-term amplitude(s) for loop diagram number 41
      CALL R2_QQ_1_R2_QQ_2_0(W(1,22),W(1,7),R2_QQQ,R2_QQT,AMPL(1,46))
      CALL R2_QQ_2_0(W(1,22),W(1,7),UV_TMASS,AMPL(1,47))
      CALL R2_QQ_2_0(W(1,22),W(1,7),UV_TMASS_1EPS,AMPL(2,48))
C     Counter-term amplitude(s) for loop diagram number 42
      CALL R2_GG_1_R2_GG_3_0(W(1,6),W(1,13),R2_GGQ,R2_GGT,AMPL(1,49))
C     Counter-term amplitude(s) for loop diagram number 43
      CALL FFV1_0(W(1,5),W(1,7),W(1,6),UV_GQQQ_1EPS,AMPL(2,50))
      CALL FFV1_0(W(1,5),W(1,7),W(1,6),UV_GQQQ_1EPS,AMPL(2,51))
      CALL FFV1_0(W(1,5),W(1,7),W(1,6),UV_GQQQ_1EPS,AMPL(2,52))
      CALL FFV1_0(W(1,5),W(1,7),W(1,6),UV_GQQQ_1EPS,AMPL(2,53))
      CALL FFV1_0(W(1,5),W(1,7),W(1,6),UV_GQQB,AMPL(1,54))
      CALL FFV1_0(W(1,5),W(1,7),W(1,6),UV_GQQQ_1EPS,AMPL(2,55))
      CALL FFV1_0(W(1,5),W(1,7),W(1,6),UV_GQQT,AMPL(1,56))
      CALL FFV1_0(W(1,5),W(1,7),W(1,6),UV_GQQQ_1EPS,AMPL(2,57))
      CALL FFV1_0(W(1,5),W(1,7),W(1,6),UV_GQQG_1EPS,AMPL(2,58))
      CALL FFV1_0(W(1,5),W(1,7),W(1,6),R2_GQQ,AMPL(1,59))
      CALL FFV1_1(W(1,4),W(1,6),GC_5,MDL_MT,MDL_WT,W(1,23))
C     Counter-term amplitude(s) for loop diagram number 45
      CALL R2_QQ_1_R2_QQ_2_0(W(1,8),W(1,23),R2_QQQ,R2_QQT,AMPL(1,60))
      CALL R2_QQ_2_0(W(1,8),W(1,23),UV_TMASS,AMPL(1,61))
      CALL R2_QQ_2_0(W(1,8),W(1,23),UV_TMASS_1EPS,AMPL(2,62))
C     Counter-term amplitude(s) for loop diagram number 46
      CALL R2_GG_1_R2_GG_3_0(W(1,6),W(1,14),R2_GGQ,R2_GGT,AMPL(1,63))
C     Counter-term amplitude(s) for loop diagram number 47
      CALL FFV1_0(W(1,8),W(1,4),W(1,6),UV_GQQQ_1EPS,AMPL(2,64))
      CALL FFV1_0(W(1,8),W(1,4),W(1,6),UV_GQQQ_1EPS,AMPL(2,65))
      CALL FFV1_0(W(1,8),W(1,4),W(1,6),UV_GQQQ_1EPS,AMPL(2,66))
      CALL FFV1_0(W(1,8),W(1,4),W(1,6),UV_GQQQ_1EPS,AMPL(2,67))
      CALL FFV1_0(W(1,8),W(1,4),W(1,6),UV_GQQB,AMPL(1,68))
      CALL FFV1_0(W(1,8),W(1,4),W(1,6),UV_GQQQ_1EPS,AMPL(2,69))
      CALL FFV1_0(W(1,8),W(1,4),W(1,6),UV_GQQT,AMPL(1,70))
      CALL FFV1_0(W(1,8),W(1,4),W(1,6),UV_GQQQ_1EPS,AMPL(2,71))
      CALL FFV1_0(W(1,8),W(1,4),W(1,6),UV_GQQG_1EPS,AMPL(2,72))
      CALL FFV1_0(W(1,8),W(1,4),W(1,6),R2_GQQ,AMPL(1,73))
C     Counter-term amplitude(s) for loop diagram number 49
      CALL FFS1_5_0(W(1,5),W(1,23),W(1,3),UV_X0TT_H,UV_X0TT_A,AMPL(1
     $ ,74))
      CALL FFS1_5_0(W(1,5),W(1,23),W(1,3),UV_X0TT_H_1EPS
     $ ,UV_X0TT_A_1EPS,AMPL(2,75))
      CALL FFS1_5_0(W(1,5),W(1,23),W(1,3),R2_TTX0_H,R2_TTX0_A,AMPL(1
     $ ,76))
C     Counter-term amplitude(s) for loop diagram number 51
      CALL VVS1_0(W(1,6),W(1,15),W(1,3),R2_GGX0_HT,AMPL(1,77))
C     Counter-term amplitude(s) for loop diagram number 55
      CALL FFS1_5_0(W(1,22),W(1,4),W(1,3),UV_X0TT_H,UV_X0TT_A,AMPL(1
     $ ,78))
      CALL FFS1_5_0(W(1,22),W(1,4),W(1,3),UV_X0TT_H_1EPS
     $ ,UV_X0TT_A_1EPS,AMPL(2,79))
      CALL FFS1_5_0(W(1,22),W(1,4),W(1,3),R2_TTX0_H,R2_TTX0_A,AMPL(1
     $ ,80))
      CALL FFS1_5_1(W(1,9),W(1,3),GC_3010H,GC_3010A,MDL_MT,MDL_WT,W(1
     $ ,24))
C     Counter-term amplitude(s) for loop diagram number 56
      CALL R2_QQ_1_R2_QQ_2_0(W(1,10),W(1,24),R2_QQQ,R2_QQT,AMPL(1,81))
      CALL R2_QQ_2_0(W(1,10),W(1,24),UV_TMASS,AMPL(1,82))
      CALL R2_QQ_2_0(W(1,10),W(1,24),UV_TMASS_1EPS,AMPL(2,83))
      CALL FFS1_5_2(W(1,10),W(1,3),GC_3010H,GC_3010A,MDL_MT,MDL_WT,W(1
     $ ,25))
C     Counter-term amplitude(s) for loop diagram number 57
      CALL R2_QQ_1_R2_QQ_2_0(W(1,25),W(1,9),R2_QQQ,R2_QQT,AMPL(1,84))
      CALL R2_QQ_2_0(W(1,25),W(1,9),UV_TMASS,AMPL(1,85))
      CALL R2_QQ_2_0(W(1,25),W(1,9),UV_TMASS_1EPS,AMPL(2,86))
C     Counter-term amplitude(s) for loop diagram number 58
      CALL FFS1_5_0(W(1,10),W(1,9),W(1,3),UV_X0TT_H,UV_X0TT_A,AMPL(1
     $ ,87))
      CALL FFS1_5_0(W(1,10),W(1,9),W(1,3),UV_X0TT_H_1EPS
     $ ,UV_X0TT_A_1EPS,AMPL(2,88))
      CALL FFS1_5_0(W(1,10),W(1,9),W(1,3),R2_TTX0_H,R2_TTX0_A,AMPL(1
     $ ,89))
C     Counter-term amplitude(s) for loop diagram number 59
      CALL FFV1_0(W(1,5),W(1,24),W(1,2),UV_GQQQ_1EPS,AMPL(2,90))
      CALL FFV1_0(W(1,5),W(1,24),W(1,2),UV_GQQQ_1EPS,AMPL(2,91))
      CALL FFV1_0(W(1,5),W(1,24),W(1,2),UV_GQQQ_1EPS,AMPL(2,92))
      CALL FFV1_0(W(1,5),W(1,24),W(1,2),UV_GQQQ_1EPS,AMPL(2,93))
      CALL FFV1_0(W(1,5),W(1,24),W(1,2),UV_GQQB,AMPL(1,94))
      CALL FFV1_0(W(1,5),W(1,24),W(1,2),UV_GQQQ_1EPS,AMPL(2,95))
      CALL FFV1_0(W(1,5),W(1,24),W(1,2),UV_GQQT,AMPL(1,96))
      CALL FFV1_0(W(1,5),W(1,24),W(1,2),UV_GQQQ_1EPS,AMPL(2,97))
      CALL FFV1_0(W(1,5),W(1,24),W(1,2),UV_GQQG_1EPS,AMPL(2,98))
      CALL FFV1_0(W(1,5),W(1,24),W(1,2),R2_GQQ,AMPL(1,99))
C     Counter-term amplitude(s) for loop diagram number 60
      CALL VVS1_0(W(1,2),W(1,16),W(1,3),R2_GGX0_HT,AMPL(1,100))
C     Counter-term amplitude(s) for loop diagram number 62
      CALL FFV1_0(W(1,8),W(1,9),W(1,2),UV_GQQQ_1EPS,AMPL(2,101))
      CALL FFV1_0(W(1,8),W(1,9),W(1,2),UV_GQQQ_1EPS,AMPL(2,102))
      CALL FFV1_0(W(1,8),W(1,9),W(1,2),UV_GQQQ_1EPS,AMPL(2,103))
      CALL FFV1_0(W(1,8),W(1,9),W(1,2),UV_GQQQ_1EPS,AMPL(2,104))
      CALL FFV1_0(W(1,8),W(1,9),W(1,2),UV_GQQB,AMPL(1,105))
      CALL FFV1_0(W(1,8),W(1,9),W(1,2),UV_GQQQ_1EPS,AMPL(2,106))
      CALL FFV1_0(W(1,8),W(1,9),W(1,2),UV_GQQT,AMPL(1,107))
      CALL FFV1_0(W(1,8),W(1,9),W(1,2),UV_GQQQ_1EPS,AMPL(2,108))
      CALL FFV1_0(W(1,8),W(1,9),W(1,2),UV_GQQG_1EPS,AMPL(2,109))
      CALL FFV1_0(W(1,8),W(1,9),W(1,2),R2_GQQ,AMPL(1,110))
      CALL FFV1_1(W(1,9),W(1,2),GC_5,MDL_MT,MDL_WT,W(1,26))
C     Counter-term amplitude(s) for loop diagram number 65
      CALL R2_QQ_1_R2_QQ_2_0(W(1,8),W(1,26),R2_QQQ,R2_QQT,AMPL(1,111))
      CALL R2_QQ_2_0(W(1,8),W(1,26),UV_TMASS,AMPL(1,112))
      CALL R2_QQ_2_0(W(1,8),W(1,26),UV_TMASS_1EPS,AMPL(2,113))
      CALL FFV1_2(W(1,8),W(1,2),GC_5,MDL_MT,MDL_WT,W(1,27))
C     Counter-term amplitude(s) for loop diagram number 66
      CALL R2_QQ_1_R2_QQ_2_0(W(1,27),W(1,9),R2_QQQ,R2_QQT,AMPL(1,114))
      CALL R2_QQ_2_0(W(1,27),W(1,9),UV_TMASS,AMPL(1,115))
      CALL R2_QQ_2_0(W(1,27),W(1,9),UV_TMASS_1EPS,AMPL(2,116))
C     Counter-term amplitude(s) for loop diagram number 68
      CALL FFS1_5_0(W(1,5),W(1,26),W(1,3),UV_X0TT_H,UV_X0TT_A,AMPL(1
     $ ,117))
      CALL FFS1_5_0(W(1,5),W(1,26),W(1,3),UV_X0TT_H_1EPS
     $ ,UV_X0TT_A_1EPS,AMPL(2,118))
      CALL FFS1_5_0(W(1,5),W(1,26),W(1,3),R2_TTX0_H,R2_TTX0_A,AMPL(1
     $ ,119))
      CALL FFS1_5_2(W(1,11),W(1,3),GC_3010H,GC_3010A,MDL_MT,MDL_WT,W(1
     $ ,28))
C     Counter-term amplitude(s) for loop diagram number 71
      CALL R2_QQ_1_R2_QQ_2_0(W(1,28),W(1,12),R2_QQQ,R2_QQT,AMPL(1,120))
      CALL R2_QQ_2_0(W(1,28),W(1,12),UV_TMASS,AMPL(1,121))
      CALL R2_QQ_2_0(W(1,28),W(1,12),UV_TMASS_1EPS,AMPL(2,122))
      CALL FFS1_5_1(W(1,12),W(1,3),GC_3010H,GC_3010A,MDL_MT,MDL_WT,W(1
     $ ,29))
C     Counter-term amplitude(s) for loop diagram number 72
      CALL R2_QQ_1_R2_QQ_2_0(W(1,11),W(1,29),R2_QQQ,R2_QQT,AMPL(1,123))
      CALL R2_QQ_2_0(W(1,11),W(1,29),UV_TMASS,AMPL(1,124))
      CALL R2_QQ_2_0(W(1,11),W(1,29),UV_TMASS_1EPS,AMPL(2,125))
C     Counter-term amplitude(s) for loop diagram number 73
      CALL FFS1_5_0(W(1,11),W(1,12),W(1,3),UV_X0TT_H,UV_X0TT_A,AMPL(1
     $ ,126))
      CALL FFS1_5_0(W(1,11),W(1,12),W(1,3),UV_X0TT_H_1EPS
     $ ,UV_X0TT_A_1EPS,AMPL(2,127))
      CALL FFS1_5_0(W(1,11),W(1,12),W(1,3),R2_TTX0_H,R2_TTX0_A,AMPL(1
     $ ,128))
C     Counter-term amplitude(s) for loop diagram number 74
      CALL VVS1_0(W(1,2),W(1,17),W(1,3),R2_GGX0_HT,AMPL(1,129))
C     Counter-term amplitude(s) for loop diagram number 75
      CALL FFV1_0(W(1,11),W(1,7),W(1,2),UV_GQQQ_1EPS,AMPL(2,130))
      CALL FFV1_0(W(1,11),W(1,7),W(1,2),UV_GQQQ_1EPS,AMPL(2,131))
      CALL FFV1_0(W(1,11),W(1,7),W(1,2),UV_GQQQ_1EPS,AMPL(2,132))
      CALL FFV1_0(W(1,11),W(1,7),W(1,2),UV_GQQQ_1EPS,AMPL(2,133))
      CALL FFV1_0(W(1,11),W(1,7),W(1,2),UV_GQQB,AMPL(1,134))
      CALL FFV1_0(W(1,11),W(1,7),W(1,2),UV_GQQQ_1EPS,AMPL(2,135))
      CALL FFV1_0(W(1,11),W(1,7),W(1,2),UV_GQQT,AMPL(1,136))
      CALL FFV1_0(W(1,11),W(1,7),W(1,2),UV_GQQQ_1EPS,AMPL(2,137))
      CALL FFV1_0(W(1,11),W(1,7),W(1,2),UV_GQQG_1EPS,AMPL(2,138))
      CALL FFV1_0(W(1,11),W(1,7),W(1,2),R2_GQQ,AMPL(1,139))
C     Counter-term amplitude(s) for loop diagram number 78
      CALL FFV1_0(W(1,28),W(1,4),W(1,2),UV_GQQQ_1EPS,AMPL(2,140))
      CALL FFV1_0(W(1,28),W(1,4),W(1,2),UV_GQQQ_1EPS,AMPL(2,141))
      CALL FFV1_0(W(1,28),W(1,4),W(1,2),UV_GQQQ_1EPS,AMPL(2,142))
      CALL FFV1_0(W(1,28),W(1,4),W(1,2),UV_GQQQ_1EPS,AMPL(2,143))
      CALL FFV1_0(W(1,28),W(1,4),W(1,2),UV_GQQB,AMPL(1,144))
      CALL FFV1_0(W(1,28),W(1,4),W(1,2),UV_GQQQ_1EPS,AMPL(2,145))
      CALL FFV1_0(W(1,28),W(1,4),W(1,2),UV_GQQT,AMPL(1,146))
      CALL FFV1_0(W(1,28),W(1,4),W(1,2),UV_GQQQ_1EPS,AMPL(2,147))
      CALL FFV1_0(W(1,28),W(1,4),W(1,2),UV_GQQG_1EPS,AMPL(2,148))
      CALL FFV1_0(W(1,28),W(1,4),W(1,2),R2_GQQ,AMPL(1,149))
      CALL FFV1_2(W(1,11),W(1,2),GC_5,MDL_MT,MDL_WT,W(1,30))
C     Counter-term amplitude(s) for loop diagram number 80
      CALL R2_QQ_1_R2_QQ_2_0(W(1,30),W(1,7),R2_QQQ,R2_QQT,AMPL(1,150))
      CALL R2_QQ_2_0(W(1,30),W(1,7),UV_TMASS,AMPL(1,151))
      CALL R2_QQ_2_0(W(1,30),W(1,7),UV_TMASS_1EPS,AMPL(2,152))
      CALL FFV1_1(W(1,7),W(1,2),GC_5,MDL_MT,MDL_WT,W(1,31))
C     Counter-term amplitude(s) for loop diagram number 81
      CALL R2_QQ_1_R2_QQ_2_0(W(1,11),W(1,31),R2_QQQ,R2_QQT,AMPL(1,153))
      CALL R2_QQ_2_0(W(1,11),W(1,31),UV_TMASS,AMPL(1,154))
      CALL R2_QQ_2_0(W(1,11),W(1,31),UV_TMASS_1EPS,AMPL(2,155))
C     Counter-term amplitude(s) for loop diagram number 83
      CALL FFS1_5_0(W(1,30),W(1,4),W(1,3),UV_X0TT_H,UV_X0TT_A,AMPL(1
     $ ,156))
      CALL FFS1_5_0(W(1,30),W(1,4),W(1,3),UV_X0TT_H_1EPS
     $ ,UV_X0TT_A_1EPS,AMPL(2,157))
      CALL FFS1_5_0(W(1,30),W(1,4),W(1,3),R2_TTX0_H,R2_TTX0_A,AMPL(1
     $ ,158))
C     Counter-term amplitude(s) for loop diagram number 86
      CALL VVS1_0(W(1,1),W(1,18),W(1,3),R2_GGX0_HT,AMPL(1,159))
C     Counter-term amplitude(s) for loop diagram number 88
      CALL FFV1_0(W(1,5),W(1,29),W(1,1),UV_GQQQ_1EPS,AMPL(2,160))
      CALL FFV1_0(W(1,5),W(1,29),W(1,1),UV_GQQQ_1EPS,AMPL(2,161))
      CALL FFV1_0(W(1,5),W(1,29),W(1,1),UV_GQQQ_1EPS,AMPL(2,162))
      CALL FFV1_0(W(1,5),W(1,29),W(1,1),UV_GQQQ_1EPS,AMPL(2,163))
      CALL FFV1_0(W(1,5),W(1,29),W(1,1),UV_GQQB,AMPL(1,164))
      CALL FFV1_0(W(1,5),W(1,29),W(1,1),UV_GQQQ_1EPS,AMPL(2,165))
      CALL FFV1_0(W(1,5),W(1,29),W(1,1),UV_GQQT,AMPL(1,166))
      CALL FFV1_0(W(1,5),W(1,29),W(1,1),UV_GQQQ_1EPS,AMPL(2,167))
      CALL FFV1_0(W(1,5),W(1,29),W(1,1),UV_GQQG_1EPS,AMPL(2,168))
      CALL FFV1_0(W(1,5),W(1,29),W(1,1),R2_GQQ,AMPL(1,169))
C     Counter-term amplitude(s) for loop diagram number 89
      CALL FFV1_0(W(1,8),W(1,12),W(1,1),UV_GQQQ_1EPS,AMPL(2,170))
      CALL FFV1_0(W(1,8),W(1,12),W(1,1),UV_GQQQ_1EPS,AMPL(2,171))
      CALL FFV1_0(W(1,8),W(1,12),W(1,1),UV_GQQQ_1EPS,AMPL(2,172))
      CALL FFV1_0(W(1,8),W(1,12),W(1,1),UV_GQQQ_1EPS,AMPL(2,173))
      CALL FFV1_0(W(1,8),W(1,12),W(1,1),UV_GQQB,AMPL(1,174))
      CALL FFV1_0(W(1,8),W(1,12),W(1,1),UV_GQQQ_1EPS,AMPL(2,175))
      CALL FFV1_0(W(1,8),W(1,12),W(1,1),UV_GQQT,AMPL(1,176))
      CALL FFV1_0(W(1,8),W(1,12),W(1,1),UV_GQQQ_1EPS,AMPL(2,177))
      CALL FFV1_0(W(1,8),W(1,12),W(1,1),UV_GQQG_1EPS,AMPL(2,178))
      CALL FFV1_0(W(1,8),W(1,12),W(1,1),R2_GQQ,AMPL(1,179))
C     Counter-term amplitude(s) for loop diagram number 92
      CALL VVS1_0(W(1,1),W(1,19),W(1,3),R2_GGX0_HT,AMPL(1,180))
C     Counter-term amplitude(s) for loop diagram number 93
      CALL FFV1_0(W(1,10),W(1,7),W(1,1),UV_GQQQ_1EPS,AMPL(2,181))
      CALL FFV1_0(W(1,10),W(1,7),W(1,1),UV_GQQQ_1EPS,AMPL(2,182))
      CALL FFV1_0(W(1,10),W(1,7),W(1,1),UV_GQQQ_1EPS,AMPL(2,183))
      CALL FFV1_0(W(1,10),W(1,7),W(1,1),UV_GQQQ_1EPS,AMPL(2,184))
      CALL FFV1_0(W(1,10),W(1,7),W(1,1),UV_GQQB,AMPL(1,185))
      CALL FFV1_0(W(1,10),W(1,7),W(1,1),UV_GQQQ_1EPS,AMPL(2,186))
      CALL FFV1_0(W(1,10),W(1,7),W(1,1),UV_GQQT,AMPL(1,187))
      CALL FFV1_0(W(1,10),W(1,7),W(1,1),UV_GQQQ_1EPS,AMPL(2,188))
      CALL FFV1_0(W(1,10),W(1,7),W(1,1),UV_GQQG_1EPS,AMPL(2,189))
      CALL FFV1_0(W(1,10),W(1,7),W(1,1),R2_GQQ,AMPL(1,190))
C     Counter-term amplitude(s) for loop diagram number 97
      CALL FFV1_0(W(1,25),W(1,4),W(1,1),UV_GQQQ_1EPS,AMPL(2,191))
      CALL FFV1_0(W(1,25),W(1,4),W(1,1),UV_GQQQ_1EPS,AMPL(2,192))
      CALL FFV1_0(W(1,25),W(1,4),W(1,1),UV_GQQQ_1EPS,AMPL(2,193))
      CALL FFV1_0(W(1,25),W(1,4),W(1,1),UV_GQQQ_1EPS,AMPL(2,194))
      CALL FFV1_0(W(1,25),W(1,4),W(1,1),UV_GQQB,AMPL(1,195))
      CALL FFV1_0(W(1,25),W(1,4),W(1,1),UV_GQQQ_1EPS,AMPL(2,196))
      CALL FFV1_0(W(1,25),W(1,4),W(1,1),UV_GQQT,AMPL(1,197))
      CALL FFV1_0(W(1,25),W(1,4),W(1,1),UV_GQQQ_1EPS,AMPL(2,198))
      CALL FFV1_0(W(1,25),W(1,4),W(1,1),UV_GQQG_1EPS,AMPL(2,199))
      CALL FFV1_0(W(1,25),W(1,4),W(1,1),R2_GQQ,AMPL(1,200))
C     Counter-term amplitude(s) for loop diagram number 100
      CALL VVV1_0(W(1,1),W(1,2),W(1,13),R2_3GQ,AMPL(1,201))
C     Counter-term amplitude(s) for loop diagram number 102
      CALL VVV1_0(W(1,1),W(1,2),W(1,14),R2_3GQ,AMPL(1,202))
C     Counter-term amplitude(s) for loop diagram number 108
      CALL FFV1_0(W(1,5),W(1,31),W(1,1),UV_GQQQ_1EPS,AMPL(2,203))
      CALL FFV1_0(W(1,5),W(1,31),W(1,1),UV_GQQQ_1EPS,AMPL(2,204))
      CALL FFV1_0(W(1,5),W(1,31),W(1,1),UV_GQQQ_1EPS,AMPL(2,205))
      CALL FFV1_0(W(1,5),W(1,31),W(1,1),UV_GQQQ_1EPS,AMPL(2,206))
      CALL FFV1_0(W(1,5),W(1,31),W(1,1),UV_GQQB,AMPL(1,207))
      CALL FFV1_0(W(1,5),W(1,31),W(1,1),UV_GQQQ_1EPS,AMPL(2,208))
      CALL FFV1_0(W(1,5),W(1,31),W(1,1),UV_GQQT,AMPL(1,209))
      CALL FFV1_0(W(1,5),W(1,31),W(1,1),UV_GQQQ_1EPS,AMPL(2,210))
      CALL FFV1_0(W(1,5),W(1,31),W(1,1),UV_GQQG_1EPS,AMPL(2,211))
      CALL FFV1_0(W(1,5),W(1,31),W(1,1),R2_GQQ,AMPL(1,212))
C     Counter-term amplitude(s) for loop diagram number 113
      CALL FFV1_0(W(1,27),W(1,4),W(1,1),UV_GQQQ_1EPS,AMPL(2,213))
      CALL FFV1_0(W(1,27),W(1,4),W(1,1),UV_GQQQ_1EPS,AMPL(2,214))
      CALL FFV1_0(W(1,27),W(1,4),W(1,1),UV_GQQQ_1EPS,AMPL(2,215))
      CALL FFV1_0(W(1,27),W(1,4),W(1,1),UV_GQQQ_1EPS,AMPL(2,216))
      CALL FFV1_0(W(1,27),W(1,4),W(1,1),UV_GQQB,AMPL(1,217))
      CALL FFV1_0(W(1,27),W(1,4),W(1,1),UV_GQQQ_1EPS,AMPL(2,218))
      CALL FFV1_0(W(1,27),W(1,4),W(1,1),UV_GQQT,AMPL(1,219))
      CALL FFV1_0(W(1,27),W(1,4),W(1,1),UV_GQQQ_1EPS,AMPL(2,220))
      CALL FFV1_0(W(1,27),W(1,4),W(1,1),UV_GQQG_1EPS,AMPL(2,221))
      CALL FFV1_0(W(1,27),W(1,4),W(1,1),R2_GQQ,AMPL(1,222))
C     Counter-term amplitude(s) for loop diagram number 118
      CALL VVS1_0(W(1,1),W(1,20),W(1,3),R2_GGX0_HT,AMPL(1,223))
      CALL FFV1_1(W(1,12),W(1,1),GC_5,MDL_MT,MDL_WT,W(1,32))
C     Counter-term amplitude(s) for loop diagram number 126
      CALL R2_QQ_1_R2_QQ_2_0(W(1,8),W(1,32),R2_QQQ,R2_QQT,AMPL(1,224))
      CALL R2_QQ_2_0(W(1,8),W(1,32),UV_TMASS,AMPL(1,225))
      CALL R2_QQ_2_0(W(1,8),W(1,32),UV_TMASS_1EPS,AMPL(2,226))
      CALL FFV1_2(W(1,8),W(1,1),GC_5,MDL_MT,MDL_WT,W(1,33))
C     Counter-term amplitude(s) for loop diagram number 127
      CALL R2_QQ_1_R2_QQ_2_0(W(1,33),W(1,12),R2_QQQ,R2_QQT,AMPL(1,227))
      CALL R2_QQ_2_0(W(1,33),W(1,12),UV_TMASS,AMPL(1,228))
      CALL R2_QQ_2_0(W(1,33),W(1,12),UV_TMASS_1EPS,AMPL(2,229))
C     Counter-term amplitude(s) for loop diagram number 129
      CALL FFS1_5_0(W(1,5),W(1,32),W(1,3),UV_X0TT_H,UV_X0TT_A,AMPL(1
     $ ,230))
      CALL FFS1_5_0(W(1,5),W(1,32),W(1,3),UV_X0TT_H_1EPS
     $ ,UV_X0TT_A_1EPS,AMPL(2,231))
      CALL FFS1_5_0(W(1,5),W(1,32),W(1,3),R2_TTX0_H,R2_TTX0_A,AMPL(1
     $ ,232))
      CALL FFV1_2(W(1,10),W(1,1),GC_5,MDL_MT,MDL_WT,W(1,34))
C     Counter-term amplitude(s) for loop diagram number 132
      CALL R2_QQ_1_R2_QQ_2_0(W(1,34),W(1,7),R2_QQQ,R2_QQT,AMPL(1,233))
      CALL R2_QQ_2_0(W(1,34),W(1,7),UV_TMASS,AMPL(1,234))
      CALL R2_QQ_2_0(W(1,34),W(1,7),UV_TMASS_1EPS,AMPL(2,235))
      CALL FFV1_1(W(1,7),W(1,1),GC_5,MDL_MT,MDL_WT,W(1,35))
C     Counter-term amplitude(s) for loop diagram number 133
      CALL R2_QQ_1_R2_QQ_2_0(W(1,10),W(1,35),R2_QQQ,R2_QQT,AMPL(1,236))
      CALL R2_QQ_2_0(W(1,10),W(1,35),UV_TMASS,AMPL(1,237))
      CALL R2_QQ_2_0(W(1,10),W(1,35),UV_TMASS_1EPS,AMPL(2,238))
C     Counter-term amplitude(s) for loop diagram number 135
      CALL FFS1_5_0(W(1,34),W(1,4),W(1,3),UV_X0TT_H,UV_X0TT_A,AMPL(1
     $ ,239))
      CALL FFS1_5_0(W(1,34),W(1,4),W(1,3),UV_X0TT_H_1EPS
     $ ,UV_X0TT_A_1EPS,AMPL(2,240))
      CALL FFS1_5_0(W(1,34),W(1,4),W(1,3),R2_TTX0_H,R2_TTX0_A,AMPL(1
     $ ,241))
C     Counter-term amplitude(s) for loop diagram number 138
      CALL FFV1_0(W(1,5),W(1,35),W(1,2),UV_GQQQ_1EPS,AMPL(2,242))
      CALL FFV1_0(W(1,5),W(1,35),W(1,2),UV_GQQQ_1EPS,AMPL(2,243))
      CALL FFV1_0(W(1,5),W(1,35),W(1,2),UV_GQQQ_1EPS,AMPL(2,244))
      CALL FFV1_0(W(1,5),W(1,35),W(1,2),UV_GQQQ_1EPS,AMPL(2,245))
      CALL FFV1_0(W(1,5),W(1,35),W(1,2),UV_GQQB,AMPL(1,246))
      CALL FFV1_0(W(1,5),W(1,35),W(1,2),UV_GQQQ_1EPS,AMPL(2,247))
      CALL FFV1_0(W(1,5),W(1,35),W(1,2),UV_GQQT,AMPL(1,248))
      CALL FFV1_0(W(1,5),W(1,35),W(1,2),UV_GQQQ_1EPS,AMPL(2,249))
      CALL FFV1_0(W(1,5),W(1,35),W(1,2),UV_GQQG_1EPS,AMPL(2,250))
      CALL FFV1_0(W(1,5),W(1,35),W(1,2),R2_GQQ,AMPL(1,251))
C     Counter-term amplitude(s) for loop diagram number 140
      CALL FFV1_0(W(1,33),W(1,4),W(1,2),UV_GQQQ_1EPS,AMPL(2,252))
      CALL FFV1_0(W(1,33),W(1,4),W(1,2),UV_GQQQ_1EPS,AMPL(2,253))
      CALL FFV1_0(W(1,33),W(1,4),W(1,2),UV_GQQQ_1EPS,AMPL(2,254))
      CALL FFV1_0(W(1,33),W(1,4),W(1,2),UV_GQQQ_1EPS,AMPL(2,255))
      CALL FFV1_0(W(1,33),W(1,4),W(1,2),UV_GQQB,AMPL(1,256))
      CALL FFV1_0(W(1,33),W(1,4),W(1,2),UV_GQQQ_1EPS,AMPL(2,257))
      CALL FFV1_0(W(1,33),W(1,4),W(1,2),UV_GQQT,AMPL(1,258))
      CALL FFV1_0(W(1,33),W(1,4),W(1,2),UV_GQQQ_1EPS,AMPL(2,259))
      CALL FFV1_0(W(1,33),W(1,4),W(1,2),UV_GQQG_1EPS,AMPL(2,260))
      CALL FFV1_0(W(1,33),W(1,4),W(1,2),R2_GQQ,AMPL(1,261))
C     Counter-term amplitude(s) for loop diagram number 143
      CALL VVS1_0(W(1,2),W(1,21),W(1,3),R2_GGX0_HT,AMPL(1,262))
C     Counter-term amplitude(s) for loop diagram number 159
      CALL R2_GG_1_R2_GG_2_0(W(1,6),W(1,13),R2_GGG_1,R2_GGG_2,AMPL(1
     $ ,263))
C     Counter-term amplitude(s) for loop diagram number 160
      CALL R2_GG_1_R2_GG_2_0(W(1,6),W(1,14),R2_GGG_1,R2_GGG_2,AMPL(1
     $ ,264))
C     Counter-term amplitude(s) for loop diagram number 161
      CALL VVV1_0(W(1,1),W(1,2),W(1,13),R2_3GG,AMPL(1,265))
C     Counter-term amplitude(s) for loop diagram number 162
      CALL VVV1_0(W(1,1),W(1,2),W(1,14),R2_3GG,AMPL(1,266))
C     At this point, all CT amps needed for (QCD=6), i.e. of split
C      order ID=1, are computed.
      IF(FILTER_SO.AND.SQSO_TARGET.EQ.1) GOTO 2000

      GOTO 1001
 2000 CONTINUE
      CT_REQ_SO_DONE=.TRUE.
 1001 CONTINUE
      END

