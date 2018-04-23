ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      written by the UFO converter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE COUP1()

      IMPLICIT NONE
      INCLUDE 'model_functions.inc'

      DOUBLE PRECISION PI, ZERO
      PARAMETER  (PI=3.141592653589793D0)
      PARAMETER  (ZERO=0D0)
      INCLUDE 'input.inc'
      INCLUDE 'coupl.inc'
      REAL*16 MP__PI, MP__ZERO
      PARAMETER (MP__PI=3.1415926535897932384626433832795E0_16)
      PARAMETER (MP__ZERO=0E0_16)
      INCLUDE 'mp_input.inc'
      INCLUDE 'mp_coupl.inc'

      GC_3009A = (MDL_SINA*MDL_KABB*MDL_YB)/MDL_SQRT__2
      MP__GC_3009A = (MP__MDL_SINA*MP__MDL_KABB*MP__MDL_YB)
     $ /MP__MDL_SQRT__2
      GC_3009H = -((MDL_COSA*MDL_COMPLEXI*MDL_KHBB*MDL_YB)/MDL_SQRT__2)
      MP__GC_3009H = -((MP__MDL_COSA*MP__MDL_COMPLEXI*MP__MDL_KHBB
     $ *MP__MDL_YB)/MP__MDL_SQRT__2)
      GC_3010A = (MDL_SINA*MDL_KATT*MDL_YT)/MDL_SQRT__2
      MP__GC_3010A = (MP__MDL_SINA*MP__MDL_KATT*MP__MDL_YT)
     $ /MP__MDL_SQRT__2
      GC_3010H = -((MDL_COSA*MDL_COMPLEXI*MDL_KHTT*MDL_YT)/MDL_SQRT__2)
      MP__GC_3010H = -((MP__MDL_COSA*MP__MDL_COMPLEXI*MP__MDL_KHTT
     $ *MP__MDL_YT)/MP__MDL_SQRT__2)
      END
