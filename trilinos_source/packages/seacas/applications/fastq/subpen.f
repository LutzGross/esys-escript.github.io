C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C    
C    See packages/seacas/LICENSE for details

C $Id: subpen.f,v 1.1 1990/11/30 11:16:52 gdsjaar Exp $
C $Log: subpen.f,v $
C Revision 1.1  1990/11/30 11:16:52  gdsjaar
C Initial revision
C
C
CC* FILE: [.QMESH]SUBPEN.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE SUBPEN (NPER, NEWPER, X, Y, NID, XSUB, YSUB, NIDSUB,
     &   NUM, M1, M2, IADD, ITRI, XCEN, YCEN)
C***********************************************************************
C
C  SUBROUTINE SUBPEN = PUTS A PENTAGON SUBREGION'S PERIMETER INTO THE
C                      NPERIM ARRAYS
C
C***********************************************************************
C
      DIMENSION X (NPER), Y (NPER), NID (NPER)
      DIMENSION XSUB (NPER), YSUB (NPER), NIDSUB (NPER)
C
C  PUT SIDE ONE AND TWO INTO THE PERIMETER LIST
C
      KOUNT = 0
      DO 100 I = 1, NUM + 1
         KOUNT = KOUNT + 1
         J = I + IADD
         IF (J .GT. NPER)J = J - NPER
         XSUB (KOUNT) = X (J)
         YSUB (KOUNT) = Y (J)
         NIDSUB (KOUNT) = NID (J)
  100 CONTINUE
C
C  PUT SIDE THREE INTO THE LIST
C
      XDIF = XCEN - XSUB (KOUNT)
      YDIF = YCEN - YSUB (KOUNT)
      XINT = XDIF / DBLE(M1)
      YINT = YDIF / DBLE(M1)
      DO 110 I = 1, M1 - 1
         KOUNT = KOUNT + 1
         XSUB (KOUNT) = XSUB (KOUNT - 1) + XINT
         YSUB (KOUNT) = YSUB (KOUNT - 1) + YINT
         NIDSUB (KOUNT) =  (ITRI * 100000) + M1 - I + 1
  110 CONTINUE
C
C  ENTER THE CENTER POINT
C
      KOUNT = KOUNT + 1
      XSUB (KOUNT) = XCEN
      YSUB (KOUNT) = YCEN
      NIDSUB (KOUNT) = 100000
C
C  PUT SIDE FOUR INTO THE LIST
C
      ITRI2 = ITRI + 2
      IF (ITRI2 .GT. 3)ITRI2 = ITRI2 - 3
      XDIF = X (IADD + 1) - XCEN
      YDIF = Y (IADD + 1) - YCEN
      XINT = XDIF / DBLE(M2)
      YINT = YDIF / DBLE(M2)
      DO 120 I = 1, M2 - 1
         KOUNT = KOUNT + 1
         XSUB (KOUNT) = XSUB (KOUNT - 1) + XINT
         YSUB (KOUNT) = YSUB (KOUNT - 1) + YINT
         NIDSUB (KOUNT) =  (100000 * ITRI2) + I + 1
  120 CONTINUE
C
      NEWPER = KOUNT
C
      RETURN
C
      END
