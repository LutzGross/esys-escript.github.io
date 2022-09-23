C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C $Id: pltpnt.f,v 1.3 2000/10/25 18:55:02 gdsjaar Exp $
C $Log: pltpnt.f,v $
C Revision 1.3  2000/10/25 18:55:02  gdsjaar
C In the pltli? functions, check for N==0 before doing any array
C accesses.
C
C Also changed all references to 'mask' to be arrays where they were
C scalars since downstream code seems to treat them as arrays.
C
C Revision 1.2  1993/07/16 18:07:57  gdsjaar
C Added external pltblk statements so that linkers would pull in block
C data subprogram to initialize constants.
C
c Revision 1.1  1993/07/16  16:49:08  gdsjaar
c Changed plt to library rather than single source file.
c
C=======================================================================
      SUBROUTINE PLTPNT(N,X,Y)
      REAL DEVCAP(23)
      REAL DEFOUT(7)
      COMMON /STATUS/DEVCAP,DEFOUT
      REAL DEVP(5)
      COMMON /DEVICE/DEVP
      REAL COLP(3)
      REAL PALETT(3,16)
      COMMON /COLOR/COLP,PALETT
      REAL TEXTP(40)
      COMMON /TEXT/TEXTP
      REAL VECTP(5)
      REAL XCUR
      REAL YCUR
      COMMON /VECTRC/VECTP,XCUR,YCUR
      INTEGER IDEX(200,2)
      INTEGER NVECT(200,2)
      REAL XSIZE(200,2)
      REAL YSIZE(200,2)
      REAL X0(2300,2)
      REAL Y0(2300,2)
      REAL X1(2300,2)
      REAL Y1(2300,2)
      COMMON /FONT/IDEX,NVECT,XSIZE,YSIZE,X0,Y0,X1,Y1
      REAL GRAPHP(100)
      COMMON /GRAPH/GRAPHP
      COMMON /MAPPAR/MAPP(11)
      REAL MAPP
      COMMON /STORAG/MEMORY(1000)
      REAL SAVLEN
      INTEGER IDSHSV
      COMMON /PLTSTY/SAVLEN,IDSHSV
      DIMENSION X(*),Y(*)
      EXTERNAL PLTBLK
      INTEGER MASK(1)

      DO 2060 I = 1,N
         IF (MAPP(10).EQ.1.) THEN
            CALL PLTP2D(X(I),Y(I),XP,YP)

         ELSE
            XP = X(I)
            YP = Y(I)
         END IF

         MASK(1) = -1
         CALL PLTVWP(MAPP(6),MAPP(8),1,MASK,XP,YP)
         IF (MASK(1).EQ.-1) THEN
            CALL VDPNTA(XP,YP)
         END IF

 2060 CONTINUE
      XCUR = X(N)
      YCUR = Y(N)
      IDSHSV = 0
      SAVLEN = 0.
      RETURN

      END
