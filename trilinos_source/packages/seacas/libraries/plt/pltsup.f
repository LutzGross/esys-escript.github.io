C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C $Id: pltsup.f,v 1.1 1993/07/16 16:49:39 gdsjaar Exp $
C $Log: pltsup.f,v $
C Revision 1.1  1993/07/16 16:49:39  gdsjaar
C Changed plt to library rather than single source file.
C
C=======================================================================
      SUBROUTINE PLTSUP(YBUMP,YCHRSZ)
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

      YBUMP = TEXTP(1)*TEXTP(33)
      TEXTP(18) = TEXTP(18) - YBUMP*TEXTP(29)
      TEXTP(19) = TEXTP(19) + YBUMP*TEXTP(28)
      DO 2110 J = 14,17
         TEXTP(J) = TEXTP(J)*TEXTP(32)
 2110 CONTINUE
      YCHRSZ = TEXTP(1)*TEXTP(32)
      RETURN

      END
