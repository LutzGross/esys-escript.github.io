C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C $Id: pltmmv.f,v 1.1 1993/07/16 16:48:52 gdsjaar Exp $
C $Log: pltmmv.f,v $
C Revision 1.1  1993/07/16 16:48:52  gdsjaar
C Changed plt to library rather than single source file.
C
C=======================================================================
      SUBROUTINE PLTMMV(DX,DY,DZ,UMAP)
      COMMON /CENBOD/XC,YC,ZC
      REAL DX,DY,DZ
      REAL UMAP(*)

      UMAP(18) = UMAP(18) + DX
      UMAP(18+1) = UMAP(18+1) + DY
      UMAP(18+2) = UMAP(18+2) + DZ
      XC = XC + DX
      YC = YC + DY
      ZC = ZC + DZ
      RETURN

      END
