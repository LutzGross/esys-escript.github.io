C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C $Id: pltmaa.f,v 1.1 1993/07/16 16:48:46 gdsjaar Exp $
C $Log: pltmaa.f,v $
C Revision 1.1  1993/07/16 16:48:46  gdsjaar
C Changed plt to library rather than single source file.
C
C=======================================================================
      SUBROUTINE PLTMAA(ALT,AZI,UMAP)
      COMMON /CENBOD/XC,YC,ZC
      REAL ALT
      REAL AZI
      REAL UMAP(*)
      REAL R

      R = (UMAP(18)-XC)**2 + (UMAP(18+1)-YC)**2 + (UMAP(18+2)-ZC)**2
      R = SQRT(R)
      UMAP(21) = 1.
      UMAP(21+1) = 0.
      UMAP(21+2) = 0.
      UMAP(24) = 0.
      UMAP(24+1) = 1.
      UMAP(24+2) = 0.
      UMAP(27) = 0.
      UMAP(27+1) = 0.
      UMAP(27+2) = -1.
      CALL PLTRTX(- (90-ALT)*3.1415927/180.,UMAP)
      CALL PLTRTZ(AZI*3.1415927/180.,UMAP)
      UMAP(18) = XC - R*UMAP(27)
      UMAP(18+1) = YC - R*UMAP(27+1)
      UMAP(18+2) = ZC - R*UMAP(27+2)
      RETURN

      END
