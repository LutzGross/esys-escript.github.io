C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C    
C    See packages/seacas/LICENSE for details

C $Id: minmax.f,v 1.1 1990/11/30 11:12:11 gdsjaar Exp $
C $Log: minmax.f,v $
C Revision 1.1  1990/11/30 11:12:11  gdsjaar
C Initial revision
C
C
CC* FILE: [.QMESH]MINMAX_FQ.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE MINMAX_FQ (NDIM, N, X, Y, XMIN, XMAX, YMIN, YMAX)
C***********************************************************************
C
C  SUBROUTINE MINMAX_FQ = COMPARES AND RECORDS X AND Y EXTREMES
C
C***********************************************************************
C
      DIMENSION X (NDIM), Y (NDIM)
      XMIN = X (1)
      XMAX = X (1)
      YMIN = Y (1)
      YMAX = Y (1)
C
      DO 100 I = 1, N
         XMIN = AMIN1 (X (I), XMIN)
         XMAX = AMAX1 (X (I), XMAX)
         YMIN = AMIN1 (Y (I), YMIN)
         YMAX = AMAX1 (Y (I), YMAX)
  100 CONTINUE
C
      RETURN
C
      END
