C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C $Id: mxzero.f,v 1.3 1993/07/19 18:08:46 gdsjaar Exp $
C $Log: mxzero.f,v $
C Revision 1.3  1993/07/19 18:08:46  gdsjaar
C Added special case for n=4 since that is how plt calls it primarily
C
c Revision 1.2  1993/07/16  19:35:48  gdsjaar
c Restructured to optimize faster
c
c Revision 1.1  1993/07/16  16:47:37  gdsjaar
c Changed plt to library rather than single source file.
c
C=======================================================================
      SUBROUTINE MXZERO(N,MAT)
      REAL MAT(N,*)

      IF (N .EQ. 4) THEN
        MAT(1,1) = 0.0
        MAT(1,2) = 0.0
        MAT(1,3) = 0.0
        MAT(1,4) = 0.0

        MAT(2,1) = 0.0
        MAT(2,2) = 0.0
        MAT(2,3) = 0.0
        MAT(2,4) = 0.0

        MAT(3,1) = 0.0
        MAT(3,2) = 0.0
        MAT(3,3) = 0.0
        MAT(3,4) = 0.0

        MAT(4,1) = 0.0
        MAT(4,2) = 0.0
        MAT(4,3) = 0.0
        MAT(4,4) = 0.0
      ELSE
        DO 3030 I= 1,N
          DO 3050 J = 1,N
            MAT(I,J) = 0.
 3050     CONTINUE
 3030   CONTINUE
      end if
      RETURN

      END
