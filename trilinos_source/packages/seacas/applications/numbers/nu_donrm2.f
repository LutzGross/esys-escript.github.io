C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C    
C    See packages/seacas/LICENSE for details

C $Id: donrm2.f,v 1.1 1991/02/21 15:42:55 gdsjaar Exp $
C $Log: donrm2.f,v $
C Revision 1.1  1991/02/21 15:42:55  gdsjaar
C Initial revision
C
      SUBROUTINE DONRM2 (COORD, LTNESS, MAP, DIRCOS, TEMP,
     *    NSEG, NUMNIQ, NUMNP)
      DIMENSION COORD(NUMNP, *), LTNESS(2,*), MAP(*), DIRCOS(4,*),
     *    TEMP(2,*)
C
      DO 10 ISEG = 1, NSEG
          XI = COORD( LTNESS(1,ISEG),1 )
          YI = COORD( LTNESS(1,ISEG),2 )
C
          XJ = COORD( LTNESS(2,ISEG),1 )
          YJ = COORD( LTNESS(2,ISEG),2 )
C
          DX = XI - XJ
          DY = YI - YJ
          RMAG = SQRT ( DX**2 + DY**2)
C
          TEMP(1,ISEG) = -DY / RMAG
          TEMP(2,ISEG) =  DX / RMAG
C
   10 CONTINUE
C
      DO 20 I=1,NUMNIQ
          DIRCOS(1,I) = 0.0
          DIRCOS(2,I) = 0.0
   20 CONTINUE
C
      DO 40 ISEG = 1, NSEG
          DO 30 J = 1, 2
              MISEG = MAP( 2 * (ISEG-1) + J )
              DIRCOS(1,MISEG) = DIRCOS(1,MISEG) + TEMP(1,ISEG)
              DIRCOS(2,MISEG) = DIRCOS(2,MISEG) + TEMP(2,ISEG)
   30     CONTINUE
   40 CONTINUE
C
C ... NORMALIZE ALL DIRECTION COSINES
C
      DO 50 I = 1, NUMNIQ
          A = DIRCOS(1,I)
          B = DIRCOS(2,I)
          R = SQRT(A**2 + B**2)
          DIRCOS(1,I) = A/R
          DIRCOS(2,I) = B/R
   50 CONTINUE
      RETURN
      END
