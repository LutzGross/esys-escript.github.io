C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C    
C    See packages/seacas/LICENSE for details

C $Id: gmdis2.f,v 1.1 1991/02/21 15:43:27 gdsjaar Exp $
C $Log: gmdis2.f,v $
C Revision 1.1  1991/02/21 15:43:27  gdsjaar
C Initial revision
C
      SUBROUTINE GMDIS2 (COORD, DIRCOS, MASSLV, NIQSLV, TDIS,
     *    NIQM, NIQS, DMAX, NUMNP)
      DIMENSION COORD (NUMNP,*), DIRCOS(4,*), MASSLV(2,*), NIQSLV(*),
     *    TDIS(3,*)
C
C... PICK MATCH BASED ON CLOSEST DISTANCE
C
      DO 30 IMAS = 1, NIQM
          X1 = COORD (MASSLV(1, IMAS),1)
          Y1 = COORD (MASSLV(1, IMAS),2)
C
          DO 10 ISLV = 1, NIQS
              X0 = COORD (NIQSLV(ISLV),1)
              Y0 = COORD (NIQSLV(ISLV),2)
C
              TDIS(3,ISLV) = (X0 - X1)**2 + (Y0 - Y1)**2
   10     CONTINUE
C
          TMIN = DMAX
          NMIN = 0
          DO 20 ISLV = 1, NIQS
              IF (TDIS(3,ISLV) .LE. TMIN) THEN
                  TMIN = TDIS(3,ISLV)
                  NMIN = ISLV
              END IF
   20     CONTINUE
C
          IF (NMIN .NE. 0) THEN
              MASSLV(2,IMAS) = NIQSLV(NMIN)
              DIRCOS(4,IMAS) = TDIS(3,NMIN)
          ELSE
              MASSLV(2,IMAS) = 0
          END IF
   30 CONTINUE
      CALL CULL2 (DIRCOS, MASSLV, NIQM)
C
C ... INITIALIZE TANGENT AND NORMAL DISTANCES
C
      DO 40 ISEG = 1, NIQM
          X1 = COORD (MASSLV(1, ISEG),1)
          Y1 = COORD (MASSLV(1, ISEG),2)
C
          X0 = COORD (MASSLV(2, ISEG),1)
          Y0 = COORD (MASSLV(2, ISEG),2)
C
C DIRCOS is average unit direction vector from surfaces at node
C
          DCS1 = DIRCOS (1, ISEG)
          DCS2 = DIRCOS (2, ISEG)
C
          T = -( DCS1 * (X1-X0) + DCS2 * (Y1-Y0) )
          DIRCOS(3,ISEG) = T
          DIRCOS(4,ISEG) = SQRT( ABS(DIRCOS(4,ISEG) - T**2) )
   40 CONTINUE
C
      RETURN
      END
