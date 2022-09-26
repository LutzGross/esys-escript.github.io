C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C    
C    See packages/seacas/LICENSE for details

C $Id: line3.f,v 1.1 1991/02/21 15:43:53 gdsjaar Exp $
C $Log: line3.f,v $
C Revision 1.1  1991/02/21 15:43:53  gdsjaar
C Initial revision
C
      SUBROUTINE LINE3 (COORD, NUMNP, DIST, T, NDIM, P1, P2, TOLER,
     *   NODEL, BOUND, SORTYP, MAP, SORUP, INUM, OPT, SELECT)
      DIMENSION COORD (NUMNP,*), DIST(*), T(*), P1(*), P2(*),
     *   TOLER(2), MAP(*)
      CHARACTER*(*) NODEL, BOUND, SORTYP, OPT
      LOGICAL SORUP, SELECT(*), ISABRT
      include 'nu_io.blk'
C
      CALL LOCOUT ('LINE', NDIM, NODEL, TOLER, SORTYP, P1, P2, BOUND)
C
      IF (BOUND(:3) .EQ. 'BOU') THEN
         BMULT = 1.0
      ELSE
         BMULT = 0.0
      END IF
C
      TEMP = TOLER(1)
      TOLER(1) = MAX(0.0, TEMP - TOLER(2))
      TOLER(2) = MAX(0.0, TEMP + TOLER(2))
C
      A  = P2(1) - P1(1)
      B  = P2(2) - P1(2)
      C  = P2(3) - P1(3)
      X1 = P1(1)
      Y1 = P1(2)
      Z1 = P1(3)
      DLINE = A**2 + B**2 + C**2
      IF (DLINE .EQ. 0.0) THEN
         CALL PRTERR ('CMDERR', 'Zero length line input')
         RETURN
      END IF
C
      DO 10 I=1, NUMNP
         IF (SELECT(I)) THEN
            X0 = COORD(I,1)
            Y0 = COORD(I,2)
            Z0 = COORD(I,3)
            T(I) = -1. * (A * (X1 - X0) + B * (Y1 - Y0) + C * (Z1 - Z0))
     *              / (A**2 + B**2 + C**2)
C
            X = X1 + A * T(I)
            Y = Y1 + B * T(I)
            Z = Z1 + C * T(I)
C
            DIST(I) = (X - X0)**2 + (Y - Y0)**2 + (Z - Z0)**2
         END IF
   10 CONTINUE
      INUM = 0
      DISMIN = 1.0E38
      DO 20 I=1, NUMNP
         IF (SELECT(I)) THEN
            DISMIN = MIN(DIST(I), ABS(DISMIN-TEMP))
            IF (DIST(I) .GE. TOLER(1)**2 .AND. DIST(I) .LE. TOLER(2)**2
     *         .AND. BMULT * T(I) .GE. 0.0 .AND. BMULT * T(I) .LE. 1.0)
     *         THEN
               INUM = INUM + 1
               MAP(INUM) = I
            END IF
         END IF
   20 CONTINUE

      IF      (SORTYP .EQ. 'X') THEN
         CALL INDEXX (COORD(1,1), MAP, INUM, .FALSE.)
      ELSE IF (SORTYP .EQ. 'Y') THEN
         CALL INDEXX (COORD(1,2), MAP, INUM, .FALSE.)
      ELSE IF (SORTYP .EQ. 'Z') THEN
         CALL INDEXX (COORD(1,3), MAP, INUM, .FALSE.)
      ELSE IF (SORTYP .EQ. 'T' .OR. SORTYP .EQ. 'PARAMETR') THEN
         CALL INDEXX (T,          MAP, INUM, .FALSE.)
      ELSE IF (SORTYP .EQ. 'DISTANCE') THEN
         CALL INDEXX (DIST,       MAP, INUM, .FALSE.)
      END IF

      IF (SORUP) THEN
         IBEG = 1
         IEND = INUM
         IINC = 1
      ELSE
         IBEG = INUM
         IEND = 1
         IINC = -1
      END IF

      IF (OPT .EQ. '*' .OR. INDEX(OPT, 'P') .GT. 0) THEN
         DO 30 IO=IOMIN, IOMAX
            WRITE (IO, 40) NODEL
   30    CONTINUE
   40    FORMAT (/,T50,'DISTANCE',/2X,A8,T16,'X',T26,'Y',T36,'Z',
     *      T45,'NORMAL',T55,'PARAMETRIC',/)
         DO 60 IN = IBEG, IEND, IINC
            IF (ISABRT()) RETURN
            I = MAP(IN)
            DO 50 IO=IOMIN, IOMAX
               WRITE (IO, 90) I, (COORD(I,J),J=1,3), SQRT(DIST(I)),
     *            T(I)
   50       CONTINUE
   60    CONTINUE
C
         IF (INUM .EQ. 0) THEN
            DO 70 IO=IOMIN, IOMAX
               WRITE (IO, 80) SQRT(DISMIN)
   70       CONTINUE
         END IF
      END IF
   80 FORMAT (/' None found within range, minimum distance = ',
     *   1PE12.3,/)
   90 FORMAT (I10, 3(F10.4), 2(1PE12.3))
      RETURN
      END
