C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C    
C    See packages/seacas/LICENSE for details

C $Id: selbox.f,v 1.2 1998/03/22 05:34:43 gdsjaar Exp $
C $Log: selbox.f,v $
C Revision 1.2  1998/03/22 05:34:43  gdsjaar
C General cleanp of unused variables. Reordered DATA statements in
C command.f so would compile with f2c.
C
C Revision 1.1.1.1  1991/02/21 15:45:24  gdsjaar
C NUMBERS: Greg Sjaardema, initial Unix release
C
c Revision 1.1  1991/02/21  15:45:23  gdsjaar
c Initial revision
c
      SUBROUTINE SELBOX (COORD, NUMNP, NDIM, P1, SELECT, NODEL)
      DIMENSION COORD (NUMNP,*), P1(*)
      LOGICAL SELECT(*)
      CHARACTER*8 NODEL
      CHARACTER*80 STRTMP
      INTEGER LENSTR
C
      CALL INILOG (NUMNP, .FALSE., SELECT)
      INUM = 0
      IF (NDIM .EQ. 2) THEN
      DO 10 I=1, NUMNP
            X0 = COORD(I,1)
            Y0 = COORD(I,2)
            IF (X0 .GE. P1(1) .AND. X0 .LE. P1(2) .AND.
     *          Y0 .GE. P1(3) .AND. Y0 .LE. P1(4)) THEN
               SELECT(I) = .TRUE.
               INUM = INUM + 1
         END IF
   10 CONTINUE
      ELSE IF (NDIM .EQ. 3) THEN
      DO 20 I=1, NUMNP
            X0 = COORD(I,1)
            Y0 = COORD(I,2)
            Z0 = COORD(I,3)
            IF (X0 .GE. P1(1) .AND. X0 .LE. P1(2) .AND.
     *          Y0 .GE. P1(3) .AND. Y0 .LE. P1(4) .AND.
     *          Z0 .GE. P1(5) .AND. Z0 .LE. P1(6)) THEN
               SELECT(I) = .TRUE.
               INUM = INUM + 1
         END IF
   20 CONTINUE
      ELSE
            CALL PRTERR ('PROGRAM', 'Illegal dimension in SELBOX')
      END IF

      IF (INUM .EQ. 0) THEN
         CALL PRTERR ('WARNING',
     *        'No '// NODEL(:LENSTR(NODEL))//' found in range.')
      ELSE
         WRITE (STRTMP, 100) INUM, NODEL
  100    FORMAT (I10,' ',A,' Selected.')
         CALL SQZSTR (STRTMP, LTMP)
         CALL PRTERR ('CMDSPEC', STRTMP(:LTMP))
      END IF
      RETURN
      END
