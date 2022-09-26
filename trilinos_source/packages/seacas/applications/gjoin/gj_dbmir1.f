C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C $Id: dbmir1.f,v 1.4 2006/02/13 20:01:58 gdsjaar Exp $
C=======================================================================
      SUBROUTINE DBMIR1 (IELB, NUMELB, NUMLNK, LINK, NAME, NDIM, NONQUD)
C=======================================================================

C   --*** DBMIR1 *** (GJOIN) Fixup element connectivity for reflections
C   --   Written by Greg Sjaardema - revised 02/10/89
C   --   Modified from DBOEB1 Written by Amy Gilkey
C   --
C   --Parameters:
C   --   IELB - IN - the element block number
C   --   NUMELB - IN - the number of elements in the block
C   --   NUMLNK - IN - the number of nodes per element
C   --   LINK - IN/OUT - the element connectivity for this block
C   --   NAME - IN - the element type for this block
C   --   NDIM - IN - the spatial dimension (1,2,3)
C   --

      include 'exodusII.inc'

      INTEGER NUMELB, NUMLNK
      INTEGER LINK(NUMLNK,*)
      CHARACTER*(MXSTLN) NAME
      LOGICAL NONQUD

      CHARACTER*132 STRING

      IF (NUMELB .GT. 0) THEN

C...8-node Hexes
        IF ((NUMLNK .EQ. 8) .AND. (NDIM .EQ. 3) .AND.
     *    NAME(:3) .EQ. 'HEX') THEN
          DO 10 NE = 1, NUMELB
            ILTMP = LINK (2,NE)
            LINK(2,NE) = LINK(4,NE)
            LINK(4,NE) = ILTMP

            ILTMP = LINK (6,NE)
            LINK(6,NE) = LINK(8,NE)
            LINK(8,NE) = ILTMP
 10       CONTINUE

C...20-node Hexes
        ELSE IF ((NUMLNK .EQ. 20) .AND. (NDIM .EQ. 3) .AND.
     *    NAME(:3) .EQ. 'HEX') THEN
          DO 15 NE = 1, NUMELB
            ILTMP = LINK (2,NE)
            LINK(2,NE) = LINK(4,NE)
            LINK(4,NE) = ILTMP

            ILTMP = LINK (6,NE)
            LINK(6,NE) = LINK(8,NE)
            LINK(8,NE) = ILTMP

            ILTMP = LINK ( 9,NE)
            LINK( 9,NE) = LINK(12,NE)
            LINK(12,NE) = ILTMP

            ILTMP = LINK (10,NE)
            LINK(10,NE) = LINK(11,NE)
            LINK(11,NE) = ILTMP

            ILTMP = LINK (14,NE)
            LINK(14,NE) = LINK(16,NE)
            LINK(16,NE) = ILTMP

            ILTMP = LINK (17,NE)
            LINK(17,NE) = LINK(20,NE)
            LINK(20,NE) = ILTMP

            ILTMP = LINK (18,NE)
            LINK(18,NE) = LINK(19,NE)
            LINK(19,NE) = ILTMP
 15       CONTINUE

C...Quads/Shells
        ELSE IF ((NUMLNK .EQ. 4) .AND.
     *      (NAME(:4) .EQ. 'QUAD' .OR. NAME(:5) .EQ. 'SHELL')) THEN
          if (name(:5) .eq. 'SHELL') NONQUD = .TRUE.
          DO 20 NE = 1, NUMELB
            ILTMP = LINK (2,NE)
            LINK(2,NE) = LINK(4,NE)
            LINK(4,NE) = ILTMP
 20       CONTINUE

C...four-node tets...
        ELSE IF ((NUMLNK .EQ. 4) .AND.
     *      (NAME(:3) .EQ. 'TET')) THEN
          NONQUD = .TRUE.
          DO 25 NE = 1, NUMELB
            ILTMP = LINK (2,NE)
            LINK(2,NE) = LINK(3,NE)
            LINK(3,NE) = ILTMP
 25       CONTINUE

C...Bars
        ELSE IF (NUMLNK .EQ. 2) THEN
          NONQUD = .TRUE.
          DO 30 NE = 1, NUMELB
            ILTMP = LINK (2,NE)
            LINK(2,NE) = LINK(1,NE)
            LINK(1,NE) = ILTMP
 30       CONTINUE

C...Triangles
        ELSE IF (NUMLNK .EQ. 3) THEN
          NONQUD = .TRUE.
          DO 40 NE = 1, NUMELB
            ILTMP = LINK (2,NE)
            LINK(2,NE) = LINK(3,NE)
            LINK(3,NE) = ILTMP
 40       CONTINUE

C...6-node Triangles
        ELSE IF (NUMLNK .EQ. 6) THEN
          NONQUD = .TRUE.
          DO 50 NE = 1, NUMELB
            ILTMP = LINK (2,NE)
            LINK(2,NE) = LINK(3,NE)
            LINK(3,NE) = ILTMP
            ILTMP = LINK (4,NE)
            LINK(4,NE) = LINK(6,NE)
            LINK(6,NE) = ILTMP
 50       CONTINUE
        ELSE
          NONQUD = .TRUE.
          WRITE (STRING, 100) IELB, NUMLNK, NAME
 100      FORMAT('Element block ',I5,' contains ',I2,'-node ',A,
     *      ' elements which are not supported for mirroring by gjoin2')
          CALL SQZSTR (STRING, LSTR)
          CALL PRTERR ('PROGRAM', STRING(:LSTR))
        END IF
      END IF
      RETURN
      END
