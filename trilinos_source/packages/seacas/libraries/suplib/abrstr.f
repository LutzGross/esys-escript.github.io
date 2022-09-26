C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE ABRSTR (RETWRD, ABBR, STRTBL)
C=======================================================================
C$Id: abrstr.f,v 1.2 2009/03/25 12:46:01 gdsjaar Exp $
C$Log: abrstr.f,v $
CRevision 1.2  2009/03/25 12:46:01  gdsjaar
CAdd copyright and license notice to all files.
C
CRevision 1.1.1.1  1990/08/14 16:11:58  gdsjaar
CTesting
C
c Revision 1.1  90/08/14  16:11:57  gdsjaar
c Initial revision
c
c Revision 1.1  90/08/09  13:39:03  gdsjaar
c Initial revision
c

C   --*** ABRSTR *** (STRLIB) Find abbreviation for string
C   --   Written by Amy Gilkey - revised 02/14/86
C   --
C   --ABRSTR returns the non-abbreviated form of the given abbreviation
C   --from the list of possible strings.  The abbreviation must either
C   --be a complete string or it must only match one string.
C   --
C   --Parameters:
C   --   RETWRD - OUT - the string for the abbreviation; ' ' if none
C   --   ABBR - IN - the abbreviation
C   --   STRTBL - IN - the table of possible strings; ended by ' '

      CHARACTER*(*) RETWRD
      CHARACTER*(*) ABBR
      CHARACTER*(*) STRTBL(*)

      RETWRD = ' '

      IF (ABBR .EQ. ' ') RETURN

      L = INDEX (ABBR, ' ') - 1
      IF (L .LT. 0) L = LEN(ABBR)

      NFOUND = 0
      I = 1
  100 CONTINUE
      IF (STRTBL(I) .NE. ' ') THEN
         IF (ABBR .EQ. STRTBL(I)(1:L)) THEN
            RETWRD = STRTBL(I)
            IF (ABBR .EQ. STRTBL(I)) GOTO 110
            NFOUND = NFOUND + 1
         END IF
         I = I + 1
         GOTO 100
      END IF

      IF (NFOUND .GT. 1) RETWRD = ' '

  110 CONTINUE
      RETURN
      END
