C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C=======================================================================
      LOGICAL FUNCTION FFMATC (IFLD, INTYP, CFIELD, MATCH, NLET)
C=======================================================================
C$Id: ffmatc.f,v 1.2 2009/03/25 12:46:02 gdsjaar Exp $
C$Log: ffmatc.f,v $
CRevision 1.2  2009/03/25 12:46:02  gdsjaar
CAdd copyright and license notice to all files.
C
CRevision 1.1.1.1  1990/08/14 16:14:28  gdsjaar
CTesting
C
c Revision 1.1  90/08/14  16:14:27  gdsjaar
c Initial revision
c
c Revision 1.1  90/08/09  13:39:25  gdsjaar
c Initial revision
c

C   --*** FFMATC *** (FFLIB) Parse free-field character string if match
C   --   Written by Amy Gilkey - revised 07/01/87
C   --
C   --FFMATC parses a character field and returns true (and increments
C   --IFLD) iff it is equal to the match string.  Only NLET letters must
C   --be in the input field to match, but if more letters are given, they
C   --must match the match string exactly.
C   --
C   --Parameters:
C   --   IFLD - IN/OUT - the index of the current field number, incremented
C   --      only if field matches
C   --   INTYP - IN - the input type from the free-field reader
C   --   CFIELD - IN - the character fields
C   --   MATCH - IN - the match string
C   --   NLET - IN - number of letters that must match

C   --Routines Called:
C   --   MATSTR - (STRLIB) Check if string matches

      INTEGER IFLD
      INTEGER INTYP(*)
      CHARACTER*(*) CFIELD(*)
      CHARACTER*(*) MATCH
      INTEGER NLET

      LOGICAL MATSTR

      IF (INTYP(IFLD) .GE. 0) THEN
         FFMATC = MATSTR (CFIELD(IFLD), MATCH, NLET)
         IF (FFMATC) IFLD = IFLD + 1
      ELSE
         FFMATC = .FALSE.
      END IF

      RETURN
      END
