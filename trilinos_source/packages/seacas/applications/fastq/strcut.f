C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C    
C    See packages/seacas/LICENSE for details

C $Id: strcut.f,v 1.1 1990/11/30 11:16:42 gdsjaar Exp $
C $Log: strcut.f,v $
C Revision 1.1  1990/11/30 11:16:42  gdsjaar
C Initial revision
C
C
      SUBROUTINE STRCUT (STRING)
C***********************************************************************
C
C  SUBROUTINE STRCUT = DELETES ALL PRECEDING BLANKS FROM A STRING
C
C***********************************************************************
C
      CHARACTER * (*) STRING, HOLD*80
      CALL STRIPB (STRING, ILEFT, IRIGHT)
      IF (IRIGHT .GT. ILEFT) THEN
         HOLD = STRING (ILEFT:)
         STRING = HOLD
      ELSE
         HOLD = ' '
      ENDIF
      RETURN
      END
