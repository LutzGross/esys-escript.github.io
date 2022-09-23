C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C $Id: chrcr.f,v 1.1 1993/07/16 16:46:17 gdsjaar Exp $
C $Log: chrcr.f,v $
C Revision 1.1  1993/07/16 16:46:17  gdsjaar
C Changed plt to library rather than single source file.
C
C=======================================================================
      LOGICAL FUNCTION CHRCR(LINE,REAL)
      CHARACTER LINE* (*),FORM*10,CL*2

      CALL CHRTRM(LINE,LL)
      CALL CHRIC(LL,CL,NL)
      FORM = '(f'//CL(1:NL)//'.0)'
      READ (LINE(1:LL),FORM,ERR=10) REAL
      CHRCR = .TRUE.
      RETURN

   10 CHRCR = .FALSE.
      RETURN

      END
