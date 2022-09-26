C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C $Id: pltuwn.f,v 1.1 1993/07/16 16:49:46 gdsjaar Exp $
C $Log: pltuwn.f,v $
C Revision 1.1  1993/07/16 16:49:46  gdsjaar
C Changed plt to library rather than single source file.
C
C=======================================================================
      SUBROUTINE PLTUWN(GMAP)
      DIMENSION GMAP(*)

      GMAP(7) = GMAP(7) - .0002
      GMAP(8) = GMAP(8) - .0002
      GMAP(9) = GMAP(9) + .0006
      GMAP(10) = GMAP(10) - .0002
      GMAP(11) = GMAP(11) + .0006
      GMAP(12) = GMAP(12) + .0006
      GMAP(13) = GMAP(13) - .0002
      GMAP(14) = GMAP(14) + .0006
      RETURN

      END
