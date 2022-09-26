C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C $Log: grview.f,v $
C Revision 1.2  2009/03/25 12:36:45  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.1  1994/04/07 20:02:58  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:52:05  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE GRVIEW (ASPECT, DVIEW0, DVIEW)
C=======================================================================

C   --*** GRVIEW *** (GRPLIB) Set graphics window
C   --   Written by Amy Gilkey - revised 02/20/86
C   --
C   --GRVIEW sets the device coordinates for the window given the X:Y
C   --axis length ratio.
C   --
C   --Parameters:
C   --   ASPECT - IN - the ratio of the X axis length to the Y axis length
C   --      = (X axis length / Y axis length) {in window units}
C   --      / (X axis length / Y axis length) {in device units}
C   --   DVIEW0 - IN - the maximum window (left, right, bottom, top)
C   --      (in device coordinates)
C   --   DVIEW - OUT - the actual window (left, right, bottom, top)
C   --      (in device coordinates)

      PARAMETER (KLFT=1, KRGT=2, KBOT=3, KTOP=4)

      REAL ASPECT
      REAL DVIEW0(KTOP), DVIEW(KTOP)

      IF (ASPECT .GE. 1.0) THEN
         DVIEW(KLFT) = DVIEW0(KLFT)
         DVIEW(KRGT) = DVIEW0(KRGT)
      ELSE
         OS = .5 * (1.0 - ASPECT) * (DVIEW0(KRGT) - DVIEW0(KLFT))
         DVIEW(KLFT) = DVIEW0(KLFT) + OS
         DVIEW(KRGT) = DVIEW0(KRGT) - OS
      END IF
      IF (ASPECT .LE. 1.0) THEN
         DVIEW(KBOT) = DVIEW0(KBOT)
         DVIEW(KTOP) = DVIEW0(KTOP)
      ELSE
         OS = .5 * (1.0 - 1.0/ASPECT) * (DVIEW0(KTOP) - DVIEW0(KBOT))
         DVIEW(KBOT) = DVIEW0(KBOT) + OS
         DVIEW(KTOP) = DVIEW0(KTOP) - OS
      END IF

      RETURN
      END
