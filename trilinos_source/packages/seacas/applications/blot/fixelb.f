C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C $Log: fixelb.f,v $
C Revision 1.3  2009/03/25 12:36:44  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.2  2009/01/22 21:34:21  gdsjaar
C There were several inline dbnums common blocks. Replaced with the
C include so they all have the same size with the added variable types.
C
C Added minor support for nodeset and sideset variables.
C
C It can print the count and the names, but that is all
C at this time.
C
C Revision 1.1  1994/04/07 20:01:02  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:50:42  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE FIXELB (IELBST, LENF, IF2EL, IF2EL2, IE2ELB, NEWELB)
C=======================================================================

C   --*** FIXELB *** (MESH) Turn on/off element blocks (3D)
C   --   Written by Amy Gilkey - revised 03/10/88
C   --
C   --FIXELB adjusts the face array to reflect the new ON element blocks.
C   --Faces in OFF element blocks are moved to the LENF(NELBLK+4) set.
C   --
C   --Parameters:
C   --   IELBST - IN - the element block status (>=0 if ON)
C   --   LENF - IN - the cumulative face counts by element block (3D only)
C   --   IF2EL - IN - the element number of each face
C   --   IF2EL2 - IN - the secondary element number of each face
C   --   IE2ELB - IN - the element block for each element
C   --   NEWELB - OUT - size = LENF(NELBLK+1)
C   --
C   --Common Variables:
C   --   Uses NELBLK of /DBNUMS/

      include 'dbnums.blk'
      COMMON /D3NUMS/ IS3DIM, NNPSUR, NUMNPF, LLNSET
      LOGICAL IS3DIM

      INTEGER IELBST(NELBLK)
      INTEGER LENF(0:NELBLK+4)
      INTEGER IF2EL(*), IF2EL2(*)
      INTEGER IE2ELB(NUMEL)
      INTEGER NEWELB(*)

C   --Check each face by its defining elements, and move face if necessary

      DO 110 IELB = 1, NELBLK+4
         DO 100 IFAC = LENF(IELB-1)+1, LENF(IELB)

C         --Determine the number of elements with ON element blocks

            NON = 0
            IF (IELBST(IE2ELB(IF2EL(IFAC))) .GE. 0) NON = NON + 1
            IEL = IF2EL2(IFAC)
            IF (IEL .GT. 0) THEN
               IF (IELBST(IE2ELB(IEL)) .GE. 0) NON = NON + 1
            END IF

            IF (NON .EQ. 0) THEN

C            --If none ON, change to an OFF face
               NEWELB(IFAC) = NELBLK+4

            ELSE IF (NON .EQ. 1) THEN

C            --If one face ON, change to a surface face
               IF (IELBST(IE2ELB(IF2EL(IFAC))) .GE. 0) THEN
                  NEWELB(IFAC) = IE2ELB(IF2EL(IFAC))
               ELSE
C               --Swap nodes to simulate surface being defined by facing element
                  NEWELB(IFAC) = - IE2ELB(IF2EL2(IFAC))
               END IF

            ELSE

C            --If both element blocks ON, change to an interior face
               NEWELB(IFAC) = NELBLK+1

            END IF

  100    CONTINUE
  110 CONTINUE

      RETURN
      END
