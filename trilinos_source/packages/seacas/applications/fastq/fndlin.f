C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C    
C    See packages/seacas/LICENSE for details

CC* FILE: [.QMESH]FNDLIN.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE FNDLIN_FQ (MXND, LXN, NODE1, NODE2, LINE, ERR)
C
C***********************************************************************
C
C  SUBROUTINE FNDLIN =  FINDS THE LINE WITH ENDS NODE1 & NODE2
C
C***********************************************************************
C
      DIMENSION LXN(4, MXND)
      DIMENSION LINES1(20), LINES2(20)
      LOGICAL ERR
C
      ERR = .FALSE.
C
      CALL GETLXN (MXND, LXN, NODE1, LINES1, NL1, ERR)
      CALL GETLXN (MXND, LXN, NODE2, LINES2, NL2, ERR)
C
      IF (.NOT.ERR) THEN
         ERR = .TRUE.
         DO 110 I = 1, NL1
            DO 100 J = 1, NL2
               IF (LINES1(I) .EQ. LINES2(J)) THEN
                  LINE = LINES1(I)
                  ERR = .FALSE.
                  RETURN
               END IF
  100       CONTINUE
  110    CONTINUE
      END IF
C
      RETURN
C
      END
