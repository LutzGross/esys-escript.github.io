C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C    
C    See packages/seacas/LICENSE for details
      SUBROUTINE SHFTI(IARAY,N,M,ISTRT,IEND,INT)
C
C     THIS SUBROUTINE SHIFTS THE ROWS IN AN INTEGER ARRAY. IF INT>0 THEN
C     ALL ROWS ISTRT TO IEND ARE SHIFTED UP "INT" ROWS. IF INT<0 THE ALL
C     ROWS ISTRT TO IEND ARE SHIFTED DOWN "INT" ROWS.
C
      DIMENSION IARAY(N,M)
C
C  CALCULATE RANGE AND INCREMENT OF DO LOOP
C
      IF(INT.LT.0)THEN
C
C  SHIFT DOWN
C
         I1=IEND
         I2=ISTRT
         ID=-1
      ELSE IF(INT.GT.0)THEN
C
C  SHIFT UP
C
         I1=ISTRT
         I2=IEND
         ID=1
      ELSE
         RETURN
      END IF
C
C  PERFORM SHIFT
C
      DO J=1,M
         DO I=I1,I2,ID
            IARAY(I-INT,J)=IARAY(I,J)
         end do
      end do
      RETURN
      END
