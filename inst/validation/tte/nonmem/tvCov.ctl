$PROBLEM tte validation: repeated events with a covariate that switches on during follow-up
; XNM is the covariate in force just before each record (X on the previous
; record): data items in $DES come from the record that ends the interval.
$INPUT ID TIME EVID AMT CMT DV MDV DVID X=DROP TSW=DROP XNM
$DATA tvCov.csv IGNORE=@
$SUBROUTINES ADVAN13 TOL=9
$MODEL COMP=(CHZ)
$PK
  LH = THETA(1) + ETA(1)
  BX = THETA(2)
$DES
  DADT(1) = EXP(LH + BX*XNM)
$ERROR
  IF (NEWIND.NE.2) OLDCHZ = 0
  CHZ = A(1) - OLDCHZ
  IF (MDV.EQ.0) OLDCHZ = A(1)
  HAZ = EXP(LH + BX*XNM)
  Y = 2*CHZ
  IF (DV.EQ.1) Y = 2*CHZ - 2*LOG(HAZ)
$THETA -1.609438 ; lh
$THETA 0 ; bx
$OMEGA 0.3 ; eta.h
$ESTIMATION METHOD=CONDITIONAL LAPLACE -2LL MAXEVAL=9999 PRINT=5 NOABORT
$COVARIANCE PRINT=E UNCONDITIONAL
