$PROBLEM tte validation: repeated events, Weibull hazard in the time since the last event
; TLAST is the time of the last event before each record.  Data items in $DES
; come from the record that ends the interval, so TLAST holds over it.
$INPUT ID TIME EVID AMT CMT DV MDV DVID TLAST
$DATA gapTime.csv IGNORE=@
$SUBROUTINES ADVAN13 TOL=9
$MODEL COMP=(CHZ)
$PK
  TE = EXP(THETA(1) + ETA(1))
  SH = EXP(THETA(2))
$DES
  TT = T - TLAST
  IF (TT.LT.1E-10) TT = 1E-10
  DADT(1) = SH/TE*(TT/TE)**(SH - 1)
$ERROR
  IF (NEWIND.NE.2) OLDCHZ = 0
  CHZ = A(1) - OLDCHZ
  IF (MDV.EQ.0) OLDCHZ = A(1)
  TT = TIME - TLAST
  IF (TT.LT.1E-10) TT = 1E-10
  HAZ = SH/TE*(TT/TE)**(SH - 1)
  Y = 2*CHZ
  IF (DV.EQ.1) Y = 2*CHZ - 2*LOG(HAZ)
$THETA 1.386294 ; lte
$THETA 0.4054651 ; lsh
$OMEGA 0.1 ; eta.te
$ESTIMATION METHOD=CONDITIONAL LAPLACE -2LL MAXEVAL=9999 PRINT=5 NOABORT
$COVARIANCE PRINT=E UNCONDITIONAL
