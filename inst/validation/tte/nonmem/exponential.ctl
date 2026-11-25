$PROBLEM tte validation: constant hazard, single event, right-censored
$INPUT ID TIME EVID AMT CMT DV MDV DVID
$DATA exponential.csv IGNORE=@
$SUBROUTINES ADVAN13 TOL=9
$MODEL COMP=(CHZ)
$PK
  LH = THETA(1) + ETA(1)
$DES
  DADT(1) = EXP(LH)
$ERROR
  ; each record scores the hazard accumulated since the previous record
  IF (NEWIND.NE.2) OLDCHZ = 0
  CHZ = A(1) - OLDCHZ
  IF (MDV.EQ.0) OLDCHZ = A(1)
  HAZ = EXP(LH)
  Y = 2*CHZ
  IF (DV.EQ.1) Y = 2*CHZ - 2*LOG(HAZ)
$THETA -2.302585 ; lh
$OMEGA 0 FIX ; no random effect
$ESTIMATION METHOD=ZERO -2LL MAXEVAL=9999 PRINT=5 NOABORT
$COVARIANCE PRINT=E UNCONDITIONAL
