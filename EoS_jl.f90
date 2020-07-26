SUBROUTINE ComputeDerivatives_Pressure & 
  ( D, E, Ne, nx, dPdD, & 
    dPdT, dPdY, dPdE, dPdDe, dPdTau, Units)

  USE KindModule, ONLY: &
    DP
  USE EquationOfStateModule_TABLE, ONLY: &
    ComputeSpecificInternalEnergy_TABLE, &
    ComputeAuxiliary_Fluid_TABLE, &
    ComputePressure_TABLE
  USE EquationOfStateModule_TABLE, only: &
    InitializeEquationOfState_TABLE, &
    FinalizeEquationOfState_TABLE
  USE UnitsModule, ONLY: &
    Gram, &
    Centimeter, &
    Kelvin, &
    Erg
  
  INTEGER,  INTENT(in)  :: nx
  REAL(DP), INTENT(in)  :: D(nx), E(nx), Ne(nx)
  REAL(DP) :: D2(nx), E2(nx), Ne2(nx) ! remove units
  
  LOGICAL,  INTENT(in)  :: Units ! if True, incoming values are in physical units

  REAL(DP), DIMENSION(nx) :: P, Cs_table
  REAL(DP), DIMENSION(nx) :: Tau, T, Y, Em, Gm, S

  REAL(DP), INTENT(out), DIMENSION(nx) :: dPdD, dPdT, dPdY
  REAL(DP), DIMENSION(nx) :: dEdD, dEdT, dEdY
  REAL(DP), INTENT(out), DIMENSION(nx) :: dPdE, dPdDe, dPdTau

  CHARACTER(128) :: EosTableName 

  EosTableName = 'wl-EOS-SFHo-25-50-100.h5'

  ! ========================================================
  
  IF ( Units ) THEN
    D2 = D * ( Gram / Centimeter**3 )
    E2 = E * ( Erg / Centimeter**3 )
    Ne2 = Ne * ( 1.0_DP / Centimeter**3 )
  ELSE
    D2 = D 
    E2 = E
    Ne2 = Ne 
  END IF
  
  CALL InitializeEquationOfState_TABLE &
    ( EquationOfStateTableName_Option &
        = TRIM( EosTableName ) )

  CALL ComputeAuxiliary_Fluid_TABLE &
        ( D2, E2, Ne2, P, T, Y, S, Em, Gm, Cs_table )
  
  CALL ComputeSpecificInternalEnergy_TABLE &
        ( D2, T, Y, Em, dEdD, dEdT, dEdY )

  CALL ComputePressure_TABLE &
        ( D2, T, Y, P, dPdD, dPdT, dPdY )
  
  Tau = 1.0_DP / D2
  
  dPdE   = ( dPdT / dEdT )
  dPdDe  = ( ( Tau ) * ( dPdY - dEdY * dPdE ) )
  dPdTau = ( (dPdDe * Y + dEdD * dPdE - dPdD) / (Tau**2) )

  IF ( Units ) THEN
    dPdE = dPdE / ( Gram / Centimeter**3 )
    dPdDe = dPdDe / ( Erg / Gram )
    dPdTau = dPdTau / ( Erg * Gram / Centimeter**6 )

    dPdD = dPdD / ( Erg / Gram )
    dPdT = dPdT / ( ( Erg / Centimeter**3 ) / Kelvin )
    dPdY = dPdY / ( Erg / Centimeter**3 )  
  END IF

  CALL FinalizeEquationOfState_TABLE

END SUBROUTINE ComputeDerivatives_Pressure

SUBROUTINE ComputeDerivatives_Pressure_Scalar & 
  ( D, E, Ne, dPdD, & 
    dPdT, dPdY, dPdE, dPdDe, dPdTau, Units)

  USE KindModule, ONLY: &
    DP
  USE EquationOfStateModule_TABLE, ONLY: &
    ComputeSpecificInternalEnergy_TABLE, &
    ComputeAuxiliary_Fluid_TABLE, &
    ComputePressure_TABLE
  USE EquationOfStateModule_TABLE, only: &
    InitializeEquationOfState_TABLE, &
    FinalizeEquationOfState_TABLE
  USE UnitsModule, ONLY: &
    Gram, &
    Centimeter, &
    Kelvin, &
    Erg

  REAL(DP), INTENT(in) :: D, E, Ne
  REAL(DP) :: D2, E2, Ne2 ! remove units

  LOGICAL,  INTENT(in)  :: Units ! if True, incoming values are in physical units

  REAL(DP) :: P, Cs_table
  REAL(DP) :: Tau, T, Y, Em, Gm, S

  REAL(DP), INTENT(out) :: dPdD, dPdT, dPdY
  REAL(DP) :: dEdD, dEdT, dEdY
  REAL(DP), INTENT(out) :: dPdE, dPdDe, dPdTau

  CHARACTER(128) :: EosTableName 

  EosTableName = 'wl-EOS-SFHo-25-50-100.h5'

  ! ========================================================

  IF ( Units ) THEN
    D2 = D * ( Gram / Centimeter**3 )
    E2 = E * ( Erg / Centimeter**3 )
    Ne2 = Ne * ( 1.0_DP / Centimeter**3 )
  ELSE
    D2 = D
    E2 = E
    Ne2 = Ne
  END IF
  
  CALL InitializeEquationOfState_TABLE &
    ( EquationOfStateTableName_Option &
        = TRIM( EosTableName ) )

  CALL ComputeAuxiliary_Fluid_TABLE &
        ( D2, E2, Ne2, P, T, Y, S, Em, Gm, Cs_table )
  
  CALL ComputeSpecificInternalEnergy_TABLE &
        ( D2, T, Y, Em, dEdD, dEdT, dEdY )

  CALL ComputePressure_TABLE &
        ( D2, T, Y, P, dPdD, dPdT, dPdY )

  Tau = 1.0_DP / D2

  dPdE   = ( dPdT / dEdT )
  dPdDe  = ( ( Tau ) * ( dPdY - dEdY * dPdE ) )
  dPdTau = ( (dPdDe * Y + dEdD * dPdE - dPdD) / (Tau**2) )

  IF ( Units ) THEN
    dPdE = dPdE / ( Gram / Centimeter**3 )
    dPdDe = dPdDe / ( Erg / Gram )
    dPdTau = dPdTau / ( Erg * Gram / Centimeter**6 )

    dPdD = dPdD / ( Erg / Gram )
    dPdT = dPdT / ( ( Erg / Centimeter**3 ) / Kelvin )
    dPdY = dPdY / ( Erg / Centimeter**3 ) 
  END IF

  CALL FinalizeEquationOfState_TABLE

END SUBROUTINE ComputeDerivatives_Pressure_Scalar

SUBROUTINE ComputeDerivatives_InternalEnergy & 
  ( D, E, Ne, nx, dEdD, dEdT, dEdY, Units)

  USE KindModule, ONLY: &
    DP
  USE EquationOfStateModule_TABLE, ONLY: &
    ComputeSpecificInternalEnergy_TABLE, &
    ComputeAuxiliary_Fluid_TABLE
  USE EquationOfStateModule_TABLE, only: &
    InitializeEquationOfState_TABLE, &
    FinalizeEquationOfState_TABLE
  USE UnitsModule, ONLY: &
    Gram, &
    Centimeter, &
    Kelvin, &
    Erg
  
  INTEGER,  INTENT(in)  :: nx
  REAL(DP), INTENT(in)  :: D(nx), E(nx), Ne(nx)
  REAL(DP) :: D2(nx), E2(nx), Ne2(nx) ! remove units
  
  LOGICAL,  INTENT(in)  :: Units ! if True, incoming values are in physical units

  REAL(DP), DIMENSION(nx) :: P, Cs_table
  REAL(DP), DIMENSION(nx) :: T, Y, Em, Gm, S

  REAL(DP), DIMENSION(nx) :: dPdD, dPdT, dPdY
  REAL(DP), DIMENSION(nx), INTENT(out) :: dEdD, dEdT, dEdY

  CHARACTER(128) :: EosTableName 

  EosTableName = 'wl-EOS-SFHo-25-50-100.h5'

  ! ========================================================
  
  IF ( Units ) THEN
    D2 = D * ( Gram / Centimeter**3 )
    E2 = E * ( Erg / Centimeter**3 )
    Ne2 = Ne * ( 1.0_DP / Centimeter**3 )
  ELSE
    D2 = D 
    E2 = E
    Ne2 = Ne 
  END IF
  
  CALL InitializeEquationOfState_TABLE &
    ( EquationOfStateTableName_Option &
        = TRIM( EosTableName ) )

  CALL ComputeAuxiliary_Fluid_TABLE &
        ( D2, E2, Ne2, P, T, Y, S, Em, Gm, Cs_table )
  
  CALL ComputeSpecificInternalEnergy_TABLE &
        ( D2, T, Y, Em, dEdD, dEdT, dEdY )

  IF ( Units ) THEN
    dEdD = dEdD / ( (Erg * Centimeter**3) / Gram**2 )
    dEdT = dEdT / ( Erg / Gram / Kelvin )
    dEdY = dEdY / ( Erg / Gram )
  END IF

  CALL FinalizeEquationOfState_TABLE

END SUBROUTINE ComputeDerivatives_InternalEnergy

SUBROUTINE ComputeDerivatives_InternalEnergy_Scalar & 
  ( D, E, Ne, dEdD, dEdT, dEdY, Units)

  USE KindModule, ONLY: &
    DP
  USE EquationOfStateModule_TABLE, ONLY: &
    ComputeSpecificInternalEnergy_TABLE, &
    ComputeAuxiliary_Fluid_TABLE
  USE EquationOfStateModule_TABLE, only: &
    InitializeEquationOfState_TABLE, &
    FinalizeEquationOfState_TABLE
  USE UnitsModule, ONLY: &
    Gram, &
    Centimeter, &
    Kelvin, &
    Erg
  
  REAL(DP), INTENT(in)  :: D, E, Ne
  REAL(DP) :: D2, E2, Ne2 ! remove units
  
  LOGICAL,  INTENT(in)  :: Units ! if True, incoming values are in physical units

  REAL(DP) :: P, Cs_table
  REAL(DP) :: T, Y, Em, Gm, S

  REAL(DP) :: dPdD, dPdT, dPdY
  REAL(DP), INTENT(out) :: dEdD, dEdT, dEdY

  CHARACTER(128) :: EosTableName 

  EosTableName = 'wl-EOS-SFHo-25-50-100.h5'

  ! ========================================================
  
  IF ( Units ) THEN
    D2 = D * ( Gram / Centimeter**3 )
    E2 = E * ( Erg / Centimeter**3 )
    Ne2 = Ne * ( 1.0_DP / Centimeter**3 )
  ELSE
    D2 = D 
    E2 = E
    Ne2 = Ne 
  END IF
  
  CALL InitializeEquationOfState_TABLE &
    ( EquationOfStateTableName_Option &
        = TRIM( EosTableName ) )

  CALL ComputeAuxiliary_Fluid_TABLE &
        ( D2, E2, Ne2, P, T, Y, S, Em, Gm, Cs_table )
  
  CALL ComputeSpecificInternalEnergy_TABLE &
        ( D2, T, Y, Em, dEdD, dEdT, dEdY )

  IF ( Units ) THEN
    dEdD = dEdD / ( (Erg * Centimeter**3) / Gram**2 )
    dEdT = dEdT / ( Erg / Gram / Kelvin )
    dEdY = dEdY / ( Erg / Gram )
  END IF

  CALL FinalizeEquationOfState_TABLE

END SUBROUTINE ComputeDerivatives_InternalEnergy_Scalar