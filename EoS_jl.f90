SUBROUTINE ComputeDerivatives_Pressure & 
  ( D, E, Ne, nx, Gmdd11, Gmdd22, Gmdd33, dPdD, & 
    dPdT, dPdY, dPdE, dPdDe, dPdTau, Units)

  USE KindModule, ONLY: &
    DP
  USE Euler_UtilitiesModule_NonRelativistic, ONLY: &
    ComputePrimitive_Euler_NonRelativistic
  USE EquationOfStateModule_TABLE, ONLY: &
    ComputeSpecificInternalEnergy_TABLE, &
    ComputeAuxiliary_Fluid_TABLE, &
    ComputePressure_TABLE
  USE EquationOfStateModule_TABLE, only: &
    InitializeEquationOfState_TABLE, &
    FinalizeEquationOfState_TABLE
  USE ProgramInitializationModule, ONLY: &
    FinalizeProgram
  USE UnitsModule, ONLY: &
    AtomicMassUnit, &
    BoltzmannConstant, &
    Gram, &
    Centimeter, &
    Kelvin, &
    Dyne, &
    Erg, &
    MeV, & 
    Second, & 
    Kilometer
  
  INTEGER,  INTENT(in)  :: nx
  REAL(DP), INTENT(in)  :: D(nx), E(nx), Ne(nx)
  REAL(DP) :: D2(nx), E2(nx), Ne2(nx) ! remove units
  REAL(DP), INTENT(in) :: Gmdd11, Gmdd22, Gmdd33
  
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
  ( D, E, Ne, nx, Gmdd11, Gmdd22, Gmdd33, dPdD, & 
    dPdT, dPdY, dPdE, dPdDe, dPdTau, Units)

  USE KindModule, ONLY: &
    DP
  USE Euler_UtilitiesModule_NonRelativistic, ONLY: &
    ComputePrimitive_Euler_NonRelativistic
  USE EquationOfStateModule_TABLE, ONLY: &
    ComputeSpecificInternalEnergy_TABLE, &
    ComputeAuxiliary_Fluid_TABLE, &
    ComputePressure_TABLE
  USE EquationOfStateModule_TABLE, only: &
    InitializeEquationOfState_TABLE, &
    FinalizeEquationOfState_TABLE
  USE ProgramInitializationModule, ONLY: &
    FinalizeProgram
  USE UnitsModule, ONLY: &
    AtomicMassUnit, &
    BoltzmannConstant, &
    Gram, &
    Centimeter, &
    Kelvin, &
    Dyne, &
    Erg, &
    MeV

  INTEGER,  INTENT(in)  :: nx
  REAL(DP), INTENT(in) :: D, E, Ne
  REAL(DP) :: D2, E2, Ne2 ! remove units
  REAL(DP), INTENT(in) :: Gmdd11, Gmdd22, Gmdd33

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