SUBROUTINE ComputeDerivatives_Pressure & 
  ( D, E, Ne, nx, Gmdd11, Gmdd22, Gmdd33, dPdD, dPdT, dPdY, dPdE, dPdDe, dPdTau)

  USE KindModule, ONLY: &
    DP
  USE Euler_UtilitiesModule_NonRelativistic, ONLY: &
    ComputePrimitive_Euler_NonRelativistic
  USE EquationOfStateModule_TABLE, ONLY: &
    ComputeSpecificInternalEnergy_TABLE, &
    ComputeAuxiliary_Fluid_TABLE, &
    ComputePressure_TABLE
  USE EquationOfStateModule, ONLY: &
    InitializeEquationOfState, &
    FinalizeEquationOfState
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

  INTEGER  :: i
  REAL(DP) :: Vu1, Vu2, Vu3, Vd1, Vd2, Vd3
  REAL(DP), DIMENSION(nx) :: P, Cs, Cs_table
  REAL(DP), DIMENSION(nx) :: Tau, T, Y, Vsq, CsSq, W, Em, Gm, S

  REAL(DP), INTENT(out), DIMENSION(nx) :: dPdD, dPdT, dPdY
  REAL(DP), DIMENSION(nx) :: dEdD, dEdT, dEdY
  REAL(DP), INTENT(out), DIMENSION(nx) :: dPdE, dPdDe, dPdTau

  CHARACTER(128) :: EosTableName 

  EosTableName = 'wl-EOS-SFHo-25-50-100.h5'

  ! ========================================================
  
  D2 = D * ( Gram / Centimeter**3 )
  E2 = E * ( Erg / Centimeter**3 )
  Ne2 = Ne * ( 1.0_DP / Centimeter**3 )
  
  CALL InitializeEquationOfState &
    ( EquationOfState_Option &
        = 'TABLE', &
      EquationOfStateTableName_Option &
        = TRIM( EosTableName ) )

  ! CALL ComputePrimitive_Euler_NonRelativistic &
  ! ( U(iCF_D ), U(iCF_S1), U(iCF_S2), &
  !   U(iCF_S3), U(iCF_E ), U(iCF_Ne), &
  !   D, Vu1, Vu2, Vu3, E, Ne, &
  !   G(iGF_Gm_dd_11), &
  !   G(iGF_Gm_dd_22), &
  !   G(iGF_Gm_dd_33) )

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

  dPdE = dPdE / ( Gram / Centimeter**3 )
  dPdDe = dPdDe / ( Erg / Gram )
  dPdTau = dPdTau / ( Erg * Gram / Centimeter**6 )

  dPdD = dPdD / ( Erg / Gram )
  dPdT = dPdT / ( ( Erg / Centimeter**3 ) / Kelvin )
  dPdY = dPdY / ( Erg / Centimeter**3 )  

  CALL FinalizeEquationOfState

END SUBROUTINE ComputeDerivatives_Pressure

SUBROUTINE ComputeSpecificInternalenergy_Output(D, E, Ne, nx, Em)
  ! Output specific internal energy
  USE KindModule, ONLY: &
    DP
  USE Euler_UtilitiesModule_NonRelativistic, ONLY: &
    ComputePrimitive_Euler_NonRelativistic
  USE EquationOfStateModule_TABLE, ONLY: &
    ComputeSpecificInternalEnergy_TABLE, &
    ComputeAuxiliary_Fluid_TABLE, &
    ComputePressure_TABLE
  USE EquationOfStateModule, ONLY: &
    InitializeEquationOfState, &
    FinalizeEquationOfState
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

  REAL(DP), DIMENSION(nx) :: P, Cs_table
  REAL(DP), DIMENSION(nx) :: T, Y, Gm, S
  REAL(DP), DIMENSION(nx), INTENT(out) :: Em

  REAL(DP), DIMENSION(nx) :: dEdD, dEdT, dEdY

  CHARACTER(128) :: EosTableName 

  EosTableName = 'wl-EOS-SFHo-25-50-100.h5'

  ! ========================================================

  D2 = D * ( Gram / Centimeter**3 )
  E2 = E * ( Erg / Centimeter**3 )
  Ne2 = Ne * ( 1.0_DP / Centimeter**3 )

  CALL InitializeEquationOfState &
    ( EquationOfState_Option &
        = 'TABLE', &
      EquationOfStateTableName_Option &
        = TRIM( EosTableName ) )

  CALL ComputeAuxiliary_Fluid_TABLE &
        ( D2, E2, Ne2, P, T, Y, S, Em, Gm, Cs_table )

  CALL ComputeSpecificInternalEnergy_TABLE &
        ( D2, T, Y, Em, dEdD, dEdT, dEdY )

  Em = Em / ( Erg / Gram )

  CALL FinalizeEquationOfState

END SUBROUTINE ComputeSpecificInternalenergy_Output

SUBROUTINE ComputeDerivatives_Pressure_Scalar & 
  ( D, E, Ne, nx, Gmdd11, Gmdd22, Gmdd33, dPdD, dPdT, dPdY, dPdE, dPdDe, dPdTau)

  USE KindModule, ONLY: &
    DP
  USE Euler_UtilitiesModule_NonRelativistic, ONLY: &
    ComputePrimitive_Euler_NonRelativistic
  USE EquationOfStateModule_TABLE, ONLY: &
    ComputeSpecificInternalEnergy_TABLE, &
    ComputeAuxiliary_Fluid_TABLE, &
    ComputePressure_TABLE
  USE EquationOfStateModule, ONLY: &
    InitializeEquationOfState, &
    FinalizeEquationOfState
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

  INTEGER  :: i
  REAL(DP) :: Vu1, Vu2, Vu3, Vd1, Vd2, Vd3
  REAL(DP) :: P, Cs, Cs_table
  REAL(DP) :: K, H, Tau, T, Y, Vsq, CsSq, W, Em, Gm, S

  REAL(DP), INTENT(out) :: dPdD, dPdT, dPdY
  REAL(DP) :: dEdD, dEdT, dEdY
  REAL(DP), INTENT(out) :: dPdE, dPdDe, dPdTau

  CHARACTER(128) :: EosTableName 

  EosTableName = 'wl-EOS-SFHo-25-50-100.h5'

  ! ========================================================

  D2 = D * ( Gram / Centimeter**3 )
  E2 = E * ( Erg / Centimeter**3 )
  Ne2 = Ne * ( 1.0_DP / Centimeter**3 )
  
  CALL InitializeEquationOfState &
    ( EquationOfState_Option &
        = 'TABLE', &
      EquationOfStateTableName_Option &
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

  dPdE = dPdE / ( Gram / Centimeter**3 )
  dPdDe = dPdDe / ( Erg / Gram )
  dPdTau = dPdTau / ( Erg * Gram / Centimeter**6 )

  dPdD = dPdD / ( Erg / Gram )
  dPdT = dPdT / ( ( Erg / Centimeter**3 ) / Kelvin )
  dPdY = dPdY / ( Erg / Centimeter**3 ) 

END SUBROUTINE ComputeDerivatives_Pressure_Scalar