#!/usr/bin/env julia

# ===============================================================================
# __    __   _______  __  .___  ___.  _______       ___       __       __      
# |  |  |  | |   ____||  | |   \/   | |       \     /   \     |  |     |  |     
# |  |__|  | |  |__   |  | |  \  /  | |  .--.  |   /  ^  \    |  |     |  |     
# |   __   | |   __|  |  | |  |\/|  | |  |  |  |  /  /_\  \   |  |     |  |     
# |  |  |  | |  |____ |  | |  |  |  | |  '--'  | /  _____  \  |  `----.|  `----.
# |__|  |__| |_______||__| |__|  |__| |_______/ /__/     \__\ |_______||_______|
# ===============================================================================

"""
Set of functions for calling the EoS routines in thornado in postprocessing in thornado output.
Primarily we're interested in obtaining the thermodynamic derivatives concistently with the EoS.

Prerequisties:
--------------

None :-)

Implemented functions:
----------------------

ComputeDerivatives_Pressure(D, E, Ne) - 2 methods
    Compute derivatives of pressure from conserved variables D, E, Ne.
    Can accept all vectors D, E, Ne or all scalar D, E, Ne

ComputeDerivatives_SpecificInternalEnergy(D, E, Ne) - 2 methods
    Compute derivatives of specific internal energy from conserved variables D, E, Ne.
    Can accept all vectors D, E, Ne or all scalar D, E, Ne

Compute_R1(D, E, Ne, Vu1, Vu2, Vu3, Ye, Em, Cs, Gmdd11, Gmdd22, Gmdd33) - 2 methods
    Compute matrix of right eigenvectors. Single point or vector of points 

Compute_invR1(D, E, Ne, Vu1, Vu2, Vu3, Ye, Em, Cs, Gmdd11, Gmdd22, Gmdd33) - 2 methods
    Compute inverse matrix of right eigenvectors. Single point or vector of points 

Compute_Characteristics(U, Vu1, Vu2, Vu3, Y, Em, Cs, Gmdd11, Gmdd22, Gmdd334) - 2 methods`
    Compute chaqracteristic fields w = R^-1 U   

RemoveUnits_U(U, u) - 2 methods
    Convert vector of conserved quantities to code units

AddUnits_U(U, u) - 2 methods
    Convert vector of conserved quantities to physical units

Compute_Cs(D, P, Y, dPdE, dPdDe, dPdTau ) - 2 methods
    Compute analytic form of sound speed from Barker et al. senior thesis.
    Accepts all vector quantities or all scalar quantities.

Quantities:
-----------
Units::Struct - struct to hold units for converting to/from code units/physical units
    u::Units implements this for some commonly used unit conversions.

TODO:
------

Implement flux jacobians
Needs testing & work for general coordinates.
Create methods for 2D data.

"""
module Heimdall

# ================================ Units Struct ================================

"""
To switch to/from code units and physical units.
Multiply physical units by corresponding units to move to code units. Divide to go back.
e.g., rho = 1e13 * ( G / Cm^3 )

Erg : Erg
G : Gram
Cm : Centimeter
Km : Kilometer
S : Second
Amu : AtomicMassUnit
"""
struct Units
    Erg :: Float64
    G   :: Float64
    Cm  :: Float64
    Km  :: Float64
    S   :: Float64
    Amu :: Float64
end

u = Units(8.2611082525066313e-052,
    7.4247138240457958E-031,
    1.0e-2,
    1000.0,
    299792458.0, 
    1.2329024849255997e-54)

# Makes u available after importing Heimdall 
export u

# ================================ Module Functions ================================

"""
Convert code units to physical units for vector of conserved quantities U.

Parameters:
-----------
U::Array{Float64,2} - array holding conserved quantities:
    U = hcat(D, S1, S2, S3, E, D * Ye)
"""
function AddUnits_U( U::Array{Float64,2} )

    nx = length(U[:,1])

    U_Units::Array{Float64,1} = 
        [ 7.424713824045795e-25, 
        2.4766179488230473e-35, 
        2.4766179488230473e-35, 
        2.4766179488230473e-35, 
        8.261108252506631e-48,
        7.424713824045795e-25 ]

    new_U::Array{Float64,2} = zero( U )

    for i in 1:nx
        @inbounds new_U[i,:] = U[i,:] ./ U_Units
    end

    return new_U

end


function AddUnits_U( U::Float64 )

    U_Units::Array{Float64,1} = 
        [ 7.424713824045795e-25, 
        2.4766179488230473e-35, 
        2.4766179488230473e-35, 
        2.4766179488230473e-35, 
        8.261108252506631e-48,
        7.424713824045795e-25 ]

    return U ./ U_Units

end


"""
Convert physical units to code units for vector of conserved quantities U.

Parameters:
-----------
U::Array{Float64,2} - array holding conserved quantities:
    U = hcat(D, S1, S2, S3, E, D * Ye)
"""
function RemoveUnits_U( U::Array{Float64,2} )

    nx = length(U[:,1])

    U_Units::Array{Float64,1} = 
        [ 7.424713824045795e-25, 
        2.4766179488230473e-35, 
        2.4766179488230473e-35, 
        2.4766179488230473e-35, 
        8.261108252506631e-48,
        7.424713824045795e-25 ]

    new_U::Array{Float64,2} = zero( U )

    for i in 1:nx
        @inbounds new_U[i,:] = U[i,:] .* U_Units
    end

    return new_U

end


function RemoveUnits_U( U::Float64 )

    U_Units::Array{Float64,1} = 
        [ 7.424713824045795e-25, 
        2.4766179488230473e-35, 
        2.4766179488230473e-35, 
        2.4766179488230473e-35, 
        8.261108252506631e-48,
        7.424713824045795e-25 ]

    return U .* U_Units

end


"""
Call ComputeDerivatives_Pressure() from EoS_jl.f90 to compute thermodynamic derivatives 
of pressure. This routine takes the conserved variables D, E, Ne as primary inputs and 
constructs the rest of the thermodynamic variables consistently from them.

Derivatives loaded into array - dPdD = array[idD]
idD = 1
idT = 2
idY = 3
idE = 4
idDe = 5
idTau = 6

Parameters:
-----------

D::Array{Float64, 1} - thornado density profile
E::Array{Float64, 1} - thornado conserved energy density profile
Ne::Array{Float64,1} - thornado conserved electron number density profile
Units_Option::Bool (default: true) - if true, inputs are given in physical units
"""
function ComputeDerivatives_Pressure( D::Array{Float64, 1}, E::Array{Float64, 1}, Ne::Array{Float64,1};
    Units_Option::Bool=true )

    # Initialize arrays to hold derivatives
    nx     :: Int32             = length( D )
    dPdD   :: Array{Float64, 1} = zeros( nx );
    dPdT   :: Array{Float64, 1} = zeros( nx );
    dPdY   :: Array{Float64, 1} = zeros( nx );
    dPdE   :: Array{Float64, 1} = zeros( nx );
    dPdDe  :: Array{Float64, 1} = zeros( nx ); 
    dPdTau :: Array{Float64, 1} = zeros( nx ); 

    # =============================================================
    # This calls the FORTRAN function :computederivatives_pressure_
    # FORTRAN compilation mangles the name. Cvoid is the return 
    # type, the next parameters are input types, followed 
    # by the arguements. 
    # =============================================================
    ccall( (:computederivatives_pressure_, "./EoS_jl.so"), Cvoid, 
    ( Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Int32}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, 
    Ref{Float64}, Ref{Int32} ), 
    D, E, Ne, nx, dPdD, dPdT, dPdY, dPdE, dPdDe, dPdTau, Units_Option )    

    return hcat( dPdD, dPdT, dPdY, dPdE, dPdDe, dPdTau )

end    


function ComputeDerivatives_Pressure( D::Float64, E::Float64, Ne::Float64; 
    Units_Option::Bool=true )

    # Initialize arrays to hold derivatives
    dPdD   :: Ref{Float64}   = 0.0;
    dPdT   :: Ref{Float64}   = 0.0;
    dPdY   :: Ref{Float64}   = 0.0;
    dPdE   :: Ref{Float64}   = 0.0;
    dPdDe  :: Ref{Float64}   = 0.0; 
    dPdTau :: Ref{Float64}   = 0.0; 

    # =============================================================
    # This calls the FORTRAN function :computederivatives_pressure_
    # FORTRAN compilation mangles the name. Cvoid is the return 
    # type, the next parameters are input types, followed 
    # by the arguements. 
    #
    # Slightly different syntax from the vector version above.
    # Probably due to my ignorance, but it's the only way I could 
    # get it to work.
    # =============================================================
    ccall( (:computederivatives_pressure_scalar_, "./EoS_jl.so"), Nothing, 
    ( Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, 
    Ref{Float64}, Ref{Int32} ), 
    D, E, Ne, dPdD, dPdT, dPdY, dPdE, dPdDe, dPdTau, Units_Option )

    return hcat( dPdD, dPdT, dPdY, dPdE, dPdDe, dPdTau )

end   


"""
Call ComputeDerivatives_InternalEnergy() from EoS_jl.f90 to compute thermodynamic derivatives 
of specific internal energy. This routine takes the conserved variables D, E, Ne as primary inputs and 
constructs the rest of the thermodynamic variables consistently from them.

Derivatives loaded into array - dEdD = array[idD]
idD = 1
idT = 2
idY = 3

Parameters:
-----------

D::Array{Float64, 1} - thornado density profile
E::Array{Float64, 1} - thornado conserved energy density profile
Ne::Array{Float64,1} - thornado conserved electron number density profile
Units_Option::Bool (default: true) - if true, inputs are given in physical units
"""
function ComputeDerivatives_SpecificInternalEnergy( D::Array{Float64, 1}, E::Array{Float64, 1}, Ne::Array{Float64,1};
    Units_Option::Bool=true )

    # Initialize arrays to hold derivatives
    nx     :: Int32             = length( D )
    dEdD   :: Array{Float64, 1} = zeros( nx );
    dEdT   :: Array{Float64, 1} = zeros( nx );
    dEdY   :: Array{Float64, 1} = zeros( nx );

    # ===================================================================
    # This calls the FORTRAN function :computederivatives_internalenergy_
    # FORTRAN compilation mangles the name. Cvoid is the return 
    # type, the next parameters are input types, followed 
    # by the arguements. 
    # ===================================================================
    ccall( (:computederivatives_internalenergy_, "./EoS_jl.so"), Cvoid, 
    ( Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Int32}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Int32} ), 
    D, E, Ne, nx, dEdD, dEdT, dEdY, Units_Option )

    return hcat( dEdD, dEdT, dEdY )

end  


function ComputeDerivatives_SpecificInternalEnergy( D::Float64, E::Float64, Ne::Float64; 
    Units_Option::Bool=true )

    # Initialize arrays to hold derivatives
    dEdD   :: Ref{Float64}   = 0.0;
    dEdT   :: Ref{Float64}   = 0.0;
    dEdY   :: Ref{Float64}   = 0.0;

    # ===================================================================
    # This calls the FORTRAN function :computederivatives_internalenergy_
    # FORTRAN compilation mangles the name. Cvoid is the return 
    # type, the next parameters are input types, followed 
    # by the arguements. 
    #
    # Slightly different syntax from the vector version above.
    # Probably due to my ignorance, but it's the only way I could 
    # get it to work.
    # ===================================================================
    ccall( (:computederivatives_internalenergy_scalar_, "./EoS_jl.so"), Nothing, 
    ( Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Int32} ), 
    D, E, Ne, dEdD, dEdT, dEdY, Units_Option )

    return hcat( dEdD, dEdT, dEdY )

end   


""" DEPRECIATED
Call ComputeDerivatives_Entropy() from EoS_jl.f90 to compute thermodynamic derivatives 
of entropy. This routine takes variables D, T, Y as primary inputs and 
constructs the rest of the thermodynamic variables consistently from them.

Derivatives loaded into array - dSdD = array[idD]
idD = 1
idT = 2
idY = 3

Parameters:
-----------

D::Array{Float64, 1} - thornado density profile
T::Array{Float64, 1} - thornado temperature profile
Y::Array{Float64,1} - thornado conserved electron fraction profile
Units_Option::Bool (default: true) - if true, inputs are given in physical units
"""
# function ComputeDerivatives_Entropy( D::Array{Float64, 1}, T::Array{Float64, 1}, Y::Array{Float64,1};
#     Units_Option::Bool=true )

#     # Initialize arrays to hold derivatives
#     nx     :: Int32             = length( D )
#     dSdD   :: Array{Float64, 1} = zeros( nx );
#     dSdT   :: Array{Float64, 1} = zeros( nx );
#     dSdY   :: Array{Float64, 1} = zeros( nx );

#     # ===================================================================
#     # This calls the FORTRAN function :computederivatives_entropy_
#     # FORTRAN compilation mangles the name. Cvoid is the return 
#     # type, the next parameters are input types, followed 
#     # by the arguements. 
#     # ===================================================================
#     ccall( (:computederivatives_entropy_, "./EoS_jl.so"), Cvoid, 
#     ( Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Int32}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Int32} ), 
#     D, T, Y, nx, dSdD, dSdT, dSdY, Units_Option )

#     return hcat( dSdD, dSdT, dSdY )

# end  


""" DEPRECIATED
Call ComputeDerivatives_Entropy() from EoS_jl.f90 to compute thermodynamic derivatives 
of entropy. This routine takes variables D, T, Y as primary inputs and 
constructs the rest of the thermodynamic variables consistently from them.

Derivatives loaded into array - dSdD = array[idD]
idD = 1
idT = 2
idY = 3

Parameters:
-----------

D::Array{Float64, 1} - thornado density profile
T::Array{Float64, 1} - thornado temperature profile
Y::Array{Float64,1} - thornado conserved electron fraction profile
Units_Option::Bool (default: true) - if true, inputs are given in physical units
"""
# function ComputeDerivatives_Entropy( D::Float64, T::Float64, Y::Float64;
#     Units_Option::Bool=true )

#     # Initialize arrays to hold derivatives
#     nx     :: Int32             = length( D )
#     dSdD   :: Float64 = zeros( nx );
#     dSdT   :: Float64 = zeros( nx );
#     dSdY   :: Float64 = zeros( nx );

#     # ===================================================================
#     # This calls the FORTRAN function :computederivatives_entropy_
#     # FORTRAN compilation mangles the name. Cvoid is the return 
#     # type, the next parameters are input types, followed 
#     # by the arguements. 
#     #
#     # Slightly different syntax from the vector version above.
#     # Probably due to my ignorance, but it's the only way I could 
#     # get it to work.
#     # ===================================================================
#     ccall( (:computederivatives_entropy_scalar_, "./EoS_jl.so"), Nothing, 
#     ( Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Int32} ), 
#     D, T, Y, dSdD, dSdT, dSdY, Units_Option )

#     return hcat( dSdD, dSdT, dSdY )

# end  


"""
Compute analytic sound speed using expression in Barker et al senior thesis.
If passed with units, this will be in cm/s.

Parameters:
-----------
D::Array{Float64, 1} - thornado density profile
P::Array{Float64, 1} - thornado pressure profile
Y::Array{Float64, 1} - thornado electron fraction profile
"""
function Compute_Cs( D::Array{Float64,1}, P::Array{Float64,1}, Y::Array{Float64,1}, 
    dPdE::Array{Float64,1}, dPdDe::Array{Float64,1}, dPdTau::Array{Float64,1} )

    Tau::Array{Float64,1} = 1.0 ./ D

    CsSq::Array{Float64,1} = Tau.^2 .* ( P .* dPdE .- dPdTau ) .+ Y .* dPdDe 
end


function Compute_Cs( D::Float64, P::Float64, Y::Float64, 
    dPdE::Float64, dPdDe::Float64, dPdTau::Float64 )

    Tau::Float64 = 1.0 / D
    CsSq::Float64 = Tau^2 * ( P * dPdE - dPdTau ) + Y * dPdDe 
end

"""
Compute analytic sound speed using standard expression.
If passed with units, this will be in cm/s.

Parameters:
-----------
dPdD::Array{Float64, 1} - derivative of pressure w.r.t density
dSdT::Array{Float64, 1} - derivative of entropy w.r.t temperature
dPdT::Array{Float64, 1} - derivative of pressure w.r.t temperature
dSdD::Array{Float64, 1} - derivative of entropy w.r.t density
"""
function Compute_Cs_2( dPdD::Array{Float64,1}, dSdT::Array{Float64,1}, dPdT::Array{Float64,1}, 
    dSdD::Array{Float64,1} )

    CsSq::Array{Float64,1} = dPdD - dPdT .* dSdD ./ dSdT
end


"""
Compute the matrix of right eigenvectors from Barker et al.

Parameters:
-----------
D::Array{Float64, 1} - thornado density profile
E::Array{Float64, 1} - thornado conserved energy density profile
Ne::Array{Float64,1} - thornado conserved electron fraction profile
Vu1::Array{Float64,1} - Velocity profile
Vu2::Array{Float64,1} - Velocity profile
Vu3::Array{Float64,1} - Velocity profile
Y::Array{Float64,1} - thornado electron fraction
Cs::Array{Float64,1} - thornado sound speed
Em::Array{Float64,1} - thornado speciic internal energy
Gmdd11::Array{Float64,1} - Diagonal components of metric
Gmdd22::Array{Float64,1} - Diagonal components of metric
Gmdd33::Array{Float64,1} - Diagonal components of metric
Units_Option::Bool (default: false) - if true, inputs are given in physical units
"""
function Compute_R1( D::Array{Float64, 1}, E::Array{Float64, 1}, Ne::Array{Float64, 1}, 
    Vu1::Array{Float64, 1}, Vu2::Array{Float64, 1}, Vu3::Array{Float64, 1}, 
    Y::Array{Float64,1}, Em::Array{Float64, 1}, Cs::Array{Float64,1},
    Gmdd11::Array{Float64,1}, Gmdd22::Array{Float64,1}, Gmdd33::Array{Float64,1}; 
    Units_Option::Bool=false )

    R::Array{Float64,3} = zeros( length(D), 6, 6 )

    if ( Units_Option == false )
        D *= (7.4247138240457958E-031 / (1.0e-2)^3.0 )
        E *= (8.2611082525066313e-52 / (1.0e-2)^3.0 )
        Ne /=  (1.0e-2)^3.0
        Vu1 *= (1.0e3 / 299792458.0)
        Vu2 *= (1.0e3 / 299792458.0)
        Vu3 *= (1.0e3 / 299792458.0)
        Cs *= (1.0e3 / 299792458.0)
        Em *= (8.2611082525066313e-52 / 7.4247138240457958E-031)
    else
        Vu1 *= 1e5
        Vu2 *= 1e5
        Vu3 *= 1e5
        Cs *= 1e5
    end

    dd::Array{Float64,2} = ComputeDerivatives_Pressure( D, E, Ne, Units_Option=Units_Option )

    Vd1::Array{Float64,1} = Gmdd11 .* Vu1
    Vd2::Array{Float64,1} = Gmdd22 .* Vu2
    Vd3::Array{Float64,1} = Gmdd33 .* Vu3

    Vsq::Array{Float64,1} = Vu1 .* Vd1 .+ Vu2 .* Vd2 .+ Vu3 .* Vd3

    Tau::Array{Float64,1} = 1.0 ./ D
    Delta::Array{Float64,1} = Vu1 .* Vd1 .- Vu2 .* Vd2 .- Vu3 .* Vd3
    B::Array{Float64,1} = 0.5 .* ( Delta .+ 2.0 .* Em .+ 
        (2.0 .* dd[:,6] .* Tau) ./ dd[:,4])
    X::Array{Float64,1} = (dd[:,4] .* ( Delta .+ 2.0 .* Em) .+ 2.0 .* dd[:,6] .* Tau )

    K::Array{Float64,1} = ( ( .- ( Y ./ Tau ) .* dd[5] .+ dd[:,4] .* ( 
          0.5 .* Vsq .+ Em ) .+ dd[:,6] .* Tau ) ./ ( dd[:,4] ) )
    H::Array{Float64,1} = ( Cs.^2 ./ ( dd[:,4] .* Tau ) ) .+ K
    # TODO: Replace H with Tau(E+P)
    # TODO: Try analytic sound speed?

    for i in 1:length(D)
        @inbounds R[i,:,1] = [ 1.0, Vd1[i] .- Cs[i] .* sqrt.( Gmdd11[i] ), Vd2[i], 
            Vd3[i], H[i] .- Cs[i] .* sqrt.( Gmdd11[i] ) .* Vu1[i], Y[i] ]
        @inbounds R[i,:,2] = [ 0.0, 0.0, 1.0, 0.0, Vu2[i], 0.0 ]
        @inbounds R[i,:,3] = [ 1.0, Vd1[i], 0.0, 0.0, B[i], 0.0 ]
        @inbounds R[i,:,4] = [ 1.0, Vd1[i], 0.0, 0.0, 0.0,
            (Tau[i] .* X[i]) ./ (2.0 * dd[:,5][i]) ]
        @inbounds R[i,:,5] = [ 0.0, 0.0, 0.0, 1.0, Vu3[i], 0.0 ]
        @inbounds R[i,:,6] = [ 1.0, Vd1[i] + Cs[i] .* sqrt.( Gmdd11[i] ), Vd2[i],
            Vd3[i], H[i] .+ Cs[i] .* sqrt.( Gmdd11[i] ) .* Vu1[i], Y[i] ]

    end

    return R

end


function Compute_R1( D::Float64, E::Float64, Ne::Float64, 
    Vu1::Float64, Vu2::Float64, Vu3::Float64, 
    Y::Float64, Em::Float64, Cs::Float64,
    Gmdd11::Float64, Gmdd22::Float64, Gmdd33::Float64; Units_Option::Bool=false )

    R::Array{Float64,2} = zeros(6,6)

    if ( Units_Option == false )
        D *= (7.4247138240457958E-031 / (1.0e-2)^3.0 )
        E *= (8.2611082525066313e-52 / (1.0e-2)^3.0 )
        Ne /=  (1.0e-2)^3.0
        Vu1 *= (1.0e3 / 299792458.0)
        Vu2 *= (1.0e3 / 299792458.0)
        Vu3 *= (1.0e3 / 299792458.0)
        Cs *= (1.0e3 / 299792458.0)
        Em *= (8.2611082525066313e-52 / 7.4247138240457958E-031)
    else
        Vu1 *= 1e5
        Vu2 *= 1e5
        Vu3 *= 1e5
        Cs *= 1e5
    end

    dd::Array{Float64,2} = ComputeDerivatives_Pressure( D, E, Ne, Units_Option=Units_Option )
    println(dd[:,1])

    Vd1::Float64 = Gmdd11 .* Vu1
    Vd2::Float64 = Gmdd22 .* Vu2
    Vd3::Float64 = Gmdd33 .* Vu3

    Vsq::Float64 = Vu1 * Vd1 .+ Vu2 * Vd2 .+ Vu3 * Vd3

    Tau::Float64 = 1.0 ./ D
    Delta::Float64 = Vu1 * Vd1 .- Vu2 * Vd2 .- Vu3 * Vd3
    B::Float64 = 0.5 .* ( Delta .+ 2.0 .* Em .+ 
        (2.0 .* dd[:,6][1] * Tau) ./ dd[:,4][1])
    X::Float64 = (dd[:,4][1] .* ( Delta .+ 2.0 * Em) .+ 2.0 .* dd[:,6][1] .* Tau )

    K::Float64 = ( ( .- ( Y ./ Tau ) .* dd[:,5][1] .+ dd[:,4][1] .* ( 
          0.5 .* Vsq .+ Em ) .+ dd[:,6][1] .* Tau ) ./ ( dd[:,4][1] ) )
    H::Float64 = ( Cs.^2 ./ ( dd[:,4][1] .* Tau ) ) .+ K

    # TODO: Replace H with Tau(E+P)
    # TODO: Try analytic sound speed?

    R[:,1] = [ 1.0, Vd1 .- Cs .* sqrt.( Gmdd11 ), Vd2, 
       Vd3, H .- Cs .* sqrt.( Gmdd11 ) .* Vu1, Y ]
    R[:,2] = [ 0.0, 0.0, 1.0, 0.0, Vu2, 0.0 ]
    R[:,3] = [ 1.0, Vd1, 0.0, 0.0, B, 0.0 ]
    R[:,4] = [ 1.0, Vd1, 0.0, 0.0, 0.0,
       (Tau .* X) ./ (2.0 * dd[:,5][1]) ]
    R[:,5] = [ 0.0, 0.0, 0.0, 1.0, Vu3, 0.0 ]
    R[:,6] = [ 1.0, Vd1 .+ Cs .* sqrt.( Gmdd11 ), Vd2,
        Vd3, H .+ Cs .* sqrt.( Gmdd11 ) .* Vu1, Y ]

    return R
end


"""
Compute the inverse matrix of right eigenvectors 1 from Barker et al.

Parameters:
-----------
D::Array{Float64,1} - thornado density profile
E::Array{Float64,1} - thornado conserved energy density profile
Ne::Array{Float64,1} - thornado conserved electron fraction profile
Vu1::Array{Float64,1} - Velocity profile
Vu2::Array{Float64,1} - Velocity profile
Vu3::Array{Float64,1} - Velocity profile
Y::Array{Float64,1} - thornado electron fraction
Cs::Array{Float64,1} - thornado sound speed
Em::Array{Float64,1} - thornado speciic internal energy
Gmdd11::Array{Float64,1} - Diagonal components of metric
Gmdd22::Array{Float64,1} - Diagonal components of metric
Gmdd33::Array{Float64,1} - Diagonal components of metric
Units_Option::Bool (default: false) - if true, inputs are given in physical units
"""
function Compute_invR1( D::Array{Float64,1}, E::Array{Float64,1}, Ne::Array{Float64,1}, 
    Vu1::Array{Float64,1}, Vu2::Array{Float64,1}, Vu3::Array{Float64,1}, 
    Y::Array{Float64,1}, Em::Array{Float64,1}, Cs::Array{Float64,1},
    Gmdd11::Array{Float64,1}, Gmdd22::Array{Float64,1}, Gmdd33::Array{Float64,1}; 
    Units_Option::Bool=false )

    invR::Array{Float64,3} = zeros( length(D), 6, 6 )
    
    if ( Units_Option == false )
        D *= (7.4247138240457958E-031 / (1.0e-2)^3.0 )
        E *= (8.2611082525066313e-52 / (1.0e-2)^3.0 )
        Ne /=  (1.0e-2)^3.0
        Vu1 *= (1.0e3 / 299792458.0)
        Vu2 *= (1.0e3 / 299792458.0)
        Vu3 *= (1.0e3 / 299792458.0)
        Cs *= (1.0e3 / 299792458.0)
        Em *= (8.2611082525066313e-52 / 7.4247138240457958E-031)
    else
        Vu1 *= 1e5
        Vu2 *= 1e5
        Vu3 *= 1e5
        Cs *= 1e5
    end

    dd::Array{Float64,2} = ComputeDerivatives_Pressure( D, E, Ne, Units_Option=Units_Option )

    Vd1::Array{Float64,1} = Gmdd11 .* Vu1
    Vd2::Array{Float64,1} = Gmdd22 .* Vu2
    Vd3::Array{Float64,1} = Gmdd33 .* Vu3

    Tau::Array{Float64,1} = 1.0 ./ D

    Phi_u1::Array{Float64,1} = dd[:,4] .* Tau .* Vu1
    Phi_u2::Array{Float64,1} = dd[:,4] .* Tau .* Vu2
    Phi_u3::Array{Float64,1} = dd[:,4] .* Tau .* Vu3

    Phi_d1::Array{Float64,1} = dd[:,4] .* Tau .* Vd1
    Phi_d2::Array{Float64,1} = dd[:,4] .* Tau .* Vd2
    Phi_d3::Array{Float64,1} = dd[:,4] .* Tau .* Vd3

    Vsq::Array{Float64,1} = Vu1 .* Vd1 .+ Vu2 .* Vd2 .+ Vu3 .* Vd3

    Delta::Array{Float64,1} = Vu1 .* Vd1 .- Vu2 .* Vd2 .- Vu3 .* Vd3
    B::Array{Float64,1} = 0.5 .* ( Delta .+ 2.0 .* Em .+ 
        (2.0 .* dd[:,6] .* Tau) ./ dd[:,4])
    X::Array{Float64,1} = (dd[:,4] .* ( Delta .+ 2.0 .* Em) .+ 2.0 .* dd[:,6] .* Tau )

    K::Array{Float64,1} = ( ( .- ( Y ./ Tau ) .* dd[:,5] .+ dd[:,4] .* ( 
          0.5 .* Vsq .+ Em ) .+ dd[:,6] .* Tau ) ./ ( dd[:,4] ) )
    H::Array{Float64,1} = ( Cs.^2 ./ ( dd[:,4] .* Tau ) ) .+ K
    Alpha::Array{Float64,1} = 2.0 .* Y .* dd[:,5] .- X .* Tau
    W = Tau .* ( dd[:,4] .* ( Vsq .- 2.0 * Em )
               .- 2.0 .* dd[:,6] .* Tau )
    invCsSq = 1.0 ./ ( Cs.^2 )

    # TODO: Replace H with Tau(E+P)
    # TODO: Try analytic sound speed?

    # TODO: These are not array operations. Could this be rewritten to be vectorized?

    for i in 1:length(D)
        @inbounds invR[i,:,1] = invCsSq[i] .*
            [ + 0.25 .* (W[i] .+ 2.0 .* Cs[i] .* sqrt.( Gmdd11[i] ) .* Vu1[i]), 
            - 0.5 .* Vd2[i] .* W[i],
            + (2.0 * Cs[i].^2 .* X[i] .+ Alpha[i] .* W[i] ./ Tau[i]) ./ (2.0 .* X[i]),
            - (Y[i]) .* dd[:,5][i] .* W[i] ./ (X[i] .* Tau[i]),
            - 0.5 .* Vd3[i] .* W[i],
            + 0.25 .* (W[i] .- 2.0 .* Cs[i] .* sqrt.( Gmdd11[i] ) .* Vu1[i]) ]

        @inbounds invR[i,:,2] = invCsSq[i] .*
            [ - 0.5 .* ( ( Cs[i] ./ sqrt.( Gmdd11[i] ) ) + Phi_u1[i] ),
            + Phi_u1[i] .* Vd2[i],
            - Phi_u1[i] .* Alpha[i] ./ (X[i] .* Tau[i]),
            + 2.0 * Y[i] .* dd[:,5][i] .* Phi_u1[i] ./ (X[i] .* Tau[i]),
            + Phi_u1[i] .* Vd3[i],
            + 0.5 .* ( ( Cs[i] ./ sqrt.( Gmdd11[i] ) ) - Phi_u1[i] ) ]

        @inbounds invR[i,:,3] = invCsSq[i] .* 
            [ - 0.5 * Phi_u2[i],
            + Cs[i].^2 + Phi_u2[i] .* Vd2[i],
            - Phi_u2[i] .* Alpha[i] ./ (X[i] .* Tau[i]),
            + 2.0 .* Y[i] .* dd[:,5][i] .* Phi_u2[i] ./ (X[i] .* Tau[i]),
            + Phi_u2[i] .* Vd3[i],
            - 0.5 .* Phi_u2[i] ]

        @inbounds invR[i,:,4] = invCsSq[i] .*
            [ - 0.5 * Phi_u3[i],
            + Phi_u3[i] .* Vd2[i],
            - Phi_u3[i] .* Alpha[i] ./ (X[i] .* Tau[i]),
            + 2.0 .* Y[i] .* dd[:,5][i] .* Phi_u3[i] ./ (X[i] .* Tau[i]),
            + Cs[i].^2 .+ Phi_u3[i] .* Vd3[i],
            - 0.5 * Phi_u3[i] ]

        @inbounds invR[i,:,5] = invCsSq[i] .*
            [ + 0.5 * dd[:,4][i] .* Tau[i],
            - Phi_d2[i],
            + dd[:,4][i] .* Alpha[i]  ./ X[i],
            - ((2.0 * Y[i] .* dd[:,5][i] .* dd[:,4][i]) ./ X[i]),
            - Phi_d3[i],
            + 0.5 * dd[:,4][i] .* Tau[i] ]

        @inbounds invR[i,:,6] = invCsSq[i] .*
            [ + 0.5 .* dd[:,5][i],
            - Vd2[i] .* dd[:,5][i],
            + (dd[:,5][i] .* (-2.0 * Cs[i].^2 + Alpha[i])) ./ (Tau[i] .* X[i]),
            + 2.0 .* dd[:,5][i] .* (Cs[i].^2 .- Y[i] .* dd[:,5][i]) ./ (Tau[i] .* X[i]),
            - Vd3[i] .* dd[:,5][i],
            + 0.5 .* dd[:,5][i] ]

    end

    return invR
end


function Compute_invR1( D::Float64, E::Float64, Ne::Float64, 
    Vu1::Float64, Vu2::Float64, Vu3::Float64, 
    Y::Float64, Em::Float64, Cs::Float64,
    Gmdd11::Float64, Gmdd22::Float64, Gmdd33::Float64; Units_Option::Bool=true )

    invR::Array{Float64,2} = zeros(6,6)

    if ( Units_Option == false )
        D *= (7.4247138240457958E-031 / (1.0e-2)^3.0 )
        E *= (8.2611082525066313e-52 / (1.0e-2)^3.0 )
        Ne /=  (1.0e-2)^3.0
        Vu1 *= (1.0e3 / 299792458.0)
        Vu2 *= (1.0e3 / 299792458.0)
        Vu3 *= (1.0e3 / 299792458.0)
        Cs *= (1.0e3 / 299792458.0)
        Em *= (8.2611082525066313e-52 / 7.4247138240457958E-031)
    else
        Vu1 *= 1e5
        Vu2 *= 1e5
        Vu3 *= 1e5
        Cs *= 1e5
    end

    dd::Array{Float64,2} = ComputeDerivatives_Pressure( D, E, Ne, Units_Option=Units_Option )

    Vd1::Float64 = Gmdd11 .* Vu1
    Vd2::Float64 = Gmdd22 .* Vu2
    Vd3::Float64 = Gmdd33 .* Vu3

    Tau::Float64 = 1.0 ./ D

    Phi_u1::Float64 = dd[:,4][1] .* Tau .* Vu1
    Phi_u2::Float64 = dd[:,4][1] .* Tau .* Vu2
    Phi_u3::Float64 = dd[:,4][1] .* Tau .* Vu3

    Phi_d1::Float64 = dd[:,4][1] .* Tau .* Vd1
    Phi_d2::Float64 = dd[:,4][1] .* Tau .* Vd2
    Phi_d3::Float64 = dd[:,4][1] .* Tau .* Vd3

    Vsq::Float64 = Vu1 .* Vd1 + Vu2 .* Vd2 + Vu3 .* Vd3

    Delta::Float64 = Vu1 .* Vd1 .- Vu2 .* Vd2 .- Vu3 .* Vd3
    B::Float64 = 0.5 .* ( Delta + 2.0 .* Em + 
        (2.0 .* dd[:,6][1] .* Tau) ./ dd[:,4][1])
    X::Float64 = (dd[:,4][1] .* ( Delta .+ 2.0 .* Em) .+ 2.0 .* dd[:,6][1] .* Tau )

    K::Float64 = ( ( - ( Y ./ Tau ) .* dd[:,5][1] .+ dd[:,4][1] .* ( 
          0.5 .* Vsq .+ Em ) .+ dd[:,6][1] .* Tau ) ./ ( dd[:,4][1] ) )
    H::Float64 = ( Cs.^2 ./ ( dd[:,4][1] .* Tau ) ) .+ K
    Alpha::Float64 = 2.0 .* Y .* dd[:,5][1] .- X .* Tau
    W = Tau .* ( dd[:,4][1] .* ( Vsq .- 2.0 .* Em )
               .- 2.0 .* dd[:,6][1] .* Tau )
    invCsSq = 1.0 ./ ( Cs.^2 )

    # TODO: Replace H with Tau(E+P)
    # TODO: Try analytic sound speed?

    invR[:,1] = invCsSq .*
        [ + 0.25 .* (W .+ 2.0 .* Cs .* sqrt.( Gmdd11 ) .* Vu1), 
          - 0.5 .* Vd2 .* W,
          + (2.0 .* Cs.^2 .* X .+ Alpha .* W ./ Tau) ./ (2.0 .* X),
          - (Y) .* dd[:,5][1] .* W / (X .* Tau),
          - 0.5 * Vd3 .* W,
          + 0.25 .* (W .- 2.0 .* Cs .* sqrt.( Gmdd11 ) .* Vu1) ]

    invR[:,2] = invCsSq .*
        [ - 0.5 .* ( ( Cs ./ sqrt.( Gmdd11 ) ) .+ Phi_u1 ),
          + Phi_u1 .* Vd2,
          - Phi_u1 .* Alpha ./ (X .* Tau),
          + 2.0 .* Y .* dd[:,5][1] .* Phi_u1 ./ (X .* Tau),
          + Phi_u1 .* Vd3,
          + 0.5 .* ( ( Cs ./ sqrt.( Gmdd11 ) ) .- Phi_u1 ) ]

    invR[:,3] = invCsSq .* 
        [ - 0.5 .* Phi_u2,
          + Cs.^2 .+ Phi_u2 .* Vd2,
          - Phi_u2 .* Alpha ./ (X .* Tau),
          + 2.0 .* Y .* dd[:,5][1] .* Phi_u2 ./ (X .* Tau),
          + Phi_u2 .* Vd3,
          - 0.5 .* Phi_u2 ]

    invR[:,4] = invCsSq .*
        [ - 0.5 * Phi_u3,
          + Phi_u3 .* Vd2,
          - Phi_u3 .* Alpha ./ (X .* Tau),
          + 2.0 .* Y .* dd[:,5][1] .* Phi_u3 ./ (X .* Tau),
          + Cs.^2 + Phi_u3 .* Vd3,
          - 0.5 * Phi_u3 ]

    invR[:,5] = invCsSq .*
        [ + 0.5 .* dd[:,4][1] .* Tau,
          - Phi_d2,
          + dd[:,4][1] .* Alpha  ./ X,
          - ((2.0 .* Y .* dd[:,5][1] .* dd[:,4][1]) ./ X),
          - Phi_d3,
          + 0.5 .* dd[:,4][1] .* Tau ]

    invR[:,6] = invCsSq .*
        [ + 0.5 .* dd[:,5][1],
          - Vd2 .* dd[:,5][1],
          + (dd[:,5][1] .* (-2.0 .* Cs.^2 .+ Alpha)) ./ (Tau .* X),
          + 2.0 .* dd[:,5][1] .* (Cs.^2 .- Y .* dd[:,5][1]) ./ (Tau .* X),
          - Vd3 .* dd[:,5][1],
          + 0.5 .* dd[:,5][1] ]

    return invR
end


"""
Compute the characteristic quantities w = invR U.
w[i,:] will access the i-th characteristic field for i in 1:length( domain ).
Calculations done in code units. If Units_Option = true, w is converted to physical units 
at the end.

Parameters:
-----------
U::Array{Float64,2} - array holding conserved quantities on the columns:
    U = hcat(uCF_D, uCF_S1, uCF_S2, uCF_S3, uCF_E, uCF_D * uAF_Ye)
Ne::Array{Float64,1} - Conserved electron density
Vu1::Array{Float64,1} - Velocity profile
Vu2::Array{Float64,1} - Velocity profile
Vu3::Array{Float64,1} - Velocity profile
Y::Array{Float64,1} - thornado electron fraction
Cs::Array{Float64,1} - thornado sound speed
Em::Array{Float64,1} - thornado speciic internal energy
Gmdd11::Array{Float64,1} - Diagonal components of metric
Gmdd22::Array{Float64,1} - Diagonal components of metric
Gmdd33::Array{Float64,1} - Diagonal components of metric
Units_Option::Bool (default: false) - if true, inputs are given in physical units
"""
function Compute_Characteristics( U::Array{Float64,2}, Ne::Array{Float64,1},
    Vu1::Array{Float64,1}, Vu2::Array{Float64,1}, Vu3::Array{Float64,1}, 
    Y::Array{Float64,1}, Em::Array{Float64,1}, Cs::Array{Float64,1},
    Gmdd11::Array{Float64,1}, Gmdd22::Array{Float64,1}, Gmdd33::Array{Float64,1}; 
    Units_Option::Bool=false)

    nx = length( U[:,1] )

    # Note that we do the calculations *without* units
    invR::Array{Float64,3} = Compute_invR1( U[:,1], U[:,5], Ne, Vu1, Vu2, Vu3, Y, Em, Cs, 
        Gmdd11, Gmdd22, Gmdd33, Units_Option=false )

    new_U::Array{Float64,2} = RemoveUnits_U( U )

    w::Array{Float64,2} = zero( U )

    for i in 1:nx
        @inbounds w[i,:] = invR[i,:,:] * new_U[i,:]
    end

    if Units_Option
        w = AddUnits_U( w )
    end

    return w

end


function Compute_Characteristics( U::Array{Float64,2}, 
    Vu1::Float64, Vu2::Float64, Vu3::Float64, 
    Y::Float64, Em::Float64, Cs::Float64,
    Gmdd11::Float64, Gmdd22::Float64, Gmdd33::Float64; 
    Units_Option::Bool=false)

    nx = length( U[:,1] )

    # Note that we do the calculations *without* units
    invR::Array{Float64,3} = Compute_invR1( U[:,1], U[:,5], Ne, Vu1, Vu2, Vu3, Y, Em, Cs, 
        Gmdd11, Gmdd22, Gmdd33, Units_Option=false )

    new_U::Array{Float64,1} = RemoveUnits_U( U )

    w::Array{Float64,2} = zero( U )

    for i in 1:nx
        @inbounds w[i,:] = invR[i,:,:] * new_U[i,:]
    end

    if Units_Option
        w = AddUnits_U( w )
    end

    return w

end

# End module
end