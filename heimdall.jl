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

DataFrames.jl - (https://juliadata.github.io/DataFrames.jl/stable/man/getting_started/#Installation-1)

Implemented:
------------

computeDerivatives_Pressure(D, E, Ne, Gmdd11, Gmdd22, Gmdd33) - 2 methods
    Compute derivatives of pressure from conserved variables D, E, Ne.
    Can accept all vectors D, E, Ne or all scalar D, E, Ne

cs(D, P, Y, dPdE, dPdDe, dPdTau ) - 2 methods
    Compute analytic form of sound speed from Barker et al. senior thesis.
    Accepts all vector quantities or all scalar quantities.

TODO:
------

Implement characteristic decomposition functions.
"""

using DataFrames

function ComputeDerivatives_pressure( D::Array{Float64, 1}, E::Array{Float64, 1}, Ne::Array{Float64,1}, 
    Gmdd11::Float64, Gmdd22::Float64, Gmdd33::Float64 )
    """
    Call ComputeDerivatives_Pressure() from EoS_jl.f90 to compute thermodynamic derivatives 
    of pressure. This routine takes the conserved variables D, E, Ne as primary inputs and 
    constructs the rest of the thermodynamic variables consistently from them.

    Derivatives are loaded into a dataframe. Accessible as, e.g., df.dPdD

    TODO: Data structures may need thought for multi-D
    
    Parameters:
    -----------

    D::Array{Float64, 1} - thornado density profile
    E::Array{Float64, 1} - thornado conserved energy density profile
    Ne::Array{Float64,1} - thornado conserved electron fraction profile
    Gmdd11::Float64 - Diagonal components of metric
    Gmdd22::Float64 - Diagonal components of metric
    Gmdd33::Float64 - Diagonal components of metric
    """

    # Initialize arrays to hold derivatives
    nx     :: Int64             = length( D )
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
    ( Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Int64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, 
    Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64} ), 
    D, E, Ne, nx, Gmdd11, Gmdd22, Gmdd33, dPdD, dPdT, dPdY, dPdE, dPdDe, dPdTau )

    df = DataFrame(dPdD=dPdD, dPdT=dPdT, dPdY=dPdY, dPdE=dPdE, dPdDe=dPdDe, dPdTau=dPdTau)

    return df

end    

function ComputeDerivatives_pressure( D::Float64, E::Float64, Ne::Float64, 
    Gmdd11::Float64, Gmdd22::Float64, Gmdd33::Float64 )
    """
    Call ComputeDerivatives_Pressure() from EoS_jl.f90 to compute thermodynamic derivatives 
    of pressure. This routine takes the conserved variables D, E, Ne as primary inputs and 
    constructs the rest of the thermodynamic variables consistently from them.
    
    Parameters:
    -----------

    D::Array{Float64, 1} - thornado density profile
    E::Array{Float64, 1} - thornado conserved energy density profile
    Ne::Array{Float64,1} - thornado conserved electron fraction profile
    Gmdd11::Float64 - Diagonal components of metric
    Gmdd22::Float64 - Diagonal components of metric
    Gmdd33::Float64 - Diagonal components of metric
    """

    # Initialize arrays to hold derivatives
    nx     :: Ref{Int64}     = 0;
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
    # SLightly different syntax from the vector version above.
    # Probably due to my ignorance, but it's the only way I could 
    # get it to work.
    # =============================================================
    ccall( (:computederivatives_pressure_scalar_, "./EoS_jl.so"), Nothing, 
    ( Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Int64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, 
    Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64} ), 
    D, E, Ne, nx, Gmdd11, Gmdd22, Gmdd33, dPdD, dPdT, dPdY, dPdE, dPdDe, dPdTau )
    
    df = DataFrame(dPdD=dPdD.x[1][1], dPdT=dPdT.x[1][1], dPdY=dPdY.x[1][1], 
                   dPdE=dPdE.x[1][1], dPdDe=dPdDe.x[1][1], dPdTau=dPdTau.x[1][1])

    return df

end   

function ComputeSpecificInternalEnergy( D::Array{Float64,1}, E::Array{Float64,1}, Ne::Array{Float64,1} )
    """
    Call ComputeSpecificInternalEnergy_Output() from EoS_jl.f90.

    !!! NOT IN USE !!!

    Parameters:
    -----------
    D::Array{Float64, 1} - thornado density profile
    E::Array{Float64, 1} - thornado conserved energy density profile
    Ne::Array{Float64,1} - thornado conserved electron fraction profile
    """

    # Initialize arrays
    nx = length( D )
    Em = zeros( nx )

    # ======================================================================
    # This calls the FORTRAN function :computespecificinternalenergy_output_
    # FORTRAN compilation mangles the name. Cvoid is the return 
    # type, the next parameters are input types, followed 
    # by the arguements. 
    # ======================================================================
    ccall( (:computespecificinternalenergy_output_, "./EoS_jl.so"), Cvoid, 
    ( Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Int64}, Ref{Float64} ), 
    D, E, Ne, nx, Em )

    return Em

end



function Cs( D::Array{Float64,1}, P::Array{Float64,1}, Y::Array{Float64,1}, 
    dPdE::Array{Float64,1}, dPdDe::Array{Float64,1}, dPdTau::Array{Float64,1} )
    """
    Compute analytic sound speed in cm/s using expression in Barker et al senior thesis.

    Parameters:
    -----------
    D::Array{Float64, 1} - thornado density profile
    P::Array{Float64, 1} - thornado pressure profile
    Y::Array{Float64, 1} - thornado electron fraction profile
    """

    Tau::Array{Float64,1} = 1.0 ./ D

    CsSq::Array{Float64,1} = Tau.^2 .* ( P .* dPdE .- dPdTau ) .+ Y .* dPdDe 
end

function Cs( D::Float64, P::Float64, dPdE::Float64, 
    dPdTau::Float64, Y::Float64, dPdDe::Float64 )
    """
    Compute analytic sound speed in cm/s using expression in Barker et al senior thesis.
    
    Parameters:
    -----------
    D::Float64 - thornado density profile
    P::Float64 - thornado pressure profile
    Y::Float64 - thornado electron fraction profile
    """

    Tau::Float64 = 1.0 / D
    CsSq::Float64 = Tau^2 * ( P * dPdE - dPdTau ) + Y * dPdDe 
end

function Compute_R1( D::Array{Float64, 1}, E::Array{Float64, 1}, Ne::Array{Float64, 1}, 
    Vu1::Array{Float64, 1}, Vu2::Array{Float64, 1}, Vu3::Array{Float64, 1}, 
    Y::Array{Float64,1}, Em::Array{Float64, 1}, Cs::Array{Float64,1},
    Gmdd11::Float64, Gmdd22::Float64, Gmdd33::Float64 )
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
    Gmdd11::Float64 - Diagonal components of metric
    Gmdd22::Float64 - Diagonal components of metric
    Gmdd33::Float64 - Diagonal components of metric
    """

    R::Array{Float64,3} = zeros( length(D), 6, 6 )

    dd::DataFrame = ComputeDerivatives_pressure( D, E, Ne, Gmdd11, Gmdd22, Gmdd33 );

    Vd1::Array{Float64,1} = Gmdd11 * Vu1
    Vd2::Array{Float64,1} = Gmdd22 * Vu2
    Vd3::Array{Float64,1} = Gmdd33 * Vu3

    Vsq::Array{Float64,1} = Vu1 .* Vd1 + Vu2 .* Vd2 + Vu3 .* Vd3

    Tau::Array{Float64,1} = 1.0 ./ D
    Delta::Array{Float64,1} = Vu1 .* Vd1 - Vu2 .* Vd2 - Vu3 .* Vd3
    B::Array{Float64,1} = 0.5 .* ( Delta + 2.0 .* Em + 
        (2.0 .* dd.dPdTau .* Tau) ./ dd.dPdE)
    X::Array{Float64,1} = (dd.dPdE .* ( Delta + 2.0 * Em) + 2.0 * dd.dPdTau .* Tau )

    K::Array{Float64,1} = ( ( - ( Y ./ Tau ) .* dd.dPdDe + dd.dPdE .* ( 
          0.5 * Vsq + Em ) + dd.dPdTau .* Tau ) ./ ( dd.dPdE ) )
    H::Array{Float64,1} = ( Cs.^2 ./ ( dd.dPdE .* Tau ) ) + K
    # TODO: Replace H with Tau(E+P)
    # TODO: Try analytic sound speed?

    # TODO: WHAT to do about R? What shape in this useage? R[nx, 6, 6]

    for i in 1:length(D)

    R[i,:,1] = [ 1.0, Vd1[i] - Cs[i] .* sqrt.( Gmdd11 ), Vd2[i], 
       Vd3[i], H[i] - Cs[i] .* sqrt.( Gmdd11 ) .* Vu1[i], Y[i] ]
    R[i,:,2] = [ 0.0, 0.0, 1.0, 0.0, Vu2[i], 0.0 ]
    R[i,:,3] = [ 1.0, Vd1[i], 0.0, 0.0, B[i], 0.0 ]
    R[i,:,4] = [ 1.0, Vd1[i], 0.0, 0.0, 0.0,
       (Tau[i] .* X[i]) ./ (2.0 * dd.dPdDe[i]) ]
    R[i,:,5] = [ 0.0, 0.0, 0.0, 1.0, Vu3[i], 0.0 ]
    R[i,:,6] = [ 1.0, Vd1[i] + Cs[i] .* sqrt.( Gmdd11 ), Vd2[i],
        Vd3[i], H[i] + Cs[i] .* sqrt.( Gmdd11 ) .* Vu1[i], Y[i] ]

    end

    return R

end

function Compute_R1( D::Float64, E::Float64, Ne::Float64, 
    Vu1::Float64, Vu2::Float64, Vu3::Float64, 
    Y::Float64, Em::Float64, Cs::Float64,
    Gmdd11::Float64, Gmdd22::Float64, Gmdd33::Float64 )
    """
    Compute the matrix of right eigenvectors from Barker et al.

    Note: we access derivatives as dd.dPdx[1] as they are returned as single element arrays.

    Parameters:
    -----------
    D::Float64 - thornado density profile
    E::Float64 - thornado conserved energy density profile
    Ne::Float64 - thornado conserved electron fraction profile
    Vu1::Float64 - Velocity profile
    Vu2::Float64 - Velocity profile
    Vu3::Float64 - Velocity profile
    Y::Float64 - thornado electron fraction
    Cs::Float64 - thornado sound speed
    Em::Float64 - thornado speciic internal energy
    Gmdd11::Float64 - Diagonal components of metric
    Gmdd22::Float64 - Diagonal components of metric
    Gmdd33::Float64 - Diagonal components of metric
    """

    R::Array{Float64,2} = zeros(6,6)

    dd::DataFrame = ComputeDerivatives_pressure( D, E, Ne, Gmdd11, Gmdd22, Gmdd33 );

    Vd1::Float64 = Gmdd11 * Vu1
    Vd2::Float64 = Gmdd22 * Vu2
    Vd3::Float64 = Gmdd33 * Vu3

    Vsq::Float64 = Vu1 * Vd1 + Vu2 * Vd2 + Vu3 * Vd3

    Tau::Float64 = 1.0 ./ D
    Delta::Float64 = Vu1 * Vd1 - Vu2 * Vd2 - Vu3 * Vd3
    B::Float64 = 0.5 .* ( Delta + 2.0 .* Em + 
        (2.0 .* dd.dPdTau[1] * Tau) ./ dd.dPdE[1])
    X::Float64 = (dd.dPdE[1] .* ( Delta + 2.0 * Em) + 2.0 * dd.dPdTau[1] .* Tau )

    K::Float64 = ( ( - ( Y ./ Tau ) .* dd.dPdDe[1] + dd.dPdE[1] .* ( 
          0.5 * Vsq + Em ) + dd.dPdTau[1] .* Tau ) ./ ( dd.dPdE[1] ) )
    H::Float64 = ( Cs.^2 ./ ( dd.dPdE[1] .* Tau ) ) + K
    # TODO: Replace H with Tau(E+P)
    # TODO: Try analytic sound speed?

    R[:,1] = [ 1.0, Vd1 - Cs .* sqrt.( Gmdd11 ), Vd2, 
       Vd3, H - Cs .* sqrt.( Gmdd11 ) .* Vu1, Y ]
    R[:,2] = [ 0.0, 0.0, 1.0, 0.0, Vu2, 0.0 ]
    R[:,3] = [ 1.0, Vd1, 0.0, 0.0, B, 0.0 ]
    R[:,4] = [ 1.0, Vd1, 0.0, 0.0, 0.0,
       (Tau .* X) ./ (2.0 * dd.dPdDe[1]) ]
    R[:,5] = [ 0.0, 0.0, 0.0, 1.0, Vu3, 0.0 ]
    R[:,6] = [ 1.0, Vd1 + Cs .* sqrt.( Gmdd11 ), Vd2,
        Vd3, H + Cs .* sqrt.( Gmdd11 ) .* Vu1, Y ]

    return R
end

function Compute_invR1( D::Float64, E::Float64, Ne::Float64, 
    Vu1::Float64, Vu2::Float64, Vu3::Float64, 
    Y::Float64, Em::Float64, Cs::Float64,
    Gmdd11::Float64, Gmdd22::Float64, Gmdd33::Float64 )
    """
    Compute the inverse matrix of right eigenvectors 1 from Barker et al.

    Note: we access derivatives as dd.dPdx[1] as they are returned as single element arrays.

    Parameters:
    -----------
    D::Float64 - thornado density profile
    E::Float64 - thornado conserved energy density profile
    Ne::Float64 - thornado conserved electron fraction profile
    Vu1::Float64 - Velocity profile
    Vu2::Float64 - Velocity profile
    Vu3::Float64 - Velocity profile
    Y::Float64 - thornado electron fraction
    Cs::Float64 - thornado sound speed
    Em::Float64 - thornado speciic internal energy
    Gmdd11::Float64 - Diagonal components of metric
    Gmdd22::Float64 - Diagonal components of metric
    Gmdd33::Float64 - Diagonal components of metric
    """

    invR::Array{Float64,2} = zeros(6,6)

    dd::DataFrame = ComputeDerivatives_pressure( D, E, Ne, Gmdd11, Gmdd22, Gmdd33 );

    Vd1::Float64 = Gmdd11 * Vu1
    Vd2::Float64 = Gmdd22 * Vu2
    Vd3::Float64 = Gmdd33 * Vu3

    Tau::Float64 = 1.0 ./ D

    Phi_u1::Float64 = dd.dPdE[1] .* Tau .* Vu1
    Phi_u2::Float64 = dd.dPdE[1] .* Tau .* Vu2
    Phi_u3::Float64 = dd.dPdE[1] .* Tau .* Vu3

    Phi_d1::Float64 = dd.dPdE[1] .* Tau .* Vd1
    Phi_d2::Float64 = dd.dPdE[1] .* Tau .* Vd2
    Phi_d3::Float64 = dd.dPdE[1] .* Tau .* Vd3

    Vsq::Float64 = Vu1 * Vd1 + Vu2 * Vd2 + Vu3 * Vd3

    Delta::Float64 = Vu1 * Vd1 - Vu2 * Vd2 - Vu3 * Vd3
    B::Float64 = 0.5 .* ( Delta + 2.0 .* Em + 
        (2.0 .* dd.dPdTau[1] * Tau) ./ dd.dPdE[1])
    X::Float64 = (dd.dPdE[1] .* ( Delta + 2.0 * Em) + 2.0 * dd.dPdTau[1] .* Tau )

    K::Float64 = ( ( - ( Y ./ Tau ) .* dd.dPdDe[1] + dd.dPdE[1] .* ( 
          0.5 * Vsq + Em ) + dd.dPdTau[1] .* Tau ) ./ ( dd.dPdE[1] ) )
    H::Float64 = ( Cs.^2 ./ ( dd.dPdE[1] .* Tau ) ) + K
    Alpha::Float64 = 2.0 * Y .* dd.dPdDe[1] - X .* Tau
    W = Tau .* ( dd.dPdE[1] .* ( Vsq - 2.0 * Em )
               - 2.0 .* dd.dPdTau[1] .* Tau )
    invCsSq = 1.0 ./ ( Cs.^2 )

    # TODO: Replace H with Tau(E+P)
    # TODO: Try analytic sound speed?

    invR[:,1] = invCsSq .*
        [ + 0.25 * (W + 2.0 * Cs .* sqrt.( Gmdd11 ) .* Vu1), 
          - 0.5 * Vd2 .* W,
          + (2.0 * Cs.^2 * X + Alpha .* W ./ Tau)./(2.0 .* X),
          - (Y) .* dd.dPdDe[1] .* W / (X .* Tau),
          - 0.5 * Vd3 .* W,
          + 0.25 * (W - 2.0 * Cs .* sqrt.( Gmdd11 ) .* Vu1) ]

    invR[:,2] = invCsSq .*
        [ - 0.5 .* ( ( Cs ./ sqrt.( Gmdd11 ) ) + Phi_u1 ),
          + Phi_u1 .* Vd2,
          - Phi_u1 .* Alpha ./ (X .* Tau),
          + 2.0 * Y .* dd.dPdDe[1] .* Phi_u1 ./ (X .* Tau),
          + Phi_u1 .* Vd3,
          + 0.5 .* ( ( Cs ./ sqrt.( Gmdd11 ) ) - Phi_u1 ) ]

    invR[:,3] = invCsSq .* 
        [ - 0.5 * Phi_u2,
          + Cs.^2 + Phi_u2 .* Vd2,
          - Phi_u2 .* Alpha ./ (X .* Tau),
          + 2.0 .* Y .* dd.dPdDe[1] .* Phi_u2 ./ (X .* Tau),
          + Phi_u2 .* Vd3,
          - 0.5 * Phi_u2 ]

    invR[:,4] = invCsSq .*
        [ - 0.5 * Phi_u3,
          + Phi_u3 .* Vd2,
          - Phi_u3 .* Alpha ./ (X .* Tau),
          + 2.0 .* Y .* dd.dPdDe[1] .* Phi_u3 ./ (X .* Tau),
          + Cs.^2 + Phi_u3 .* Vd3,
          - 0.5 * Phi_u3 ]

    invR[:,5] = invCsSq .*
        [ + 0.5 * dd.dPdE[1] .* Tau,
          - Phi_d2,
          + dd.dPdE[1] .* Alpha  ./ X,
          - ((2.0 * Y .* dd.dPdDe[1] .* dd.dPdE[1]) ./ X),
          - Phi_d3,
          + 0.5 * dd.dPdE[1] .* Tau ]

    invR[:,6] = invCsSq .*
        [ + 0.5 * dd.dPdDe[1],
          - Vd2 .* dd.dPdDe[1],
          + (dd.dPdDe[1] .* (-2.0 * Cs.^2 + Alpha)) ./ (Tau .* X),
          + 2.0 * dd.dPdDe[1] .* (Cs.^2 .- Y .* dd.dPdDe[1]) ./ (Tau .* X),
          - Vd3 .* dd.dPdDe[1],
          + 0.5 * dd.dPdDe[1] ]

    return invR
end