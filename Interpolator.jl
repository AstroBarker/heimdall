"""
Set of table interpolation routines for Heimdall.
"""

using HDF5;

include("Matrix.jl")
include("Utils.jl")

# """
# Read the h5 EoS table. Currently only implemented for pressure, but extension is simple.
# """
# function read_table( fn::String )

#     fid :: HDF5File = h5open(fn,"r")

#     iRho :: Int32 = fid["ThermoState/iRho"][1]
#     iT   :: Int32 = fid["ThermoState/iT"][1]
#     iYe  :: Int32 = fid["ThermoState/iYe"][1]

#     dimensions  :: Array{Int32, 1} = fid["ThermoState/Dimensions"][:]

#     density     :: Array{Float64,1} = fid["ThermoState/Density"][:]
#     temperature :: Array{Float64,1} = fid["ThermoState/Temperature"][:]
#     ye          :: Array{Float64,1} = fid["ThermoState/Electron Fraction"][:]

#     dims_table :: Array{Int32, 1} = fid["DependentVariables/Dimensions"][:]
#     iPressure  :: Int32 = fid["DependentVariables/iPressure"][1]

#     pressure :: Array{Float64, 3} = fid["DependentVariables/Pressure"][:,:,:]

#     close( fid )

#     return pressure

# end

function read_table( fn::String )
    fid :: HDF5File = h5open(fn,"r")
    density     :: Array{Float64,1} = fid["ThermoState/Density"][:]
    temperature :: Array{Float64,1} = fid["ThermoState/Temperature"][:]
    ye          :: Array{Float64,1} = fid["ThermoState/Electron Fraction"][:]

    return fid, density, temperature, ye
end


function close_table( table::HDF5File )
    close( table )
end


function load_pressure( table::HDF5File )
    return table["DependentVariables/Pressure"][:,:,:]
end


function TriLinear( p000::Float64, p100::Float64, p010::Float64, p110::Float64, 
                    p001::Float64, p101::Float64, p011::Float64, p111::Float64, 
                    dX1::Float64,  dX2::Float64,  dX3::Float64 )

    ddX1 :: Float64 = 1.0 - dX1
    ddX2 :: Float64 = 1.0 - dX2
    ddX3 :: Float64 = 1.0 - dX3

    val  :: Float64 = ddX3 * ( ddX2 * ( ddX1 * p000 + dX1 * p100 )   
                     + dX2 * ( ddX1 * p010 + dX1 * p110 ) )
                     + dX3 * ( ddX2 * ( ddX1 * p001 + dX1 * p101 )   
                     + dX2 * ( ddX1 * p011 + dX1 * p111 ) )

    # xd = (x - x0)/(x1 - x0);
    # yd = (y - y0)/(y1 - y0);
    # zd = (z - z0)/(z1 - z0);

    xd :: Float64 = dX1
    yd :: Float64 = dX2
    zd :: Float64 = dX3

    # Interpolate along x
    c00 :: Float64 = p000 * (1 - xd) + p100 * xd;
    c01 :: Float64 = p001 * (1 - xd) + p101 * xd;
    c10 :: Float64 = p010 * (1 - xd) + p110 * xd;
    c11 :: Float64 = p011 * (1 - xd) + p111 * xd;

    #interpolate y
    c0  :: Float64 = c00 * (1 - yd) + c10 * yd;
    c1  :: Float64 = c01 * (1 - yd) + c11 * yd;

    # interpolate z
    val :: Float64 = c0 * (1 - zd) + c1 * zd;
    
    return val;
end


function tricubic_get_coeff_stacked( x::Array{Float64,1} )

    a :: Array{Float64, 1} = zeros( 64 )
    for i in 1:64
        for j in 1:64
            a[i]+=A[i,j] * x[j];
        end
    end

  return a
end

function tricubic_get_coeff( f::Array{Float64,1}, dfdx::Array{Float64,1}, dfdy::Array{Float64,1}, 
                             dfdz::Array{Float64,1}, d2fdxdy::Array{Float64,1}, d2fdxdz::Array{Float64,1}, 
                             d2fdydz::Array{Float64,1}, d3fdxdydz::Array{Float64,1} )

  x :: Array{Float64, 1} = zeros( 64 )
    for i in 1 : 8
        x[0+i]=f[i];
        x[8+i]=dfdx[i];
        x[16+i]=dfdy[i];
        x[24+i]=dfdz[i];
        x[32+i]=d2fdxdy[i];
        x[40+i]=d2fdxdz[i];
        x[48+i]=d2fdydz[i];
        x[56+i]=d3fdxdydz[i];
    end
  tricubic_get_coeff_stacked( x );
end


function TriCubic( a::Array{Float64, 2}, D::Float64, T::Float64, Y::Float64 )

    Interpolant :: Float64 = 0.0

    for i in 1 : 3
        for j in 1 : 3
            for k in 1 : 3
                Interpolant += a[ijk2n(i,j,k)] * x^i * y^j * z^k
            end
        end
    end

    return Interpolant
end


function LogInterpolate_Linear( D::Float64, T::Float64, Y::Float64, 
                                Ds::Array{Float64,1}, Ts::Array{Float64,1}, 
                                Ys::Array{Float64,1}, OS::Float64, 
                                Table::Array{Float64,3} )

    SizeDs :: Int64 = length( Ds )
    SizeTs :: Int64 = length( Ts )
    SizeYs :: Int64 = length( Ys )

    iD :: Int64 = Index1D( D, Ds )
    iT :: Int64 = Index1D( T, Ts )
    iY :: Int64 = Index1D( Y, Ys )

    dD :: Float64 = log10( D / Ds[iD] ) / log10( Ds[iD+1] / Ds[iD] )
    dT :: Float64 = log10( T / Ts[iT] ) / log10( Ts[iT+1] / Ts[iT] )
    dY :: Float64 = ( Y - Ys[iY] ) / ( Ys[iY+1] - Ys[iY] )

    p000 :: Float64 = Table[ iD  , iT  , iY   ]
    p100 :: Float64 = Table[ iD+1, iT  , iY   ]
    p010 :: Float64 = Table[ iD  , iT+1, iY   ]
    p110 :: Float64 = Table[ iD+1, iT+1, iY   ]
    p001 :: Float64 = Table[ iD  , iT  , iY+1 ]
    p101 :: Float64 = Table[ iD+1, iT  , iY+1 ]
    p011 :: Float64 = Table[ iD  , iT+1, iY+1 ]
    p111 :: Float64 = Table[ iD+1, iT+1, iY+1 ]

    Interpolant :: Float64 = 10.0^( TriLinear( p000, p100, p010, p110, 
              p001, p101, p011, p111, dD, dT, dY ) ) - OS

    return Interpolant
end