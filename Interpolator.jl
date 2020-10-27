"""
Set of table interpolation routines for Heimdall.
"""

using HDF5;

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
                    dX1::Float64, dX2::Float64, dX3::Float64 )

    ddX1 :: Float64 = 1.0 - dX1
    ddX2 :: Float64 = 1.0 - dX2
    ddX3 :: Float64 = 1.0 - dX3

    val :: Float64 =  ddX3 * ( ddX2 * ( ddX1 * p000 + dX1 * p100 )   
                     + dX2 * ( ddX1 * p010 + dX1 * p110 ) )
                     + dX3 * ( ddX2 * ( ddX1 * p001 + dX1 * p101 )   
                     + dX2 * ( ddX1 * p011 + dX1 * p111 ) )
    
    return val;
end

function TriCubic()
    return 0.0
end

function Index1D( val::Float64, array::Array{Float64,1} )
   
    index :: Int64 = - 1 # error val
    size :: Int64 = length( array )

    il :: Int64 = 0
    iu :: Int64 = size + 1
    im :: Int64 = 0
    while ( iu - il > 1 )
      im = floor(Int, (iu+il) / 2 ) # round down
      if ((array[size] > array[1]) & (val > array[im] ))
        il = im
      else
        iu = im
      end
    end

    if ( val == array[1] )
      index = 1
    elseif ( val == array[size] )
      index = size - 1
    else
      index = il
    end

    # only works for monatonically increasing array
    # index = maximum( findall( array .<= val ) )

    return index

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