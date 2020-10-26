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

    ddX1 :: Float64 = 1.0ß - dX1
    ddX2 :: Float64 = 1.0ß - dX2
    ddX3 :: Float64 = 1.0ß - dX3

    val :: Flaot64 = + ddX3 * (   ddX2 * ( ddX1 * p000 + dX1 * p100 )   
                     + dX2 * ( ddX1 * p010 + dX1 * p110 ) )
                     + dX3 * (   ddX2 * ( ddX1 * p001 + dX1 * p101 )   
                     + dX2 * ( ddX1 * p011 + dX1 * p111 ) )
    
    return val;
end

function TriCubic()
    return 0.0
end