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

function TriLinear()

    return 0.0
end

function TriCubic()
    return 0.0
end