#!/usr/bin/env julia
"""
Read Thornado data and load into a JuliaDB
Simple, only implemented for 1D native thornado

Prerequisites:
--------------
HDF5.jl (https://github.com/JuliaIO/HDF5.jl)
DataFrames.jl (https://juliadata.github.io/DataFrames.jl/stable/man/getting_started/#Installation-1)
"""

using HDF5;
using DataFrames;

"""
Read 1D native thornado data.

Parameters:
-----------
Dir::String 
    Directory containing .h5 files
filenumber::String 
    String of output to read -- must contain leading zeros. e.g., 000165
run::String
    Simulation identifier: e.g., GravitationalCollapse, RiemannProblem
"""
function load_thornado_single( Dir::AbstractString, filenumber::AbstractString; run::AbstractString = "RiemannProblem" )

    # ================ Fluid Fields ================

    fn::String = string(Dir, "/", run, "_FluidFields_", filenumber, ".h5")
    
    fid = h5open(fn,"r")

    t      :: Float65          = fid["Time"][:][1]
    x1     :: Array{Float64,1} = fid["/Spatial Grid/X1"][:]   
    uCF_D  :: Array{Float64,1} = fid["/Fluid Fields/Conserved/Conserved Baryon Density"][:,1,1]
    uCF_S1 :: Array{Float64,1} = fid["/Fluid Fields/Conserved/Conserved Momentum Density (1)"][:,1,1]
    uCF_S2 :: Array{Float64,1} = fid["/Fluid Fields/Conserved/Conserved Momentum Density (2)"][:,1,1]
    uCF_S3 :: Array{Float64,1} = fid["/Fluid Fields/Conserved/Conserved Momentum Density (3)"][:,1,1]
    uPF_V1 :: Array{Float64,1} = fid["/Fluid Fields/Primitive/Three-Velocity (1)"][:,1,1]
    uPF_V2 :: Array{Float64,1} = fid["/Fluid Fields/Primitive/Three-Velocity (2)"][:,1,1]
    uPF_V3 :: Array{Float64,1} = fid["/Fluid Fields/Primitive/Three-Velocity (3)"][:,1,1]
    uPF_Ev :: Array{Float64,1} = fid["/Fluid Fields/Primitive/Internal Energy Density"][:,1,1]
    uAF_T  :: Array{Float64,1} = fid["/Fluid Fields/Auxiliary/Temperature"][:,1,1]
    uAF_Em :: Array{Float64,1} = fid["/Fluid Fields/Auxiliary/Specific Internal Energy"][:,1,1]
    uAF_Ye :: Array{Float64,1} = fid["/Fluid Fields/Auxiliary/Electron Fraction"][:,1,1]
    uAF_P  :: Array{Float64,1} = fid["/Fluid Fields/Auxiliary/Pressure"][:,1,1]
    uAF_S  :: Array{Float64,1} = fid["/Fluid Fields/Auxiliary/Entropy Per Baryon"][:,1,1]
    uCF_E  :: Array{Float64,1} = fid["/Fluid Fields/Conserved/Conserved Energy Density"][:,1,1]
    uCF_Ne :: Array{Float64,1} = fid["/Fluid Fields/Conserved/Conserved Electron Density"][:,1,1]
    uAF_Cs :: Array{Float64,1} = fid["/Fluid Fields/Auxiliary/Sound Speed"][:,1,1]

    # ================ Geometry Fields ================

    println("Time: $t" )

    return DataFrame( x1=x1, uCF_D=uCF_D, uCF_E=uCF_E, uCF_S1=uCF_S1, uCF_S2=uCF_S2, uCF_S3=uCF_S3,
        uAF_T=uAF_T, uAF_Ye=uAF_Ye, uAF_P=uAF_P, uAF_Cs=uAF_Cs, uCF_Ne=uCF_Ne, 
        uPF_V1=uPF_V1, uPF_V2=uPF_V2, uPF_V3=uPF_V3, uPF_Ev=uPF_Ev, uAF_Em=uAF_Em )
end

"""
Compute cell averages,

Parameters:
-----------
df::DataFrame 
    data, contained in a dataframe. Output of load_thornado_single()
nNodes::Int 
    number of nodes. Determines weights.
"""
function cell_average(df::DataFrame, nNodes::Int)

    if ( nNodes == 3 )
        wG :: Array{Float64,1} = [ 5.0/(1*18.0), 8.0/(1*18.0), 5.0/(1*18.0) ]
    elseif ( nNodes == 2)
        wG:: Array{Float64,1} = [0.5, 0.5]
    else
        println("Assuming nNodes = 1. Not implemented for > 3. Returning original DataFrame.")
        return df
    end

    N_N::Int64 = length( df[!,:x1] )

    N_K::Float64 = N_N / nNodes

    df2 :: DataFrame = DataFrame()
    k::Int64 = 1
    # Loop over columns of df
    for col in propertynames(df)
        avg::Array{Float64,1} = zeros( floor(Int, N_K) )
        for j in 1:floor(Int, N_K)
            avg[j] = sum( wG .* df[!,col][(j-1)*nNodes+1:(j)*nNodes] )
        end
        df2[!,col] = avg

    end
        
    return df2   


end