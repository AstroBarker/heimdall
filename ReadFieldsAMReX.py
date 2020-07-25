#!/usr/local/bin/python3

import yt
import numpy as np
import subprocess
from os import listdir
from os import environ
from os.path import isfile
from sys import argv, exit

"""
Python function for loading thornado-AMReX data.
The output is such that a field may be accessed as Data[i]

Parameters:
-----------
DataDirectory: str
    Directory containing simulation outputs

Default usage, plots last plotfile in DataDirectory:

  $ python3 PlotFieldsAMReX.py

Alernate usage, plot specific file in DataDirectory:

  $ python3 PlotFieldsAMReX.py thornado_00000010

TODO:
  - Add SymLogNorm scaling
"""

def Load_AMReX( DataDirectory: str ):

    # https://yt-project.org/doc/faq/index.html#how-can-i-change-yt-s-log-level
    yt.funcs.mylog.setLevel(40) # Suppress yt warnings

    # --- Get user's THORNADO_DIR directory ---

    THORNADO_DIR = environ["THORNADO_DIR"]
    # THORNADO_DIR = THORNADO_DIR[:-1].decode( "utf-8" ) + '/'

    #### ========== User Input ==========

    # Load Vorticity? [Not Implemented]
    Vorticity = False
    # Polytropic Constant? [Not Implemented]
    PolytropicConstant = False

    # Specify plot file base name
    PlotFileBaseName = 'thornado'

    #### ====== End of User Input =======

    # Append "/" to DataDirectory, if not present
    if( not DataDirectory[-1] == '/' ): DataDirectory += '/'

    if( len( argv ) == 1 ):
        # Get last plotfile in directory

        FileArray \
        = np.sort( np.array( [ file for file in listdir( DataDirectory ) ] ) )

        FileList = []

        for iFile in range( FileArray.shape[0] ):

            sFile = FileArray[iFile]

            if( sFile[0:len(PlotFileBaseName)+1] == PlotFileBaseName + '_' \
                and sFile[len(PlotFileBaseName)+1].isdigit() ):
                FileList.append( sFile )

        FileArray = np.array( FileList )
        File      = FileArray[-1]

    elif( len( argv ) == 2 ):
        File = argv[1]

    else:
        print( 'Invalid number of optional parameters' )
        exit( 'Exiting...' )

    # Remove "/" at end of filename, if present
    if ( File[-1] == '/' ): File = File[:-1]

    ds = yt.load( '{:}'.format( DataDirectory + File ) )

    print( 'Reading from file: {:}'.format( File ) )
    MaxLevel = ds.index.max_level
    Time     = ds.current_time
    nX       = ds.domain_dimensions
    xL       = ds.domain_left_edge
    xH       = ds.domain_right_edge

    # Get dimensionality of problem
    if  ( nX[1] == 1 and nX[2] == 1 ):
        nDims = 1
    elif( nX[2] == 1 ):
        nDims = 2
    else:
        nDims = 3

    """
    https://yt-project.org/doc/reference/api/
    yt.data_objects.construction_data_containers.html#yt.data_objects.
    construction_data_containers.YTCoveringGrid
    """
    CoveringGrid \
    = ds.covering_grid \
        ( level           = MaxLevel, \
            left_edge       = xL, \
            dims            = nX * 2**MaxLevel, \
            num_ghost_zones = nX[0] )

    # XXX.to_ndarray() strips yt array of units

    DataUnit = ''

    # === Begin loading data ===

    #nx = 
    #ny = 
    numFields = 15 # UPDATE if you add/remove from Data !!! 
    temp = CoveringGrid['PF_D' ].to_ndarray()

    if nDims == 1:
        nx = len( temp[0] )
        ny = 1
        nz = 1
        Data = np.zeros( (numFields, nx, ny, nz) )
    if nDims == 2:
        nx = len( temp[0] )
        ny = len( temp[1] )
        nz = 1
        Data = np.zeros( (numFields, nx, ny, nz) )
    if nDims == 3:
        nx = len( temp[0] )
        ny = len( temp[1] )
        nz = len( temp[2])
        Data = np.zeros( (numFields, nx, ny, nz) )

    Data[0,:,:] = temp
    iPF_D = 1

    # TODO: Make "fields = [...]" and loop over fields.
    Data[1,:,:,:] = CoveringGrid['PF_V1'].to_ndarray()
    iPF_V1 = 2

    Data[2,:,:,:] = CoveringGrid['PF_V2'].to_ndarray()
    iPF_V2 = 3

    Data[3,:,:,:] = CoveringGrid['PF_V3'].to_ndarray()
    iPF_V3 = 4

    Data[4,:,:] = CoveringGrid['PF_E' ].to_ndarray()
    iPF_E = 5

    Data[5,:,:,:] = CoveringGrid['CF_D' ].to_ndarray()
    iCF_D = 6

    Data[6,:,:,:] = CoveringGrid['CF_S1'].to_ndarray()
    iCF_S1 = 7

    Data[7,:,:,:] = CoveringGrid['CF_S2'].to_ndarray()
    iCF_S2 = 8

    Data[8,:,:,:] = CoveringGrid['CF_S3'].to_ndarray()
    iCF_S3 = 9

    Data[9,:,:,:] = CoveringGrid['CF_E' ].to_ndarray()
    iCF_E = 10

    Data[10,:,:,:] = CoveringGrid['AF_P' ].to_ndarray()
    iAF_P = 11

    Data[11,:,:,:] = CoveringGrid['AF_Cs'].to_ndarray()
    iAF_Cs = 12

    Data[12,:,:,:] = CoveringGrid['DF_Sh_X1'].to_ndarray()
    iDF_Sh_X1 = 13

    Data[13,:,:,:] = CoveringGrid['DF_Sh_X2'].to_ndarray()
    iDF_Sh_X2 = 14

    Data[14,:,:,:] = CoveringGrid['DF_Sh_X3'].to_ndarray()
    iDF_Sh_X3 = 15

    return Data, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iCF_D, iCF_S1, iCF_S2, iCF_S3, \
        iCF_E, iAF_P, iAF_Cs, iDF_Sh_X1, iDF_Sh_X2, iDF_Sh_X3

    # --- Derived Fields --- [None Implemented]

    # elif( PolytropicConstant ):
    #     PF_D  = CoveringGrid['PF_D' ].to_ndarray()
    #     AF_P  = CoveringGrid['AF_P' ].to_ndarray()
    #     AF_Gm = CoveringGrid['AF_Gm'].to_ndarray()
    #     Data  = AF_P / PF_D**AF_Gm
    #     if( UsePhysicalUnits ):
    #         if( round( AF_Gm[0][0][0], 1 ) == 1.4 ):
    #             DataUnit = 'erg/cm**3/(g/cm**3)**({:})'.format( \
    #                         AF_Gm[0][0][0] )
    #         else:
    #             DataUnit = 'erg/cm**3/(g/cm**3)**({:}/3)'.format( \
    #                         int( 3 * AF_Gm[0][0][0] ) )

    # elif( Vorticity ):

    #     dX1 = ( xH[0].to_ndarray() - xL[0].to_ndarray() ) / nX[0]
    #     dX2 = ( xH[1].to_ndarray() - xL[1].to_ndarray() ) / nX[1]

    #     XL = xL.to_ndarray() + 0.5 * np.array( [ dX1, dX2, 0.0 ] )
    #     XH = xH.to_ndarray() - 0.5 * np.array( [ dX1, dX2, 0.0 ] )

    #     X1 = np.linspace( XL[0], XH[0], nX[0] )
    #     X2 = np.linspace( XL[1], XH[1], nX[1] )

    #     PF_V1 = CoveringGrid['PF_V1'].to_ndarray()
    #     PF_V2 = CoveringGrid['PF_V2'].to_ndarray()
    #     indX1 = np.linspace( 1, nX[0]-2, nX[0]-2, dtype = int )
    #     indX2 = np.linspace( 1, nX[1]-2, nX[1]-2, dtype = int )
    #     Data = np.zeros( (nX[0],nX[1],1), float )
    #     for j in indX1:
    #         for k in indX2:
    #             Data[j,k,0] \
    #             = 1.0 / X1[j] \
    #                 * ( ( X1[j+1]**2 * PF_V2[j+1,k] \
    #                         - X1[j-1]**2 * PF_V2[j-1,k] ) \
    #                     / ( 2.0 * dX1 ) \
    #                     - ( PF_V1[j,k+1] \
    #                         - PF_V1[j,k-1] ) \
    #                     / ( 2.0 * dX2 ) )

    #     # Apply boundary conditions to theta elements
    #     for j in indX1:

    #         # North pole
    #         k = 0
    #         Data[j,k,0] \
    #         = 1.0 / X1[j] \
    #             * ( ( X1[j+1]**2 * PF_V2[j+1,k] \
    #                     - X1[j-1]**2 * PF_V2[j-1,k] ) \
    #                 / ( 2.0 * dX1 ) \
    #                 - ( PF_V1[j,k+1] \
    #                     - PF_V1[j,k] ) \
    #                 / ( 2.0 * dX2 ) )

    #         # South pole
    #         k = nX[1]-1
    #         Data[j,k,0] \
    #         = 1.0 / X1[j] \
    #             * ( ( X1[j+1]**2 * PF_V2[j+1,k] \
    #                     - X1[j-1]**2 * PF_V2[j-1,k] ) \
    #                 / ( 2.0 * dX1 ) \
    #                 - ( PF_V1[j,k] \
    #                     - PF_V1[j,k-1] ) \
    #                 / ( 2.0 * dX2 ) )

    #     if( UsePhysicalUnits ):
    #         DataUnit = '1/s'

test, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iCF_D, iCF_S1, iCF_S2, iCF_S3, \
        iCF_E, iAF_P, iAF_Cs, iDF_Sh_X1, iDF_Sh_X2, iDF_Sh_X3 = \
        Load_AMReX( "/Users/barker/ornl_ut/codes/thornado/SandBox/AMReX/temp", ProblemName = "implosion", CoordinateSystem = "cartesian", UsePhysicalUnits = True )
print(test[0])