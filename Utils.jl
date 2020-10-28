"""
Utility functions
"""

"""
Convert an (i,j,k) index to 1D index n.
"""
function ijk2n( i::Int64, j::Int64, k::Int64 ) 
  return (i+4*j+16*k)
end

function point2xyz( p::Int64 )

    if (p == 0)
        return 0, 0, 0
    elseif (p == 1)
        return 1, 0, 0
    elseif (p == 2)
        return 0, 1, 0
    elseif (p == 3)
        return 1, 1, 0
    elseif (p == 4)
        return 0, 0, 1
    elseif (p == 5)
        return 1, 0, 1
    elseif (p == 6)
        return 0, 1, 1
    elseif (p == 7)
        return 1, 1, 1
    else
        return 0,0,0
    end

end

"""
Search algorithm.
"""
function Index1D( val::Float64, array::Array{Float64,1} )
   
    index :: Int64 = - 1 # error val
    size  :: Int64 = length( array )

    il :: Int64 = 0
    im :: Int64 = 0
    iu :: Int64 = size + 1
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