module Functools

export lenscat_load

using DelimitedFiles
using FITSIO
#using ..Cosmology

"""
constructs the jackknife patches masks using Lambert cylindrical equal-area projection
"""
function construct_K(ra, dec; NK=100)
        
    K = zeros(Bool, (length(ra),NK+1))
    K[:, 1] .= true
    
    sqrt_NK = floor(Int64, sqrt(NK))
    cdec = sind.(dec)
    ramin, ramax = extrema(ra)
    decmin, decmax = extrema(dec)
    dra = (ramax - ramin)/sqrt_NK
    ddec = (decmax - decmin)/sqrt_NK
    
    c = 2
    for a in 1:sqrt_NK
        for d in 1:sqrt_NK
            mra = (ra .>= ramin + a*dra) .&& (ra .< ramin + (a+1)*dra)
            mdec = (cdec .>= decmin + d*ddec) .&& (cdec .< decmin + (d+1)*ddec)
            K[:, c] .= .~(mra.&&mdec)
            c += 1
        end
    end
    return K
end

"""
Loads a void lens catalogue given the range of radii, redshifts and void type and
constucts the jackknife masks. Optionally can be set to do MICE (w/ id as first column), return
only the Rv, RA, DEC, Z cols or fullshape, and separate the cat in NCHUNCKS.
"""
function lenscat_load(name::String, 
    Rv_min::Real, Rv_max::Real, z_min::Real, z_max::Real,
    delta_min::Real, delta_max::Real;
    rho1_min::Real=-1.0,rho1_max::Real=0.0, flag::Real=2.0, 
    NCHUNCKS::Integer=1, NK::Integer=100, octant::Bool=false, MICE::Bool=false, fullshape::Bool=true)
    
    if MICE
        RV, RA, DEC, Z, R1, R2 = 2,3,4,5,9,10
    else
        RV, RA, DEC, Z, R1, R2 = 1,2,3,4,8,9
    end

    L = readdlm(name)
    
    if octant
        println("not implemented")
        return nothing
    end
    
    mask = ((L[:,RV] .>= Rv_min).&&(L[:,RV] .< Rv_max).&&(L[:,Z] .>= z_min).&&(L[:,Z] .< z_max).&&(L[:,R1] .>= rho1_min)
    .&&(L[:,R1] .< rho1_max).&&(L[:,R2] .>= delta_min).&&(L[:,R2] .< delta_max).&&(L[:,12] .>= flag))
    
    nvoids = sum(mask)
    #K = K[1:nvoids,:]
    if fullshape
        L = L[mask, :]
    else
        L = L[mask, [RV, RA, DEC, Z]]
    end
    K = construct_K(L[:,RA], L[:, DEC]; NK=NK)
    
    if NCHUNCKS!=1
        if NCHUNCKS > nvoids
            NCHUNCKS = nvoids
        end
        lbins = ceil(Int64, nvoids/NCHUNCKS)
        slices = collect(NCHUNCKS:NCHUNCKS:(lbins-1)*NCHUNCKS)
        #slices = slices[(slices < nvoids)]
        L = [ L[i:min(i+NCHUNCKS-1, nvoids),:] for i in slices ]
        K = [ K[i:min(i+NCHUNCKS-1, nvoids),:] for i in slices ]
    end

    return L, K, nvoids

end

function sourcecat_load(name::String)
    f = FITS(name)
    cat = readdlm(f[1]; use_mmap=true)
    return cat
end

end #module