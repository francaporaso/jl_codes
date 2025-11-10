using DataFrames
#using DelimitedFiles
using FITSIO
using NearestNeighbors

include("Cosmology.jl")
using .Cosmology

cosmo = Planck18

function cos_angular_separation(α1, α2, δ1, δ2)
    return sin(δ1)*sin(δ2) + cos(δ1)*cos(δ2)*cos(α1-α2)
end

function distance(χ1, χ2, α1, α2, δ1, δ2)
    return sqrt(χ1^2+χ2^2- 2*χ1*χ2*cos_angular_separation(α1, α2, δ1, δ2))
end

function sphere2box(ra::Real, dec::Real, χ::Real)
    x = χ*cos(ra)*cos(dec) 
    y = χ*sin(ra)*cos(dec) 
    z = χ*sin(dec) 
    return x,y,z
end

function sphere2box(ra::Vector{T}, dec::Vector{T}, χ::Vector{T}) where {T<:Real}
    x = @. χ*cos(ra)*cos(dec) 
    y = @. χ*sin(ra)*cos(dec) 
    z = @. χ*sin(dec) 
    return hcat(x,y,z)
end

gals = FITS("mockgals.fits")
gals = DataFrame(gals[2])

randgals = FITS("random_gals.fits")
randgals = DataFrame(randgals[2])

voids = FITS("mockvoids.fits")
voids = DataFrame(voids[2])

randvoids = FITS("random_voids.fits")
randvoids = DataFrame(randvoids[2])

gals_box = sphere2box(gals.ra, gals.dec, gals.chi)

ngals = length(gals.ra)

tree = KDTree(gals_box[1:100000])

function peebles_bruteforce(gals, voids, randgals, randvoids, N=20, R_min=0.1, R_max=2.0)
    dR = (R_max - R_min)/N

    dd_count = zeros(N)
    rr_count = zeros(N)

    for v in eachrow(voids)
        for g in eachrow(gals)
            d = distance(g.chi, v.chi, g.ra, v.ra, g.dec, v.dec)/v.Rv
            if d>R_max || d<R_min
                continue
            else
                jbin = floor(Int64,(d-R_min)/dR) + 1
                #println(jbin)
                dd_count[jbin] += 1
            end
        end
    end

    for v in eachrow(randvoids)
        for g in eachrow(randgals)
            d = distance(g.chi, v.chi, g.ra, v.ra, g.dec, v.dec)/v.Rv
            if d>R_max || d<R_min
                continue
            else
                jbin = floor(Int64,(d-R_min)/dR) + 1
                #println(jbin)
                rr_count[jbin] += 1
            end
        end
    end

    ξ = dd_count./rr_count .- 1.0
    return ξ
end

println(dd_count)
println(rr_count)

function calculate_vgcf()
    return    
end