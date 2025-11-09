using DataFrames
using DelimitedFiles

include("Cosmology.jl")
using .Cosmology

cosmo = Planck18

gals = readdlm("mockgalaxies.dat")
gals = DataFrame(gals[2:end,:], gals[1,:])

randgals = readdlm("random_gals.dat")
randgals = DataFrame(randgals[2:end,:], randgals[1,:])

voids = readdlm("mockvoids.dat")
voids = DataFrame(voids[2:end,:], voids[1,:])

randvoids = readdlm("random_voids.dat")
randvoids = DataFrame(randvoids[2:end,:], randvoids[1,:])

N = 20
R_min = 0.1
R_max = 5.0

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

println(dd_count)
println(rr_count)

function distance(χ1, χ2, α1, α2, δ1, δ2)
    return sqrt(χ1^2+χ2^2- 2*χ1*χ2*cos_angular_separation(α1, α2, δ1, δ2))
end

function cos_angular_separation(α1, α2, δ1, δ2)
    return sin(δ1)*sin(δ2) + cos(δ1)*cos(δ2)*cos(α1-α2)
end

function calculate_vgcf()
    return    
end