using DataFrames
#using Plots
using DelimitedFiles
using Random
#using Distributions

include("Constants.jl")
include("Cosmology.jl")
using .Constants
using .Cosmology

cosmo = Planck18

function make_mock_galaxies(filename::AbstractString; nrows=1000, seed=0)
    Random.seed!(seed)
    ra = 360.0.*rand(Float32, nrows)
    dec = 90.0.*(2.0.*rand(Float32, nrows).-1)
    z = 1.4.*rand(Float32, nrows)
    χ = comoving_distance.(cosmo, z)
    randgals = DataFrame(ra=ra, dec=dec, z=z, chi=χ)
    
    open(filename, "w") do IO
        writedlm(IO, Iterators.flatten(([names(randgals)], eachrow(randgals))))
    end

end

function make_mock_voids(filename::AbstractString; nrows=100, seed=0)
    Random.seed!(seed)
    ra = 360.0.*rand(Float32, nrows)
    dec = 90.0.*(2.0.*rand(Float32, nrows).-1)
    z = 0.6.*rand(Float32, nrows)
    χ = comoving_distance.(cosmo, z)
    Rv = 8.0 .+ (30.0-8.0).*rand(Float32, nrows)

    randvoids = DataFrame(Rv=Rv, ra=ra, dec=dec, z=z, chi=χ)
    
    open(filename, "w") do IO
        writedlm(IO, Iterators.flatten(([names(randvoids)], eachrow(randvoids))))
    end

end

Ngals = 100_000
Nvoids = 1_000

make_mock_galaxies("mockgals.dat", seed=2, nrows=Ngals)
make_mock_voids("mockvoids.dat", seed=3, nrows=Nvoids)

make_mock_galaxies("random_gals.dat", seed=2, nrows=10*Ngals)
make_mock_voids("random_voids.dat", seed=3, nrows=10*Nvoids)

using Plots
function plot_mocks()

    randgals = readdlm("mockgals.dat")
    randgals = DataFrame(randgals[2:end, :], randgals[1,:])
    scatter(randgals.ra, randgals.dec)
    
    randvoids = readdlm("mockvoids.dat")
    randvoids = DataFrame(randvoids[2:end, :], randvoids[1,:])
    scatter!(randvoids.ra, randvoids.dec, markersize=randvoids.Rv, alpha=0.5)
    
end

plot_mocks()