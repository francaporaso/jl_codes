include("Constants.jl")
include("Cosmology.jl")
using .Constants
using .Cosmology

cosmo = Planck18

ncols = 3 #ra,dec,z
nrows = 1_000
randgals = rand(Float32, (nrows, ncols))

reduce(hcat, [randgals, comoving_distance.(cosmo, randgals[:,3])])