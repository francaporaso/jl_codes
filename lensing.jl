include("Cosmology.jl")
include("Constants.jl")
using .Cosmology
using .Constants: sc_const

const cosmo = Planck18

function sigma_crit(cosmo::FlatΛCDM, z_l::Real, z_s::Real)
    d_l = angular_diameter_distance(cosmo, z_l)
    d_s = angular_diameter_distance(cosmo, z_s)
    d_ls = angular_diameter_distance_z1z2(cosmo, z_l, z_s)
    return sc_const*d_s/(d_ls*d_l)
end

function void_stacking()
    N = 10
    NK = 100
    NCORES = 16

    N_inbin = zeros(Float64, (NK+1, N))
    N_inbin = zeros(Float64, (NK+1, N))
    N_inbin = zeros(Float64, (NK+1, N))
    N_inbin = zeros(Float64, (NK+1, N))
end