using NearestNeighbors
using Statistics, LinearAlgebra
using Random
using Plots

include("Constants.jl")
include("Cosmology.jl")
include("Functools.jl")
using .Cosmology
using .Functools

function count_pairs_in_bins(tree, galaxy_pos::Matrix, void_centers::Matrix, void_radii::Vector;
                             RIN::Real=0.1, ROUT::Real=2.0, NBINS::Integer=20)
    
    counts = zeros(Int64, NBINS)
    nvoids = length(void_radii)
    dr = (ROUT-RIN)/NBINS

    #@inbounds for i in 1:nvoids
    for i in 1:nvoids
        rmin = RIN*void_radii[i]
        rmax = ROUT*void_radii[i]
        idx = inrange(tree, void_centers[:, i], rmax)
        if isempty(idx) continue end

        for j in idx
            dist = norm(galaxy_pos[:, j] - void_centers[:, i])
            jbin = floor(Int64, ((dist-rmin)/(dr*void_radii[i]))) + 1
            if jbin<1 continue end
            counts[jbin] += 1
        end
    end
    
    return counts
end

function galaxy_void_xcorr_peebles(galaxy_pos::Matrix, void_centers::Matrix, 
                                   void_radii::Vector, random_pos::Matrix;
                                   NBINS::Integer=20, RIN::Float64=0.5, ROUT::Float64=3.0)
    
    # Construct trees
    tree_true = KDTree(galaxy_pos)
    tree_rand = KDTree(random_pos)

    # Normalized bin edges
    r_norm_edges = range(RIN, ROUT, length=NBINS+1)
    r_norm_centers = [(r_norm_edges[i] + r_norm_edges[i+1]) / 2 for i in 1:NBINS]
    
    # Count DD pairs (Data-Data, i.e., galaxy-void pairs)
    DD = count_pairs_in_bins(
        tree_true, galaxy_pos, void_centers, void_radii, 
        RIN=RIN, ROUT=ROUT, NBINS=NBINS
    )
    
    # Count RR pairs (Random-Random, i.e., random-void pairs)
    RR = count_pairs_in_bins(
        tree_rand, random_pos, void_centers, void_radii,
        RIN=RIN, ROUT=ROUT, NBINS=NBINS
    )
    
    # Normalization factors
    n_gal = size(galaxy_pos, 2)
    n_rand = size(random_pos, 2)
    norm_factor = n_rand / n_gal
    
    # Peebles estimator: ξ = (DD/RR) * norm_factor - 1
    xi = zeros(Float64, NBINS)
    for i in 1:NBINS
        if RR[i] > 0
            xi[i] = (DD[i] / RR[i]) * norm_factor - 1.0
        else
            xi[i] = NaN
        end
    end
    
    return r_norm_centers, xi, DD, RR
end

function test_voidcat_random()
    Random.seed!(42)
    Rv_min, Rv_max, z_min, z_max, delta_min, delta_max = 10.0, 15.0, 0.2, 0.22, -1.0, -0.1
    RIN, ROUT, NBINS = 0.1, 2.0, 20
    cosmo = Cosmology.Planck18
    
    L, _, _ = lenscat_load(
        "/home/fmcaporaso/cats/L768/voids_LCDM_09.dat",
        Rv_min, Rv_max,
        z_min, z_max,
        delta_min, delta_max;
        fullshape=false
        )

    void_centers = sphere2box(cosmo, L[:,2], L[:,3], L[:,4])
    void_radii = L[:, 1]
    
    boxsize = sum(abs.(extrema(void_centers)))

    n_gals = 10000
    n_rand = n_gals*10
    galaxy_pos = rand(3, n_rand) .* boxsize
    random_pos = rand(3, n_rand) .* boxsize

    r, xi, _, _ = galaxy_void_xcorr_peebles(
        galaxy_pos, void_centers, void_radii, random_pos; RIN=RIN, ROUT=ROUT, NBINS=NBINS
    )

    plot(r, xi)

end

function test_random_random()
    Random.seed!(42)

    # Simulate data
    n_gal = 50000
    n_rand = 50000 * 5  # 5x more randoms for better statistics
    boxsize = 3000.0  # Mpc/h

    galaxy_pos = rand(3, n_gal) .* boxsize
    random_pos = rand(3, n_rand) .* boxsize

    n_voids = 500
    void_centers = rand(3, n_voids) .* boxsize
    void_radii = rand(n_voids) .* 40 .+ 10  # 10-50 Mpc/h

    RIN, ROUT, NBINS = 0.1, 2.0, 20

    r, xi, _, _ = galaxy_void_xcorr_peebles(
        galaxy_pos, void_centers, void_radii, random_pos; RIN=RIN, ROUT=ROUT, NBINS=NBINS
    )

    plot(r, xi)
end

test()