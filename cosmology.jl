using QuadGK

struct LCDM
    Om0::Float32
    Ode0::Float32
    H0::Float32
end

function E(cosmo::LCDM, redshift::Float32)
    return sqrt(cosmo.Om0*(1+redshift)^3 + cosmo.Ode0)
end

function E(cosmo::LCDM, redshift::Vector{Float32})
    return @. sqrt(cosmo.Om0 * (1f0+redshift)^3 + cosmo.Ode0)
end

function comoving_distance(cosmo::LCDM, redshift::Float32)
    DH = 299792.458f0/cosmo.H0 # Mpc or Mpc/h if H0 = 100
    integral, _ = quadgk( z->1/E(cosmo, z), 0f0, redshift)
    return DH*integral
end

function comoving_distance(cosmo::LCDM, redshift::Vector{Float32})
    DH = 299792.458f0/cosmo.H0 # Mpc or Mpc/h if H0 = 100
    individual = zeros(Float32, redshift.size)
    for (i,zi) in enumerate(redshift)
        individual[i] = quadgk( z->1/E(cosmo, z), 0f0, zi)[1]
    end 
    return DH.*individual
end