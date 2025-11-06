module Cosmology

export ΛCDM, FlatΛCDM, comoving_distance, angular_diameter_distance, luminosity_distance,
    angular_diameter_distance_z1z2

using QuadGK

const c = 299792.458 #km/s

abstract type AbstractCosmology end
abstract type FLRW <: AbstractCosmology end
abstract type AbstractΛCDM <: FLRW end
abstract type AbstractwCDM <: FLRW end

struct ΛCDM{T<:Real} <: AbstractΛCDM
    Ω_m0::T
    Ω_Λ0::T
    h::T
    Ω_r0::T
    Ω_k0::T
end

struct FlatΛCDM{T<:Real} <: AbstractΛCDM
    Ω_m0::T
    h::T
    Ω_Λ0::T
    Ω_r0::T
    Ω_k0::T
    
    function FlatΛCDM(Ω_m0::T, h::T; Ω_Λ0::T=zero(T), Ω_r0::T=zero(T), Ω_k0::T=zero(T)) where {T<:Real}
        Ω_Λ0 = 1 - Ω_m0 - Ω_r0
        new{T}(Ω_m0, h, Ω_Λ0, zero(h), zero(h))
    end
end

Base.broadcastable(c::AbstractCosmology) = Ref(c)

function invE(cosmo::AbstractΛCDM, z::Real)
    return 1.0/sqrt(cosmo.Ω_r0*(1.0+z)^4 + cosmo.Ω_m0*(1.0+z)^3 
                    + cosmo.Ω_k0*(1.0+z)^2 + cosmo.Ω_Λ0)
end

function comoving_distance(cosmo::AbstractΛCDM, z::Real)
    integral, _ = quadgk(z->invE(cosmo, z), 0.0, z)
    return 0.01*c/cosmo.h * integral
end

function angular_diameter_distance(cosmo::FlatΛCDM, z::Real)
    return comoving_distance(cosmo, z)*(1.0+z)^(-1)
end

function luminosity_distance(cosmo::FlatΛCDM, z::Real)
    return comoving_distance(cosmo, z)*(1.0+z)
end

function angular_diameter_distance_z1z2(cosmo::FlatΛCDM, z1::Real, z2::Real)
    D_1 = comoving_distance(cosmo, z1)
    D_2 = comoving_distance(cosmo, z2)
    D_H = 0.01*c/cosmo.h
    D_12 = 1.0/(1.0+z2) * (D_2*sqrt(1.0+cosmo.Ω_r0*(D_1/D_H)^2) - D_1*sqrt(1.0+cosmo.Ω_r0*(D_2/D_H)^2))
    return D_12
end

end #Cosmology