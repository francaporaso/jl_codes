include("cosmology.jl")

cosmo = LCDM(0.3089f0, 0.6911f0, 100f0)

function ang2xyz(ra::Float32, dec::Float32, redshift::Float32; cosmo::LCDM=cosmo)
    comdist = comoving_distance(cosmo, redshift)
    x = comdist*cosd(dec)*cosd(ra)
    y = comdist*cosd(dec)*sind(ra)
    z = comdist*sind(dec)
    return x,y,z
end

function ang2xyz(ra::Vector{Float32}, dec::Vector{Float32}, redshift::Vector{Float32}; cosmo::LCDM=cosmo)
    comdist = comoving_distance(cosmo, redshift)
    x = comdist.* cosd.(dec) .* cosd.(ra)
    y = comdist.* cosd.(dec) .* sind.(ra)
    z = comdist.* sind.(dec)
    return x,y,z
end

function mean_density_comovingshell(xh::Vector{Float32},yh::Vector{Float32},zh::Vector{Float32},logm::Vector{Float32},
                                    RMAX::Int32,rv::Float32,xv::Float32,yv::Float32,zv::Float32)

    dist_void = hypot(xv,yv,zv)
    dist = hypot.(xh,yh,zh)
    χ_min = dist_void - RMAX*rv
    χ_max = dist_void + RMAX*rvx

    lmh = @. logm[(dist>χ_min) && (dist<χ_max)]
    vol = Float32(4*pi/3)*(χ_max^3 - χ_min^3)

    return sum(x->1.0f1^x, lmh)/vol, length(lmh)/vol
end

function density(N::Int32,RMAX::Int32,xh::Vector{Float32},yh::Vector{Float32},zh::Vector{Float32},logm::Vector{Float32},
                 rv::Float32,xv::Float32,yv::Float32,zv::Float32)
    number_gx = zeros(Float32, N)
    mass_bin = zeros(Float32, N)
    vol = zeros(Float32, N)
    dist = @. hypot(xh-xv, yh-yv, zh-zv)
    c = RMAX*rv/N
    
    mean_den, mean_gx = mean_density_comovingshell(xh,yh,zh,logm,RMAX,rv,xv,yv,zv)
    
    mask_mean = dist .< 1.1f0*RMAX*rv
    logmass = logm[mask_mean]
    dist = dist[mask_mean]
    
    for k=0:N-1
        mask = (dist .< (k+1)*c) .&& (dist .>= k*c)
        number_gx[k+1] = sum(mask)
        mass_bin[k+1] = sum(x->1f1^x, logmass[mask], init=0f0)
        vol[k+1] = 3k*(k+1) + 1 # expansion de (k+1)^3 - k^3
    end

    vol *= Float32(4/3*pi)*c^3
    return number_gx, mass_bin, mean_den, mean_gx
end