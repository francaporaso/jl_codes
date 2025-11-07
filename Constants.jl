module Constants
    
export c, G, pc, M_sun, sc_const

const c = 299792458.0 #m/s
const G = 6.6743e-11 #m^3/(kg s^2)
const pc = 3.085677581491367e+16 #m
const M_sun = 1.988409870698051e+30 #kg
const sc_const = (c^2/(4π*G)) * (pc/M_sun)*1e-6

end