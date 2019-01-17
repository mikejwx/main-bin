# -*- coding: utf-8 -*-
"""
Gibbs Moist Air Library: Functions to calculate various thermodynamic
properties of moist air. Formulae according to Ambaum "Thermal Physics of
the Atmosphere" (Wiley-Blackwell, 2010).

Inputs to all functions should be given in units:
    Temperature: C
    Pressures: mb
    Mixing ratios/specific humidities: kg/kg

Remi Tailleux/Bethan Harris, University of Reading, 11th November 2016
"""
import numpy as np

#==============================================================================
# Define useful thermodynamic constants
#==============================================================================

t0 = 273.15 #0 degrees Celsius in K
t25_k = 25. + 273.15 #25 degrees Celsius in K
cl = 4179.9 #specific heat capacity for liquid water (J/kg/K)
cpd = 1005.7 #specific heat at constant pressure for dry air (J/kg/K)
cpv = 1865.1 #specific heat capacity at constant pressure for water vapour (J/kg/K)
ci = 1960. #specific heat capacity at constant pressure for ice (J/kg/K)
e0_mb = 10. #reference partial vapour pressure for water vapour (mb)
es0_mb = 31.6743 #reference saturation vapour pressure(mb)
alv0 = 2.444e6 #latent heat of evaporation at T=25C (J/kg)
sub0 = 2.826e6 #latent heat of sublimation at T=0C (J/kg)
p0_mb = 1000. #reference pressure (mb)
pd0_mb = p0_mb - e0_mb #reference partial pressure of dry air (mb)
rd = 287.04 #specific gas constant for dry air (J/kg/K)
rv = 461.5 #specific gas constant for water vapour (J/kg/K)
eeps = rd/rv
chi = rd/cpd
cpvmcl = cl-cpv
alpha = cpvmcl/rv

def saturation_vapour_pressure(t_c):
    #Compute saturation water vapour pressure (mb) as a function of temperature
    t_k = t_c + 273.15
    dt_k = t_k - t25_k
    es_mb = (es0_mb*np.exp((alv0/rv)*dt_k/(t_k*t25_k))
    *(t25_k/t_k)**alpha*np.exp(alpha*dt_k/t_k))
    return es_mb
    
def saturation_mixing_ratio(t_c,p_mb):
    #Compute saturation mixing ratio (kg/kg) as a function of temperature and pressure
    es_mb = saturation_vapour_pressure(t_c)
    #Compute saturation mixing ratio as if everything was fine
    rvs_kgkg = eeps*es_mb/(p_mb - es_mb)
    #Set negative values to a very high values that should guarantee 
    #lack of saturation
    rvs_kgkg[rvs_kgkg < 0.] = 0.
    return rvs_kgkg

def latent_heat(t_c):
    #Compute latent heat of evaporation (J/kg)
    return alv0 - cpvmcl*(t_c-25.)
    
def latent_heat_sub(t_c):
    #Compute latent heat of sublimation (J/kg)
    return sub0 - (ci-cpv)*t_c
    
def rt_from_qt(qt_kgkg):
    #Convert total water content (specific humidity, kg/kg) into mixing ratio (kg/kg)
    return qt_kgkg/(1-qt_kgkg)
    
def liquid_mixing_ratio_old(t_c,p_mb,qt_kgkg):
    #Compute the liquid water mixing ratio (kg/kg) for saturated moist air
    rt_kgkg = rt_from_qt(qt_kgkg)
    rs_kgkg = saturation_mixing_ratio(t_c,p_mb)
    dr = rt_kgkg - rs_kgkg
    rl_kgkg = 0.5*(dr + abs(dr))
    return rl_kgkg
    
def partial_vapour_pressure(t_c,p_mb,qt_kgkg):
    #Compute the partial pressure of water vapour (mb)
    rt_kgkg = rt_from_qt(qt_kgkg)
    e_mb = rt_kgkg*p_mb/(eeps+rt_kgkg)
    es_mb = saturation_vapour_pressure(t_c)
    e_mb = np.minimum(e_mb,es_mb)
    return e_mb

def entropy_dry_air(t_c,pd_mb):
    #Compute the partial specific entropy for dry air (J/kg/K)
    t_k = t_c + t0
    eta_d = cpd*np.log(t_k/t0) - rd*np.log(pd_mb/pd0_mb)
    return eta_d
    
def entropy_water_vapour(t_c,e_mb):
    #Compute the partial specific entropy for water vapour (J/kg/K)
    t_k = t_c + t0
    e_mb_term = rv*np.log(e_mb/e0_mb)
    eta_v = cpv*np.log(t_k/t0) - e_mb_term
    #Set -Inf values to 0
    eta_v[np.isinf(e_mb_term)]= 0
    return eta_v    

def moist_entropy(t_c,p_mb,qt_kgkg):
    #Compute the moist specific entropy (J/kg/K) assuming phase liquid equilibrium
    t_k = t_c + t0
    rl_kgkg = liquid_mixing_ratio(t_c,p_mb,qt_kgkg)
    ql_kgkg = rl_kgkg*(1-qt_kgkg)
    e_mb = partial_vapour_pressure(t_c,p_mb,qt_kgkg)
    pd_mb = p_mb - e_mb
    lv = latent_heat(t_c)
    eta_ma = ((1-qt_kgkg)*entropy_dry_air(t_c,pd_mb)
    + qt_kgkg*entropy_water_vapour(t_c,e_mb) - ql_kgkg*(lv/t_k))
    return eta_ma
    
def moist_entropy_ice(t_c,p_mb,qt_kgkg):
    #Compute the moist specific entropy (J/kg/K) assuming phase liquid equilibrium

    t_k = t_c + t0
    entropy = np.zeros_like(t_c,dtype='float32')
    rl_kgkg = liquid_mixing_ratio(t_c,p_mb,qt_kgkg)
    ql_kgkg = rl_kgkg*(1-qt_kgkg)
    e_mb = partial_vapour_pressure(t_c,p_mb,qt_kgkg)
    pd_mb = p_mb - e_mb
    
    lv = latent_heat(t_c)
    eta_ma = ((1-qt_kgkg)*entropy_dry_air(t_c,pd_mb)
    + qt_kgkg*entropy_water_vapour(t_c,e_mb) - ql_kgkg*(lv/t_k))
    
    #Compute entropy if all liquid is frozen to ice
    ls = latent_heat_sub(t_c)
    eta_ice = ((1-qt_kgkg)*entropy_dry_air(t_c,pd_mb)
    + qt_kgkg*entropy_water_vapour(t_c,e_mb) - ql_kgkg*(ls/t_k))
    
    ice = np.where(t_c<0.)
    liquid = np.where(t_c>=0.)
    entropy[ice] = eta_ice[ice]
    entropy[liquid] = eta_ma[liquid]
    
    return entropy
    
def saturation_vapour_pressure_deriv(t_c):   
    #Compute T derivative of Ambaum's saturation vapour pressure formula
    #d(es)/dT = aT^(-d-2)*exp(-(c-T)(b+cd)/cT)*...
    #(bT^d*(c/T)^d-d*Tc^d + cdx^d*(c/T)^d)
    #where a = es0_mb
    #b = alv0/rv
    #c = 298.15 K
    #d = alpha = (cpl-cpv)/rv
    t_k = t_c + t0
    a = es0_mb
    b = alv0/rv
    c = t25_k
    d = alpha
    T = t_k
    term1 = a*T**(-d-2)
    term2 = ((c-T)*(b+c*d))/(c*T)
    term3 = b*T**d*(c/T)**d - d*T*c**d + c*d*T**d*(c/T)**d
    deriv = term1*np.exp(-term2)*term3
    return deriv
    
def cp_moist_air_exact(t_c,p_mb,qt_kgkg):
    #Compute effective specific heat capacity (J/kg/K) for moist air
    rl_kgkg = liquid_mixing_ratio(t_c,p_mb,qt_kgkg);
    ql_kgkg = rl_kgkg*(1-qt_kgkg)
    qv_kgkg = qt_kgkg - ql_kgkg
    
    e_mb = partial_vapour_pressure(t_c,p_mb,qt_kgkg)
    pd_mb = p_mb - e_mb
    
    tk = t_c + t0
    lv = latent_heat(t_c)
    qvlvdrdt = (qv_kgkg/rv)*(p_mb/pd_mb)*(lv/tk)**2
    
    cp_moist = cpd*(1-qt_kgkg) + cpv*qv_kgkg
    isat = rl_kgkg > 0.
    cp_moist[isat] = cp_moist[isat] + cl*ql_kgkg[isat] + qvlvdrdt[isat]
    return cp_moist
    
def liquid_mixing_ratio(t_c,p_mb,qt_kgkg):
    #Compute the liquid water mixing ratio (kg/kg) for saturated moist air
    #Compute partial pressure assuming that air is unsaturated
    rt_kgkg = rt_from_qt(qt_kgkg)
    e_mb = rt_kgkg*p_mb/(eeps+rt_kgkg)
    #Compute saturation water pressure
    es_mb = saturation_vapour_pressure(t_c)
    #Define partial water vapour pressure as the mininum of the two pressure
    e_mb = np.minimum(e_mb,es_mb)
    #Compute water vapour and liquid mixing ratio
    rv_kgkg = eeps*e_mb/(p_mb-e_mb)
    rl_kgkg = rt_kgkg - rv_kgkg
    return rl_kgkg
    
def spec_vol(t_c,p_mb,qt_kgkg):
    #Compute the specific volume (m^3/kg) of moist air
    t_k = t_c + t0
    #Compute partial pressure of dry air
    pd_mb = p_mb - partial_vapour_pressure(t_c,p_mb,qt_kgkg)
    pd_ppa = pd_mb*100;
    #Compute specific volume according to alpha = (1-qt)*alpha_d
    spec_vol = (1-qt_kgkg)*rd*t_k/pd_ppa;
    return spec_vol
    
def moist_enthalpy(t_c,p_mb,qt_kgkg):
    #Compute specific moist enthalpy (J/kg/K) assuming phase liquid equilibrium
    #Compute mixing ratio and specific humidity for the liquid part
    rl_kgkg = liquid_mixing_ratio(t_c,p_mb,qt_kgkg)
    ql_kgkg = rl_kgkg*(1-qt_kgkg)
    #Compute latent heat
    lv = latent_heat(t_c)
    #Estimate moist entropy as weighted partial entropies, accounting for
    #liquid mixing ratio contribution separately for saturated case
    cp_eff = cpd*(1-qt_kgkg) + cpv*qt_kgkg
    h_ma = cp_eff*t_c - ql_kgkg*lv
    return h_ma
    
def spec_hum_from_rel_hum(t_c, p_mb, rh):
    #Calculate specific humidity (kg/kg) from relative humidity (frac)
    es_mb = saturation_vapour_pressure(t_c)
    e_mb = rh*es_mb
    spec_hum = e_mb/(e_mb+(p_mb-e_mb)/eeps)
    return spec_hum
    
def dmudq(t_c,qt_kgkg):
    t_k = t_c + t0
    dmudq = rv*t_k/qt_kgkg + rd*t_k/(1-qt_kgkg)
    return dmudq
    
def diffusion_coeff(t_c):
    Dv = 21.2e-6*(1+0.0071*t_c)
    return Dv
    
def gibbs_dry_air(t_c, pd_mb):
    t_k = t_c + t0
    gd = cpd*(t_k-t0-t_k*np.log(t_k/t0)) + rd*t_k*np.log(pd_mb/pd0_mb)
    return gd
    
def gibbs_water_vapour(t_c, e_mb):
    t_k = t_c + t0
    gv = cpv*(t_k-t0-t_k*np.log(t_k/t0)) + rv*t_k*np.log(e_mb/e0_mb)
    return gv
    
def chemical_potential(t_c, p_mb,qt_kgkg):
    e_mb = partial_vapour_pressure(t_c,p_mb,qt_kgkg)
    pd_mb = p_mb - e_mb
    mu = gibbs_water_vapour(t_c,e_mb) - gibbs_dry_air(t_c,pd_mb)
    return mu
    
def dry_potential_temperature(t_c,p_mb):
    theta = (t_c+t0)*(p0_mb/p_mb)**(rd/cpd)
    return theta
    
def equivalent_potential_temperature(t_c, p_mb, qt_kgkg):
    theta = dry_potential_temperature(t_c, p_mb)
    #rvs = saturation_mixing_ratio(t_c, p_mb)
    rt = rt_from_qt(qt_kgkg)
    rw = liquid_mixing_ratio(t_c, p_mb, qt_kgkg)
    rv = rt - rw
    es = saturation_vapour_pressure(t_c)
    a = np.exp(alv0*rv/(cpd*(t_c+t0)))
    b = ((t_c+t0)/t0)**(rw*cl/cpd)
    c = (1 - es/p_mb)**(-rd/cpd)
    thetae = theta*a*b*c
    return thetae
    
def saturated_equivalent_potential_temperature(t_c, p_mb, qt_kgkg):
    theta = dry_potential_temperature(t_c, p_mb)
    rvs = saturation_mixing_ratio(t_c, p_mb)
   # rt = rt_from_qt(qt_kgkg)
    rw = liquid_mixing_ratio(t_c, p_mb, qt_kgkg)
    #rv = rt - rw
    es = saturation_vapour_pressure(t_c)
    a = np.exp(alv0*rvs/(cpd*(t_c+t0)))
    b = ((t_c+t0)/t0)**(rw*cl/cpd)
    c = (1 - es/p_mb)**(-rd/cpd)
    thetaes = theta*a*b*c
    return thetaes