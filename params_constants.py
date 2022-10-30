
v = 1.557           # conversion from mol/s to m3/day for liquid water
a = 1.6             # ratio of diffusivity water vapor/CO2
P_baro = 97.42      # kPa, atmospheric pressure at 330 m

nu = 10**6          # unit conversion Pa/MPa
rhoH2O = 1000.0     # density of water, kg/m3
g = 9.81            # acceleration due to gravity, m/s2
Mw = 0.018          # molecular weight of water, kg/mol
rhoCO2 = 1.98       # density of CO2, kg/m3
Mco2 = 0.044        # molecular weight of CO2, kg/mol, ...so 45mol/m3
#Sd = 36000          # day lengths in seconds, assume 10 hrs of usable daylight
Sd = 60*60*24          # For shorter temporal resolution use total number of seconds in day, irr response will take care of daylight
CO2conc = 0.000400  # in mol/mol
g_aero = 0.02       #aerodynamic conductance to water vapor transport from Oleson, KW et al (2008) via Liu et al (2017) m/s
irr_b_o = 2.0       # decay coefficient for shortwave radiation as a function of LAI [-]
cp_air = 1005.0     # Constant pressure specific heat capacity of air [J/Kg C]
lambda_w = 2450000.0 # Latent heat of vaporization of water [J/Kg]
Mair = 0.02897       #molecular mass of air [Kg/mol]
rhoAir = 1.2250     # Density of air [Kg /m3]