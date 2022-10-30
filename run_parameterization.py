from utility_functions import import_traits_data_Seb, initialize_plant, weather_interp, run_full_param
from params_soil import soil_dict
import observation_data
import numpy as np
import pandas as pd
import ipyparallel as ipp
import csv

rc=ipp.Client()
print(rc.ids)
print('checkpoint')
view=rc.load_balanced_view()
view.apply(lambda : "Hello, World")
async_results=[]

start_param = 45000
filename='elep_twolayer_parameters_full_111121.txt'
param_factorial = np.loadtxt(filename)

#run info
tmax=40
s0=.53
dt=.04
obs_start_day=37.0
days_of_obs = 3

# set soil conditions
soil_type = 'tuned'
#smin = soil_dict[soil_type]['sh']
#sfc = soil_dict[soil_type]['sfc']
#sst = soil_dict[soil_type]['sst']
#n = soil_dict[soil_type]['n']
#Ks = soil_dict[soil_type]['Ksat']

# set species traits
traits = import_traits_data_Seb(filepath='WA_hydraulic_traits_SebData.xls')
plant = initialize_plant('Elep.3', traits, soil_type)

# load observation data
dfelep3 = pd.read_csv('Elep_3_obs.csv')
obs_g = dfelep3.Mean_gs
obs_g_t = dfelep3.Timestamp
obs_A = dfelep3.Mean_A
obs_A_t = dfelep3.Timestamp
obs_sfinal_upper = observation_data.elep_3_sfinal_upper
obs_sfinal_lower = observation_data.elep_3_sfinal_lower
LWP = observation_data.elep_3_LWP
LWP_t = observation_data.elep_3_LWP_t

##Post-rain meteorological data
dfweather = pd.read_csv('weather_station.csv')
VPD_m = dfweather.VPD
VPD_mt = dfweather.raw_t
Irr_m = dfweather.Irr_W
Irr_mt = dfweather.raw_t
T_m = dfweather.AirTemp_Avg
T_mt = dfweather.raw_t

new_sets = len((param_factorial[start_param:]))
param_scores = np.zeros((new_sets, 16))
tRun = np.arange(0, tmax + dt, dt)
Amax = plant.canopy.Amax
R = plant.canopy.R()
depth_re = np.zeros(len(tRun))  # no rainfall
original_gsmax = plant.canopy.Gs_leaf

#interpolate climate data
VPD_t = weather_interp(tRun, VPD_m, VPD_mt)
Irr_t = weather_interp(tRun, Irr_m, Irr_mt)
T_t = weather_interp(tRun, T_m, T_mt)

#placeholders for daily data for bare soil evap - only uses if plants die
Irr_day = Irr_t
T_day = T_t

print(new_sets)

for j in range(new_sets):
    set_num = start_param + j
    ar=view.apply_async(run_full_param, plant, param_factorial, obs_g, obs_g_t, obs_A, obs_A_t, obs_sfinal_upper, obs_sfinal_lower, LWP, LWP_t, depth_re, tRun, dt, s0,
              VPD_t, Irr_t, T_t, Irr_day, T_day, original_gsmax, set_num, obs_start_day, days_of_obs)
    async_results.append(ar)



results = []

with open('parameterization_scores_elep3_twolayer.csv', 'wb') as out:
    csv_out=csv.writer(out)
    for ar in async_results:
        res = ar.get()
        results.append(res)
        csv_out.writerow(res)
        out.flush() #updates csv file with each set of results



