from utility_functions import import_traits_data_Seb, initialize_plant, weather_interp, run_trajectory, historical_rainfall
import numpy as np
import pandas as pd
import ipyparallel as ipp
import csv
import os

#gets engines running and checks to make sure they started up
rc=ipp.Client()
print(rc.ids)
print('checkpoint')
view=rc.load_balanced_view()
view.apply(lambda : "Hello, World")
async_results=[]


start_param = 0 #can adjust starting parameter set if running in batches
filename='elep_top_sets.txt' #file that contains parameter sets to run
param_factorial = np.loadtxt(filename)
param_sets = len(param_factorial)

#run info
start_doy = 121 #starting day - May 1 day 121
tmax = 488 #End Sept 1 following year 488
dt = 0.04 #time step [d] set for approx hourly
s0_upper = 0.363 #initial soil moisture in upper layer
s0_lower = 0.473 #initial soil moisture in lower layer

# set soil conditions
soil_type = 'tuned' #uses 'tuned' soil as base and overwrites with parameter values in each set

# set species traits
traits = import_traits_data_Seb(filepath='WA_hydraulic_traits_SebData.xls') #species traits, some are overwritten by parameters
plant = initialize_plant('Elep.3', traits, soil_type) #initialize plant with traits and soil


#import and format climate data to match simulation time series
dfWeather=pd.read_csv('weather_station_DOY_avg.csv')
VPD_raw=dfWeather.VPD
Irr_raw=dfWeather.Irr_W
Time_raw=dfWeather.Time
Irr_day_raw=dfWeather.Sum_MJ_day
T_raw=dfWeather.T_C
T_day_raw=dfWeather.T_avg_day_C
measurement_dt=Time_raw[1]-Time_raw[0]#input data on hourly basis
last_ind=int(len(Time_raw)-1)
start_ind = int((float(start_doy)-1.0)/(measurement_dt))

msmt_tRun=np.arange(0,tmax+measurement_dt,measurement_dt)
VPD_adj=np.zeros(len(msmt_tRun))
Irr_adj=np.zeros(len(msmt_tRun))
T_adj = np.zeros(len(msmt_tRun))
Irr_day_adj = np.zeros(len(msmt_tRun))
T_day_adj = np.zeros(len(msmt_tRun))

raw_index=start_ind

for i in range(len(msmt_tRun)):
    VPD_adj[i]=VPD_raw[raw_index]
    Irr_adj[i]=Irr_raw[raw_index]
    T_adj[i] = T_raw[raw_index]
    Irr_day_adj[i] = Irr_day_raw[raw_index]
    T_day_adj[i] = T_day_raw[raw_index]

    raw_index=raw_index+1
    if raw_index > last_ind:
        raw_index=0 #starts looping from beginning of year

new_sets = len((param_factorial[start_param:]))
tRun = np.arange(0, tmax + dt, dt)
slice_size = int(1.0/dt)
tRun_daily = tRun[::slice_size]
tRun_daily_sum = tRun_daily[:-1]
Amax = plant.canopy.Amax
R = plant.canopy.R()
original_gsmax = plant.canopy.Gs_leaf

#interpolate forcing data to dt of simulation
VPD_t = weather_interp(tRun, VPD_adj, msmt_tRun)
Irr_t = weather_interp(tRun, Irr_adj, msmt_tRun)
Irr_day_t = weather_interp(tRun, Irr_day_adj, msmt_tRun)
T_t = weather_interp(tRun, T_adj, msmt_tRun)
T_day_t = weather_interp(tRun, T_day_adj, msmt_tRun)

#get historical rainfall trajectories
Zr, n = plant.soil_root.L_root, plant.soil.n
depth_re, start_years=historical_rainfall(tmax, dt, start_doy, n, Zr)
n_trajectories = len(depth_re)

directory = "elep_3_trajectories"

parent_dir = "/global/scratch/users/jvwilkening/ET_cover_output"
path = os.path.join(parent_dir, directory)

np.savetxt(os.path.join(path, "tRun.txt"), tRun, delimiter=',')
np.savetxt(os.path.join(path, "tRun_daily.txt"), tRun_daily, delimiter=',')
np.savetxt(os.path.join(path, "tRun_daily_sum.txt"), tRun_daily_sum, delimiter=',')
np.savetxt(os.path.join(path, "normalized_rain_depths.txt"), depth_re, delimiter=',')
np.savetxt(os.path.join(path, "rain_start_years.txt"), start_years, delimiter=',')

for param_counter in range(new_sets):
    set_num = start_param + param_counter
    for yr_counter in range(n_trajectories):
        ar=view.apply_async(run_trajectory, plant, param_factorial, depth_re, start_years, tRun, dt, s0_upper, s0_lower, VPD_t,
                       Irr_t, T_t, Irr_day_t, T_day_t, original_gsmax, set_num, yr_counter, path)
        async_results.append(ar)


results = []

with open('completed_sets.csv', 'wb') as out:
    csv_out=csv.writer(out)
    for ar in async_results:
        res = ar.get()
        results.append(res)
        csv_out.writerow(res)
        out.flush() #updates csv file with each set of results
