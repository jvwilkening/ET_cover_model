import xlrd
import numpy as np
import pandas as pd
import os
from soil_plant_model import Canopy, Stem, Soil_root, Whole_plant, Soil
from params_soil import soil_dict


def import_traits_data_Seb(filepath=None, species=None, sp_coln=None):
    if filepath == None: filepath = '../WA_hydraulic_traits_SebData.xls'
    if species == None: species = ['Elep.3', 'Elox.3', 'Elep.5', 'Elox.5', 'Elep.7', 'Elox.7', 'Elox.3.small', 'Elox.3.large']
    if sp_coln == None: sp_coln = [2, 3, 4, 5, 6, 7, 8, 9]
    book = xlrd.open_workbook(filepath)
    sheet = book.sheet_by_name('parameters')
    keys = np.asarray(filter(None, sheet.col_values(0)), dtype='str')
    canopy_keys = ['A_canopy', 'Gs_leaf', 'Amax', 'rho', 'lamABA', 'm']
    stem_keys = ['L_stem', 'A_stem', 'Ksat_stem', 'a_stem', 'P50_stem']
    root_keys = ['L_root_tot', 'L_root_1', 'f_root_1', 'A_root_mult', 'd_root', 'RAI', 'A_crown', 'LAI']
    chap_dict = {}
    for sp, spc in zip(species, sp_coln):
        chap_dict[sp] = {}
        for part_keys, part_dict in zip([canopy_keys, stem_keys, root_keys], ['canopy_dict', 'stem_dict', 'root_dict']):
            chap_dict[sp][part_dict] = {}
            for key in part_keys:
                j = np.where(keys == key)[0][0] + 1  # to account for the column of the species
                chap_dict[sp][part_dict][key] = sheet.col_values(spc)[j]
    return chap_dict

def initialize_plant(sp, params, soil_type):
    # initializing each species with its name (in strings)
    canopy_dict = params[sp]['canopy_dict']
    stem_dict = params[sp]['stem_dict']
    root_dict = params[sp]['root_dict']
    plant = Whole_plant(species=sp)
    plant.canopy = Canopy(**canopy_dict)
    plant.stem = Stem(**stem_dict)
    plant.soil_root = Soil_root(soil_type=soil_type, **root_dict)
    plant.soil = Soil(soil_type, **root_dict)
    return plant

def simulate_rainfall(n_trajectories, tRun, dt, lam, gam):
    size = len(tRun) * n_trajectories
    depthExp = -np.log(1.0 - np.random.random(size=size)) / gam
    freqUnif = np.random.random(size=size)

    depth = np.zeros(size)
    # the occurence of rainfall in any independent interval is lam*dt
    depth[freqUnif < np.tile(lam, size) * dt] = depthExp[
        freqUnif < np.tile(lam, size) * dt]  # rain falls according to prob within an increment
    depth_re = np.reshape(depth, (n_trajectories, len(tRun)))
    return depth_re

def simulate_s_t_norefilling_varVPD(depths, tRun, dt, sInit_upper, sInit_lower, plant, VPD_t, Irr_t, T_t, Irr_Day, T_Day):
    # need to incorporate the case of NO refilling (refilling==False)
    Px_death = np.log((1.0/0.8)-1.0)/plant.stem.a_stem + plant.stem.P50_stem  # Stem water potential that corresponds to 80% PLC
    L_sum = 0.0  # sum of leakage
    s_t_upper = np.zeros(len(tRun))
    s_t_lower = np.zeros(len(tRun))
    s0_upper = sInit_upper
    s0_lower = sInit_lower
    assm_t = np.zeros_like(s_t_upper)
    assm_t[0] = 0
    px_t = np.zeros_like(s_t_upper)
    pxmin_t = np.zeros_like(s_t_upper)
    E_t = np.zeros_like(s_t_upper)
    gs_t = np.zeros_like(s_t_upper)
    pl_t = np.zeros_like(s_t_upper)
    Arate_t = np.zeros_like(s_t_upper)
    ABA_t = np.zeros_like(s_t_upper)
    Lsum_t = np.zeros_like(s_t_upper)
    k_loss = np.zeros_like(s_t_upper)
    VPD0 = VPD_t[0]
    Irr0 = Irr_t[0]

    P50 = plant.stem.P50_stem
    a_stem = plant.stem.a_stem
    R = plant.canopy.R()  # respiration rate as fraction of Amax
    R_sec = plant.canopy.R_rate()  # respiration date on second
    Ar = plant.soil_root.A_root_mult * plant.soil_root.A_crown
    Zr_tot, n, Zr_1, Zr_2 = plant.soil_root.L_root, plant.soil.n, plant.soil_root.L_root_1, plant.soil_root.L_root_2
    Ksat, b = plant.soil.Ksat_soil, plant.soil.b

    Px_min = plant.get_fluxes_scalar(VPD0, s0_upper, s0_lower, 0.0, Irr0)[1]
    px_t[0] = Px_min
    plant_alive = 0  # switch for plant death, 0=alive, 1=dead
    first_daybreak = 0  # first daybreak solution fed 1=yes 0=no
    first_daybreak2 = 0  # second time step daybreak 1=yes 0=no
    first_sunset = 0
    Irr_prev = 0.0
    Irr_2prev = 0.0
    time_of_death = max(tRun)

    for i in range(len(tRun)):
        if plant_alive == 0: #when plant is alive assume 22% interception under canopy
            R_normed_1 = (depths[i] * Zr_tot / Zr_1)*((0.78*(1.0/plant.soil_root.A_root_mult))+1.0*(1.0-(1.0/plant.soil_root.A_root_mult)))  # change to account for nsim, depths are normed to total depth
        else: #when plant is dead assume no interception
            R_normed_1 = depths[i] * Zr_tot / Zr_1
        if R_normed_1 > 1.0-s0_upper: #if rainfall exceeds storage capacity in upper layer
            Infil_normed_1 = 1.0 - s0_upper #upper layer infiltration
            R_excess_1 = R_normed_1 - Infil_normed_1 #excess from layer one in terms of layer 1 saturation
            R_excess_1_depth = R_excess_1 * n * Zr_1 #excess from layer one in terms of m water
            R_normed_2 = R_excess_1_depth / (n * Zr_2) #excess to layer 2 in terms of layer 2 saturation
            if R_normed_2 > 1.0-s0_lower: #still exceeds storage in bottom layer
                Infil_normed_2 = (1.0-s0_lower) #what stays in layer 2
                immediate_drainage = (R_normed_2 - Infil_normed_2) * n * Zr_2 #excess goes to deep drainage, convert to m
            else:
                Infil_normed_2 = R_normed_2
                immediate_drainage = 0.0
        else: #incoming precip does not exceed upper level storage
            Infil_normed_1 = R_normed_1
            Infil_normed_2 = 0.0
            immediate_drainage = 0.0

        if i > 1:
            Irr_prev = Irr_t[i - 1]
            Irr_2prev = Irr_t[i - 2]
        VPD = VPD_t[i]
        Irr = Irr_t[i]
        T = T_t[i]

        if (Irr_prev < 0.0001) & (Irr > 0.0) & (
                first_daybreak == 1):  # first light uses previous first light solution to feed solver
            prev_sol = prev_sol_daybreak
        if (Irr_2prev < 0.0001) & (Irr > 0.0) & (
                first_daybreak2 == 1):  # first light uses previous first light solution to feed solver
            prev_sol = prev_sol_daybreak2
        if (Irr_prev > 0.0) & (Irr < 0.0001) & (
                first_sunset == 1):  # sunset uses previous sunset solution to feed solver
            prev_sol = prev_sol_sunset
        # if tRun[i] >1.92:
        # pdb.set_trace()

        # Translate s into soil water potential P_soil
        s_upper_1 = s0_upper + Infil_normed_1  # update with rainfall
        s_lower_1 = s0_lower + Infil_normed_2
        P_soil_upper = plant.soil.P_soil_solver(s_upper_1)
        P_soil_lower = plant.soil.P_soil_solver(s_lower_1)

        if i > 0:
            if prev_sol[1] > P_soil_lower and prev_sol[1] > P_soil_upper:
                if Irr < 1.0: #if dark guess plant in equilibrium with most positive soil layer
                    prev_sol[1] = max(P_soil_lower, P_soil_upper)
                    _, prev_sol[2], prev_sol[4], prev_sol[3], prev_sol[0] = plant.calculate_from_Px(prev_sol[1], VPD, Px_min, P_soil_upper, P_soil_lower, Irr)
                else: #if not dark guess 10% offfset from most positive soil layer
                    prev_sol[1] = 1.1*max(P_soil_lower, P_soil_upper)
                    _, prev_sol[2], prev_sol[4], prev_sol[3], prev_sol[0] = plant.calculate_from_Px(prev_sol[1], VPD, Px_min,
                                                                                       P_soil_upper, P_soil_lower, Irr)
        if plant_alive == 0:  # plant still alive
            # Calculate fluxes and plant water potentials
            if i > 0:
                Trans, Px, _, _, _, _ = plant.get_fluxes_scalar(VPD, s_upper_1, s_lower_1, Px_min, Irr, prev_sol)
            else:
                Trans, Px, _, _, _, _ = plant.get_fluxes_scalar(VPD, s_upper_1, s_lower_1, Px_min, Irr)
            if Trans < 0.000001: #ie plant not actively transpiring
                J=0.0
                J_upper = 0.0
                J_lower = 0.0
                Px = max(P_soil_upper, P_soil_lower)
                if i > 0:
                    #Pl = prev_sol[2]
                    Pl = Px
                    ABA = prev_sol[4]
                else:
                    _, _, Pl, _, ABA, _ = plant.get_fluxes_scalar(VPD, s_upper_1, s_lower_1, Px_min, Irr)
                gs = 0.0
            else:
                J_soil = plant.ET_soil_func_2layer(P_soil_upper, P_soil_lower, Px)
                J= Trans
                if J > 0.95*J_soil and J < 1.05*J_soil: #check that transpiration based on full system and soil are equal
                    J_upper = min(plant.ET_soil_func_upper_layer(P_soil_upper, Px), Trans)
                    J_lower = J - J_upper
                else: #if discrepany default to defining fluxes based on soil and stem potential
                    if J_soil > 0.25: #unreasonable J_soil estimate
                        if P_soil_upper > P_soil_lower: #prioritize transpiration from wetter layer
                            J_upper = min(plant.ET_soil_func_upper_layer(P_soil_upper, Px), Trans)
                            J_lower = J - J_upper
                        else:
                            J_lower = min(plant.ET_soil_func_lower_layer(P_soil_lower, Px), Trans)
                            J_upper = J - J_lower
                    else:
                        J = J_soil
                        J_upper = plant.ET_soil_func_upper_layer(P_soil_upper, Px)
                        J_lower = plant.ET_soil_func_lower_layer(P_soil_lower, Px)
                Pl = plant.Pl_ET_func(J, Px, Px_min)
                if Pl > 0.0:
                    Pl = prev_sol[2]

                # calculate the rest of the variables
                ABA = plant.ABA_func(Pl, J)
                gs = plant.stomata_func2(ABA, VPD, Pl, Irr)   # stomatal conductance for CO2
            # use soil moisture water balance
            E_soil = 0.0 #no bare soil evaporation
            L_upper = Ksat * s0_upper ** (2 * b + 3) #Leakage rate upper layer to lower layer[m/d]
            L_lower = Ksat * s0_lower ** (2 * b + 3) #Leakage rate lower layer to deep drainage[m/d]
            E_upper = J_upper / (n * Zr_1 * Ar)
            E_lower = J_lower / (n * Zr_2 * Ar)
            s_upper_out = max(s_upper_1 - dt * (E_upper + E_soil) - dt * (L_upper / (n * Zr_1)), 10 ** (-20))  # drawdown
            s_lower_out = max(s_lower_1 - dt * (E_lower) - dt *((L_lower - L_upper) / (n * Zr_2)), 10 ** (-20))  # drawdown
            if s_lower_out > 1.0: #if storage capacity is exceeded
                s_lower_excess = s_lower_out - 1.0
                s_lower_out =1.0
            else:
                s_lower_excess = 0.0
            ASM = plant.assimilation_daily(gs)
            A_rate = plant.assimilation_rate(gs)  # instantaneous assimilation rate
            L_sum = L_sum + (L_lower * dt) + (s_lower_excess * (n * Zr_2) + immediate_drainage)  # cumulative leakage m
        if plant_alive == 1:  # plant has died
            Px = Px_death
            J_upper = 0.0
            J_lower = 0.0
            Pl = 0.0
            ABA = 0.0
            gs = 0.0
            ASM = 0.0
            A_rate = 0.0
            R = 0.0  # no respiration once dead
            R_sec = 0.0
             #soil evaporation is calculated on daily scale so only do at beginning of day
            cur_time=np.floor(tRun[i])
            if cur_time == 0.0:
                prev_time = -1.0
            else:
                prev_time = np.floor(tRun[i-1])
            if cur_time > prev_time: #first time step of new day
                T_d = T_Day[i] #Average temp over whole day C
                Irr_d = Irr_Day[i] #Total daily solar irradiation MJ/m2
                E_soil = plant.soil.soil_evap(T_d, Irr_d, s0_upper) / (n * Zr_1) #soil evaporation (top layer only) output m/d
            else:
                E_soil=0.0
            L_upper = Ksat * s0_upper ** (2 * b + 3)
            L_lower = Ksat * s0_lower ** (2 * b + 3)
            E_upper = J_upper / (n * Zr_1 * Ar)
            E_lower = J_lower / (n * Zr_2 * Ar)
            s_upper_out = max(s_upper_1 - dt * (E_upper) - E_soil - dt * (L_upper / (n * Zr_1)), 10 ** (-20))# drawdown
            s_lower_out = max(s_lower_1 - dt * (E_lower) - dt *((L_lower - L_upper) / (n * Zr_2)), 10 ** (-20))# drawdown

            if s_lower_out > 1.0:
                s_lower_excess = s_lower_out - 1.0
                s_lower_out=1.0
            else:
                s_lower_excess = 0.0
            L_sum = L_sum + (L_lower * dt) + (s_lower_excess * (n * Zr_2) + immediate_drainage)  # cumulative leakage m
        # update Px_min
        Px_min = min(Px, Px_min)

        E = J_upper + J_lower
        # use previous step as initial conditions for next step
        prev_sol = [E, Px, Pl, gs, ABA]

        if (Irr_prev < 0.0001) & (Irr > 0.0):  # first light uses previous first light solution to feed solver
            prev_sol_daybreak = prev_sol
            first_daybreak = 1
        if (Irr_2prev < 0.0001) & (Irr > 0.0):  # first light uses previous first light solution to feed solver
            prev_sol_daybreak2 = prev_sol
            first_daybreak2 = 1

        if (Irr < 0.0001) & (Irr_prev > 0.0):  # sunset uses previous sunset solution to feed solver
            prev_sol_sunset = prev_sol
            first_sunset = 1

        if Px_min <= Px_death:  # if reached 80% PLC changes plant death switch
            plant_alive = 1
            Px_min = Px_death
            time_of_death = min(time_of_death, tRun[i])
        else:
            # update Px_min
            Px_min = min(Px, Px_min)
        # update to next step
        k_loss_step = (1.0 / (1.0 +np.exp(a_stem*(Px_min-P50))))
        s_t_upper[i] = s_upper_out;
        s_t_lower[i] = s_lower_out;
        s0_upper = s_upper_out
        s0_lower = s_lower_out
        assm_t[i] = (ASM - R) #mol m^-2 d^-1
        px_t[i] = Px
        pxmin_t[i] = Px_min
        k_loss[i]=k_loss_step
        E_t[i] = E
        gs_t[i] = gs*1.6  # save as conductance to H2O
        pl_t[i] = Pl
        Arate_t[i] = max((A_rate - R_sec) * (10 ** 6), 0.0)  # instantaneous net assimilation rate umol/m2-s
        ABA_t[i] = ABA
        Lsum_t[i] = L_sum

    # return variables
    return s_t_upper, s_t_lower, assm_t, px_t, pxmin_t, k_loss, E_t, gs_t, pl_t, Arate_t, ABA_t, Lsum_t, plant_alive, time_of_death

def weather_interp(tRun, measured, times):  # interpolate VPD measurements to fit time of run
    weather_t = np.interp(tRun, times, measured)
    return weather_t

def doy_convert(month, day):  # converts month and day to day of year
    if month == 1:
        doy = day
    if month == 2:
        doy = day + 31
    if month == 3:
        doy = day + 59
    if month == 4:
        doy = day + 90
    if month == 5:
        doy = day + 120
    if month == 6:
        doy = day + 151
    if month == 7:
        doy = day + 181
    if month == 8:
        doy = day + 212
    if month == 9:
        doy = day + 243
    if month == 10:
        doy = day + 273
    if month == 11:
        doy = day + 304
    if month == 12:
        doy = day + 334
    return doy

def find_output_index(obs_time, obs_start_day, dt):
    true_obs_time = obs_time + obs_start_day
    lower_index = np.floor(true_obs_time / dt)
    upper_index = np.ceil(true_obs_time / dt)

    return lower_index.astype(np.int), upper_index.astype(np.int)

def param_fact_setup(a_seq_vals, D0_vals, beta_vals, lamABA_vals, m_vals, Ksat_stem_vals, a_stem_vals, P50_stem_vals,
                     A_root_mult_vals, gs_mult_vals, filename, newfile=True):
    if newfile == False:
        param_fact = np.loadtxt(filename)
        starting_set = len(param_fact)
        first_set = False
        repeat_check = False
    elif newfile == True:
        param_fact = np.zeros(10)
        first_set = True
        repeat_check = False
        starting_set = 0

    for a_seq_i in range(len(a_seq_vals)):
        a_seq = a_seq_vals[a_seq_i]
        for D0_i in range(len(D0_vals)):
            D0 = D0_vals[D0_i]
            for beta_i in range(len(beta_vals)):
                beta = beta_vals[beta_i]
                for lamABA_i in range(len(lamABA_vals)):
                    lamABA = lamABA_vals[lamABA_i]
                    for m_i in range(len(m_vals)):
                        m = m_vals[m_i]
                        for Ksat_i in range(len(Ksat_stem_vals)):
                            Ksat = Ksat_stem_vals[Ksat_i]
                            for a_stem_i in range(len(a_stem_vals)):
                                a_stem = a_stem_vals[a_stem_i]
                                for P50_i in range(len(P50_stem_vals)):
                                    P50 = P50_stem_vals[P50_i]
                                    for A_root_i in range(len(A_root_mult_vals)):
                                        Ar_mult = A_root_mult_vals[A_root_i]
                                        for gs_mult_i in range(len(gs_mult_vals)):
                                            gs_mult = gs_mult_vals[gs_mult_i]
                                            new_set = [a_seq, D0, beta, lamABA, m, Ksat, a_stem, P50, Ar_mult, gs_mult]
                                            if first_set == False:  # check doesn't work for single row array
                                                repeat_check = (param_fact == new_set).all(
                                                    1).any()  # check if set already in factorial
                                            if repeat_check == False or first_set == True:
                                                param_fact = np.vstack((param_fact, new_set))
                                                first_set = False
    if newfile == True:
        param_fact = np.delete(param_fact, 0, 0)  # deleted first row of zeros

    np.savetxt(filename, param_fact, fmt="%5.7f")

    return starting_set, param_fact

def rand_param_set_soil_gen(a_seq_vals, D0_vals, beta_vals, lamABA_vals, m_vals, Ksat_stem_vals, a_stem_vals, P50_stem_vals,
                     A_root_mult_vals, gs_mult_vals, Ksat_soil_vals, Psat_vals, b_vals, sh_vals, new_sets, filename, newfile=True):
    if newfile == False:
        param_fact = np.loadtxt(filename)
        starting_set = len(param_fact)
        first_set = False
        repeat_check = False
    elif newfile == True:
        param_fact = np.zeros(14)
        first_set = True
        repeat_check = False
        starting_set = 0

    for i in range(new_sets):
        a_seq = np.random.uniform(a_seq_vals[0], a_seq_vals[1])
        D0 = np.random.uniform(D0_vals[0], D0_vals[1])
        beta = np.random.uniform(beta_vals[0], beta_vals[1])
        lamABA = np.random.uniform(lamABA_vals[0], lamABA_vals[1])
        m = np.random.uniform(m_vals[0], m_vals[1])
        Ksat = np.random.uniform(Ksat_stem_vals[0], Ksat_stem_vals[1])
        a_stem = np.random.uniform(a_stem_vals[0], a_stem_vals[1])
        P50 = np.random.uniform(P50_stem_vals[0], P50_stem_vals[1])
        Ar_mult = np.random.uniform(A_root_mult_vals[0], A_root_mult_vals[1])
        gs_mult= np.random.uniform(gs_mult_vals[0], gs_mult_vals[1])
        K_soil = np.random.uniform(Ksat_soil_vals[0], Ksat_soil_vals[1])
        Psat = np.random.uniform(Psat_vals[0], Psat_vals[1])
        b = np.random.uniform(b_vals[0], b_vals[1])
        sh = np.random.uniform(sh_vals[0], sh_vals[1])


        new_set = [a_seq, D0, beta, lamABA, m, Ksat, a_stem, P50, Ar_mult, gs_mult, K_soil, Psat, b, sh]
        if first_set == False:  # check doesn't work for single row array
            repeat_check = (param_fact == new_set).all(
                1).any()  # check if set already in factorial
        if repeat_check == False or first_set == True:
            param_fact = np.vstack((param_fact, new_set))
            first_set = False
        elif repeat_check == True:
            i = i - 1;

    if newfile == True:
        param_fact = np.delete(param_fact, 0, 0)  # deleted first row of zeros

    np.savetxt(filename, param_fact, fmt="%5.7f")

    return starting_set, param_fact

def rand_param_set_full_gen(a_seq_vals, D0_vals, beta_vals, lamABA_vals, m_vals, Ksat_stem_vals, a_stem_vals, P50_stem_vals,
                     A_root_mult_vals, gs_mult_vals, Psat_vals, rho_vals, Amax_vals, new_sets, filename, newfile=True):
    if newfile == False:
        param_fact = np.loadtxt(filename)
        starting_set = len(param_fact)
        first_set = False
        repeat_check = False
    elif newfile == True:
        param_fact = np.zeros(16)
        first_set = True
        repeat_check = False
        starting_set = 0

    for i in range(new_sets):
        a_seq = np.random.uniform(a_seq_vals[0], a_seq_vals[1])
        D0 = np.random.uniform(D0_vals[0], D0_vals[1])
        beta = np.random.uniform(beta_vals[0], beta_vals[1])
        lamABA = np.random.uniform(lamABA_vals[0], lamABA_vals[1])
        m = np.random.uniform(m_vals[0], m_vals[1])
        Ksat = np.random.uniform(Ksat_stem_vals[0], Ksat_stem_vals[1])
        a_stem = np.random.uniform(a_stem_vals[0], a_stem_vals[1])
        P50 = np.random.uniform(P50_stem_vals[0], P50_stem_vals[1])
        Ar_mult = np.random.uniform(A_root_mult_vals[0], A_root_mult_vals[1])
        gs_mult= np.random.uniform(gs_mult_vals[0], gs_mult_vals[1])
        Psat = np.random.uniform(Psat_vals[0], Psat_vals[1])
        K_soil= 60.83709996868046*Psat + 0.15420007648383935
        b = 284.5378935340771*Psat + 6.646111280274093
        sh = 0.95*(0.09/0.47) #estimated as 95% of min value measured at site
        rho = np.random.uniform(rho_vals[0], rho_vals[1])
        Amax = np.random.uniform(Amax_vals[0], Amax_vals[1])


        new_set = [a_seq, D0, beta, lamABA, m, Ksat, a_stem, P50, Ar_mult, gs_mult, K_soil, Psat, b, sh, rho, Amax]
        if first_set == False:  # check doesn't work for single row array
            repeat_check = (param_fact == new_set).all(
                1).any()  # check if set already in factorial
        if repeat_check == False or first_set == True:
            param_fact = np.vstack((param_fact, new_set))
            first_set = False
        elif repeat_check == True:
            i = i - 1;

    if newfile == True:
        param_fact = np.delete(param_fact, 0, 0)  # deleted first row of zeros

    np.savetxt(filename, param_fact, fmt="%5.7f")

    return starting_set, param_fact

def run_param(plant, param_factorial, obs_g, obs_g_t, obs_A, obs_A_t, obs_sfinal_upper, obs_sfinal_lower, LWP, LWP_t, depth_re, tRun, dt, s0, VPD_t, Irr_t, T_t, Irr_day, T_day, original_gsmax, set_num, obs_start_day, days_of_obs):
    #update plant with parameters
    plant.a_seq = param_factorial[set_num, 0]
    plant.D0 = param_factorial[set_num, 1]
    plant.beta = param_factorial[set_num, 2]
    plant.canopy.lamABA = param_factorial[set_num, 3]
    plant.canopy.m = param_factorial[set_num, 4]
    plant.stem.Ksat_stem = param_factorial[set_num, 5]
    plant.stem.a_stem = param_factorial[set_num, 6]
    plant.stem.P50_stem = param_factorial[set_num, 7]
    plant.soil_root.A_root_mult = param_factorial[set_num, 8]
    plant.canopy.A_root_mult = param_factorial[set_num, 8]
    plant.canopy.Gs_leaf = param_factorial[set_num, 9] * original_gsmax

    #run simulation
    s_t_upper, s_t_lower, assm_t, px_t, pxmin_t, _, E_t, gs_t, pl_t, Arate_t, ABA_t, L_sum, _, time_of_death = simulate_s_t_norefilling_varVPD(
        depth_re, tRun, dt, s0, s0, plant, VPD_t, Irr_t, T_t, Irr_day, T_day)

    #get scores
    g_diff = np.zeros_like(obs_g)
    for i in range(len(obs_g)):
        obs_time = obs_g_t[i]
        lower_index, upper_index = find_output_index(obs_time, obs_start_day, dt)
        model_avg = (gs_t[lower_index] + gs_t[upper_index]) / 2.0
        g_diff[i] = abs(obs_g[i] - model_avg) / obs_g[i]
    mean_g_diff = np.mean(g_diff)

    lwp_diff = np.zeros_like(LWP)
    for i in range(len(LWP)):
        obs_time = LWP_t[i]
        lower_index, upper_index = find_output_index(obs_time, obs_start_day, dt)
        model_avg = (pl_t[lower_index] + pl_t[upper_index]) / 2.0
        lwp_diff[i] = abs(LWP[i] - model_avg) / abs(LWP[i])
    mean_lwp_diff = np.mean(lwp_diff)

    soil_start_index, _ = find_output_index(0.0, obs_start_day, dt)
    _, soil_end_index = find_output_index(days_of_obs, obs_start_day, dt)
    model_soil_mean_upper = np.mean(s_t_upper[soil_start_index:soil_end_index])
    model_soil_mean_lower = np.mean(s_t_lower[soil_start_index:soil_end_index])
    soil_diff = ((abs(model_soil_mean_upper - obs_sfinal_upper) / obs_sfinal_upper) + (abs(model_soil_mean_lower - obs_sfinal_lower) / obs_sfinal_lower))/2.0

    total_score = (mean_g_diff + mean_lwp_diff + soil_diff) / 3.0

    return set_num, param_factorial[set_num, 0], param_factorial[set_num, 1], param_factorial[set_num, 2], param_factorial[set_num, 3], param_factorial[set_num, 4], param_factorial[set_num, 5], param_factorial[set_num, 6], param_factorial[set_num, 7], param_factorial[set_num, 8], param_factorial[set_num, 9], mean_g_diff, mean_lwp_diff, soil_diff, total_score, time_of_death

def run_param_tuned_soil(plant, param_factorial, obs_g, obs_g_t, obs_A, obs_A_t, obs_sfinal_upper, obs_sfinal_lower, LWP, LWP_t, depth_re, tRun, dt, s0, VPD_t, Irr_t, T_t, Irr_day, T_day, original_gsmax, set_num, obs_start_day, days_of_obs):
    #update plant with parameters
    plant.a_seq = param_factorial[set_num, 0]
    plant.D0 = param_factorial[set_num, 1]
    plant.beta = param_factorial[set_num, 2]
    plant.canopy.lamABA = param_factorial[set_num, 3]
    plant.canopy.m = param_factorial[set_num, 4]
    plant.stem.Ksat_stem = param_factorial[set_num, 5]
    plant.stem.a_stem = param_factorial[set_num, 6]
    plant.stem.P50_stem = param_factorial[set_num, 7]
    plant.soil_root.A_root_mult = param_factorial[set_num, 8]
    plant.canopy.A_root_mult = param_factorial[set_num, 8]
    plant.canopy.Gs_leaf = param_factorial[set_num, 9] * original_gsmax
    plant.soil.Ksat_soil = param_factorial[set_num, 10]
    plant.soil.Psat_soil = param_factorial[set_num, 11]
    plant.soil.b = param_factorial[set_num, 12]
    plant.soil.sh = param_factorial[set_num, 13]

    parameter_set_soil = {'Ps': param_factorial[set_num, 11], 'b':param_factorial[set_num, 12], 'n':0.47, 'Ksat': param_factorial[set_num, 10], 'sh':param_factorial[set_num, 13], 'sfc':0.6, 'sst':0.57}
    soil_dict['tuned'] = parameter_set_soil

    #run simulation
    s_t_upper, s_t_lower, assm_t, px_t, pxmin_t, _, E_t, gs_t, pl_t, Arate_t, ABA_t, L_sum, _, time_of_death = simulate_s_t_norefilling_varVPD(
        depth_re, tRun, dt, s0, s0, plant, VPD_t, Irr_t, T_t, Irr_day, T_day)

    #get scores
    g_diff = np.zeros_like(obs_g)
    for i in range(len(obs_g)):
        obs_time = obs_g_t[i]
        lower_index, upper_index = find_output_index(obs_time, obs_start_day, dt)
        model_avg = (gs_t[lower_index] + gs_t[upper_index]) / 2.0
        g_diff[i] = abs(obs_g[i] - model_avg) / obs_g[i]
    mean_g_diff = np.mean(g_diff)

    lwp_diff = np.zeros_like(LWP)
    lwp_md_diff = np.zeros_like(LWP)
    for i in range(len(LWP)):
        obs_time = LWP_t[i]
        if obs_time % 1.0 > 0.5: #checks if LWP is a midday msmt
            midday=True
        else:
            midday =False
        lower_index, upper_index = find_output_index(obs_time, obs_start_day, dt)
        model_avg = (pl_t[lower_index] + pl_t[upper_index]) / 2.0
        lwp_diff[i] = abs(LWP[i] - model_avg) / abs(LWP[i])
        if midday == True:
            lwp_md_diff[i] = abs(LWP[i] - model_avg) / abs(LWP[i])
        else:
            lwp_md_diff[i] = np.nan
    mean_lwp_diff = np.mean(lwp_diff)
    mean_md_lwp_diff = np.nanmean(lwp_md_diff)

    soil_start_index, _ = find_output_index(0.0, obs_start_day, dt)
    _, soil_end_index = find_output_index(days_of_obs, obs_start_day, dt)
    model_soil_mean_upper = np.mean(s_t_upper[soil_start_index:soil_end_index])
    model_soil_mean_lower = np.mean(s_t_lower[soil_start_index:soil_end_index])
    soil_diff = ((abs(model_soil_mean_upper - obs_sfinal_upper) / obs_sfinal_upper) + (abs(model_soil_mean_lower - obs_sfinal_lower) / obs_sfinal_lower))/2.0

    total_score = (mean_g_diff + mean_md_lwp_diff + soil_diff) / 3.0

    return set_num, param_factorial[set_num, 0], param_factorial[set_num, 1], param_factorial[set_num, 2], param_factorial[set_num, 3], param_factorial[set_num, 4], param_factorial[set_num, 5], param_factorial[set_num, 6], param_factorial[set_num, 7], param_factorial[set_num, 8], param_factorial[set_num, 9], param_factorial[set_num, 10], param_factorial[set_num, 11], param_factorial[set_num, 12], param_factorial[set_num, 13], mean_g_diff, mean_lwp_diff, mean_md_lwp_diff, soil_diff, total_score, time_of_death

def rand_historical_rainfall(n_trajectories, tmax, dt, start_doy, n, Zr):
    dfN = pd.read_csv('Data_Downloads/Ninghan_Station.csv')
    tRun = np.arange(0, tmax, dt)
    dfN['Date'] = pd.to_datetime(dfN[['Year', 'Month', 'Day']])
    upper_dates=dfN.Date.size-(2.0*tmax)
    depth_re=np.zeros((n_trajectories, len(tRun)))
    start_year = np.zeros(n_trajectories)
    traj_num=0
    while traj_num < n_trajectories:
        start_seed = np.random.randint(0, upper_dates)
        start_match=0
        scan_index=0
        year = 0.0
        while start_match==0:
            check_index=start_seed+scan_index
            rain_doy=doy_convert(dfN.Month[check_index], dfN.Day[check_index])
            if rain_doy == start_doy:
                start_match=1
                year = dfN.Year[check_index]
            else:
                scan_index += 1
        next_day=0
        missing_val=0
        for i in range(len(tRun)):
            current_day=np.floor(tRun[i])
            if current_day == next_day:
                rain_index=check_index + next_day
                meas_rain=dfN.Rainfall[rain_index]
                if meas_rain ==0 or meas_rain >0:
                    depth_re[traj_num, i]=meas_rain/(n*Zr)/1000.0
                    next_day=next_day + 1
                else: #missing rainfall value
                    missing_val=1
                    break
            else:
                depth_re[traj_num, i]=0.0
        if missing_val == 0:
            start_year[traj_num] = year
            traj_num = traj_num + 1
    return depth_re, start_year

def historical_rainfall(tmax, dt, start_doy, n, Zr):
    #finds all complete trajectories with a given starting date and length
    n_max_trajectories=150
    dfN = pd.read_csv('Data_Downloads/Ninghan_Station.csv')
    tRun = np.arange(0,tmax+dt,dt)
    dfN['Date'] = pd.to_datetime(dfN[['Year', 'Month', 'Day']])
    upper_dates=dfN.Date.size-(2.0*tmax)
    depth_re=np.zeros((n_max_trajectories, len(tRun)))
    start_year = np.zeros(n_max_trajectories)
    traj_num=0

    for i in range(len(dfN.Year)-(tmax)):
        rain_doy = doy_convert(dfN.Month[i], dfN.Day[i])
        if rain_doy == start_doy:
            year = dfN.Year[i]
            next_day = 0
            missing_val = 0
            for j in range(len(tRun)):
                current_day = np.floor(tRun[j])
                if current_day == next_day:
                    rain_index = i + next_day
                    meas_rain = dfN.Rainfall[rain_index]
                    if meas_rain == 0 or meas_rain > 0:
                        depth_re[traj_num, j] = meas_rain / (n * Zr) / 1000.0
                        next_day = next_day + 1
                    else:  # missing rainfall value
                        missing_val = 1
                        break
                else:
                    depth_re[traj_num, j] = 0.0
            if missing_val == 0: #no missing values, keeps trajectory
                start_year[traj_num] = year
                traj_num = traj_num + 1
    depth_re_out=depth_re[0:(traj_num), :] #get rid of empty rows
    start_year_out = start_year[0:(traj_num)]
    return depth_re_out, start_year_out

def run_full_param(plant, param_factorial, obs_g, obs_g_t, obs_A, obs_A_t, obs_sfinal_upper, obs_sfinal_lower, LWP, LWP_t, depth_re, tRun, dt, s0, VPD_t, Irr_t, T_t, Irr_day, T_day, original_gsmax, set_num, obs_start_day, days_of_obs):
    #update plant with parameters
    plant.a_seq = param_factorial[set_num, 0]
    plant.D0 = param_factorial[set_num, 1]
    plant.beta = param_factorial[set_num, 2]
    plant.canopy.lamABA = param_factorial[set_num, 3]
    plant.canopy.m = param_factorial[set_num, 4]
    plant.stem.Ksat_stem = param_factorial[set_num, 5]
    plant.stem.a_stem = param_factorial[set_num, 6]
    plant.stem.P50_stem = param_factorial[set_num, 7]
    plant.soil_root.A_root_mult = param_factorial[set_num, 8]
    plant.canopy.Gs_leaf = param_factorial[set_num, 9] * original_gsmax
    plant.soil.Ksat_soil = param_factorial[set_num, 10]
    plant.soil.Psat_soil = param_factorial[set_num, 11]
    plant.soil.b = param_factorial[set_num, 12]
    plant.soil.sh = param_factorial[set_num, 13]
    plant.canopy.rho = param_factorial[set_num, 14]
    plant.canopy.Amax = param_factorial[set_num, 15]

    parameter_set_soil = {'Ps': param_factorial[set_num, 11], 'b': param_factorial[set_num, 12], 'n': 0.47,
                          'Ksat': param_factorial[set_num, 10], 'sh': param_factorial[set_num, 13], 'sfc': 0.6,
                          'sst': 0.57}
    soil_dict['tuned'] = parameter_set_soil

    #run simulation
    s_t_upper, s_t_lower, assm_t, px_t, pxmin_t, _, E_t, gs_t, pl_t, Arate_t, ABA_t, L_sum, _, time_of_death = simulate_s_t_norefilling_varVPD(
        depth_re, tRun, dt, s0, s0, plant, VPD_t, Irr_t, T_t, Irr_day, T_day)

    #get scores
    g_diff = np.zeros_like(obs_g)
    for i in range(len(obs_g)):
        obs_time = obs_g_t[i]
        lower_index, upper_index = find_output_index(obs_time, obs_start_day, dt)
        model_avg = (gs_t[lower_index] + gs_t[upper_index]) / 2.0
        g_diff[i] = abs(obs_g[i] - model_avg) / obs_g[i]
    mean_g_diff = np.mean(g_diff)

    A_diff = np.zeros_like(obs_A)
    for i in range(len(obs_A)):
        obs_time = obs_A_t[i]
        lower_index, upper_index = find_output_index(obs_time, obs_start_day, dt)
        model_avg = (Arate_t[lower_index] + Arate_t[upper_index]) / 2.0
        A_diff[i] = abs(obs_A[i] - model_avg) / obs_A[i]
    mean_A_diff = np.mean(A_diff)

    lwp_diff = np.zeros_like(LWP)
    for i in range(len(LWP)):
        obs_time = LWP_t[i]
        lower_index, upper_index = find_output_index(obs_time, obs_start_day, dt)
        model_avg = (pl_t[lower_index] + pl_t[upper_index]) / 2.0
        lwp_diff[i] = abs(LWP[i] - model_avg) / abs(LWP[i])
    mean_lwp_diff = np.mean(lwp_diff)

    soil_start_index, _ = find_output_index(0.0, obs_start_day, dt)
    _, soil_end_index = find_output_index(days_of_obs, obs_start_day, dt)
    model_soil_mean_upper = np.mean(s_t_upper[soil_start_index:soil_end_index])
    model_soil_mean_lower = np.mean(s_t_lower[soil_start_index:soil_end_index])
    soil_diff = ((abs(model_soil_mean_upper - obs_sfinal_upper) / obs_sfinal_upper) + (abs(model_soil_mean_lower - obs_sfinal_lower) / obs_sfinal_lower))/2.0

    total_score = (mean_g_diff + mean_lwp_diff + soil_diff + mean_A_diff) / 4.0

    return set_num, param_factorial[set_num, 0], param_factorial[set_num, 1], param_factorial[set_num, 2], \
           param_factorial[set_num, 3], param_factorial[set_num, 4], param_factorial[set_num, 5], param_factorial[
               set_num, 6], param_factorial[set_num, 7], param_factorial[set_num, 8], param_factorial[set_num, 9], \
           param_factorial[set_num, 10], param_factorial[set_num, 11], param_factorial[set_num, 12],param_factorial[
               set_num, 13], param_factorial[set_num, 14],param_factorial[set_num, 15], mean_g_diff, mean_lwp_diff,\
           soil_diff, mean_A_diff, total_score, time_of_death

def run_trajectory(plant, param_factorial,rainfall_data, rainfall_dates, tRun, dt, s0_upper, s0_lower, VPD_t, Irr_t, T_t, Irr_day, T_day, original_gsmax, param_set_num, rain_year, main_folder):
    save_shorter_data = True #if false saves hourly resolution data, otherwise daily

    rain_date = rainfall_dates[rain_year]

    directory = "year_{0}_param_{1}".format(int(rain_date), int(param_set_num))

    #update plant with parameters
    plant.a_seq = param_factorial[param_set_num, 0]
    plant.D0 = param_factorial[param_set_num, 1]
    plant.beta = param_factorial[param_set_num, 2]
    plant.canopy.lamABA = param_factorial[param_set_num, 3]
    plant.canopy.m = param_factorial[param_set_num, 4]
    plant.stem.Ksat_stem = param_factorial[param_set_num, 5]
    plant.stem.a_stem = param_factorial[param_set_num, 6]
    plant.stem.P50_stem = param_factorial[param_set_num, 7]
    plant.soil_root.A_root_mult = param_factorial[param_set_num, 8]
    plant.canopy.Gs_leaf = param_factorial[param_set_num, 9] * original_gsmax
    plant.soil.Ksat_soil = param_factorial[param_set_num, 10]
    plant.soil.Psat_soil = param_factorial[param_set_num, 11]
    plant.soil.b = param_factorial[param_set_num, 12]
    plant.soil.sh = param_factorial[param_set_num, 13]
    plant.canopy.rho = param_factorial[param_set_num, 14]
    plant.canopy.Amax = param_factorial[param_set_num, 15]

    parameter_set_soil = {'Ps': param_factorial[param_set_num, 11], 'b': param_factorial[param_set_num, 12], 'n': 0.47,
                          'Ksat': param_factorial[param_set_num, 10], 'sh': param_factorial[param_set_num, 13], 'sfc': 0.6,
                          'sst': 0.57}
    soil_dict['tuned'] = parameter_set_soil

    depth_re = rainfall_data[rain_year]

    #run simulation
    s_t_upper, s_t_lower, assm_t, px_t, pxmin_t, k_loss, E_t, gs_t, pl_t, Arate_t, ABA_t, L_sum, _, time_of_death = simulate_s_t_norefilling_varVPD(
        depth_re, tRun, dt, s0_upper, s0_lower, plant, VPD_t, Irr_t, T_t, Irr_day, T_day)

    path = os.path.join(main_folder, directory)
    os.makedirs(path)

    if save_shorter_data == False:

        np.savetxt(os.path.join(path, "s_t_upper.txt"), s_t_upper, delimiter=',')
        np.savetxt(os.path.join(path, "s_t_lower.txt"), s_t_lower, delimiter=',')
        np.savetxt(os.path.join(path, "assm_t.txt"), assm_t, delimiter=',')
        np.savetxt(os.path.join(path, "px_t.txt"), px_t, delimiter=',')
        np.savetxt(os.path.join(path, "pxmin_t.txt"), pxmin_t, delimiter=',')
        np.savetxt(os.path.join(path, "k_loss.txt"), k_loss, delimiter=',')
        np.savetxt(os.path.join(path, "E_t.txt"), E_t, delimiter=',')
        np.savetxt(os.path.join(path, "gs_t.txt"), gs_t, delimiter=',')
        np.savetxt(os.path.join(path, "pl_t.txt"), pl_t, delimiter=',')
        np.savetxt(os.path.join(path, "Arate_t.txt"), Arate_t, delimiter=',')
        np.savetxt(os.path.join(path, "L_sum.txt"), L_sum, delimiter=',')

    else:
        #saves these variables at midnight each day
        slice_size = int(1.0/dt)
        s_t_upper_daily = s_t_upper[::slice_size]
        s_t_lower_daily = s_t_lower[::slice_size]
        pxmin_daily = pxmin_t[::slice_size]
        Lsum_daily = L_sum[::slice_size]
        k_loss_daily = k_loss[::slice_size]

        np.savetxt(os.path.join(path, "s_t_upper_daily.txt"), s_t_upper_daily, delimiter=',')
        np.savetxt(os.path.join(path, "s_t_lower_daily.txt"), s_t_lower_daily, delimiter=',')
        np.savetxt(os.path.join(path, "running_pxmin_daily.txt"), pxmin_daily, delimiter=',')
        np.savetxt(os.path.join(path, "Lsum_daily.txt"), Lsum_daily, delimiter=',')
        np.savetxt(os.path.join(path, "k_loss_daily.txt"), k_loss_daily, delimiter=',')

        #saves these as daily totals
        E_t_sub = E_t[:-1]
        assm_t_sub = assm_t[:-1]
        daily_Px_min_sub = px_t[:-1]

        E_daily = dt*np.sum(E_t_sub.reshape(-1,slice_size), axis=1) #total daily transpiration [m3]
        assm_daily = dt*np.sum(assm_t_sub.reshape(-1,slice_size), axis=1) #total daily assimilation [mol m^-2]
        daily_Px_min = np.min(daily_Px_min_sub.reshape(-1, slice_size), axis=1)  # daily minimum water potential [MPa]

        np.savetxt(os.path.join(path, "E_daily_total.txt"), E_daily, delimiter=',')
        np.savetxt(os.path.join(path, "assm_daily.txt"), assm_daily, delimiter=',')
        np.savetxt(os.path.join(path, "daily_Px_min.txt"), daily_Px_min, delimiter=',')

    return param_set_num, rain_date, time_of_death

def nan_helper(y):
    """Helper to handle indices and logical indices of NaNs.

    Input:
        - y, 1d numpy array with possible NaNs
    Output:
        - nans, logical indices of NaNs
        - index, a function, with signature indices= index(logical_indices),
          to convert logical indices of NaNs to 'equivalent' indices
    Example:
        >>> # linear interpolation of NaNs
        >>> nans, x= nan_helper(y)
        >>> y[nans]= np.interp(x(nans), x(~nans), y[~nans])
    """

    return np.isnan(y), lambda z: z.nonzero()[0]

