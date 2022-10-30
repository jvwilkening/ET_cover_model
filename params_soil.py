"""

Ps:saturated soil water potential; MPa
    b: exponent on the water retention curve
    n: porosity
    Ksat: saturated soil hydraulic conductivity; m/d;

For the cover soil, the hygroscopic and field capacity values were estimated from the minimum and maximum measured soil
water contents, respectively, during the course of the experiment
http://www.nrcs.usda.gov/wps/portal/nrcs/detail/soils/survey/office/ssr10/tr/?cid=nrcs144p2_074846
Parameters for loamy_sandy derived from measurements at field sites in Morongo Valley

cover2 comes from fitting soil water retention curve to measurements of soil moisture and predawn leaf water potentials

"""
soil_dict = {'sand':        {'Ps': -0.34*10**(-3), 'b':4.05, 'n':0.35, 'Ksat':2.0, 'sh':0.08, 'sfc':0.35, 'sst':0.33},
             'loamy_sand':  {'Ps': -0.17*10**(-3), 'b':4.38, 'n':0.42, 'Ksat':1.0, 'sh':0.08, 'sfc':0.52, 'sst':0.31},
             'sandy_loam':  {'Ps': -0.70*10**(-3), 'b':4.90, 'n':0.43, 'Ksat':0.80, 'sh':0.14, 'sfc':0.56, 'sst':0.46},
             'loam':        {'Ps': -1.43*10**(-3), 'b':5.39, 'n':0.45, 'Ksat':0.20, 'sh':0.19, 'sfc':0.65, 'sst':0.57},
             'clay':        {'Ps': -1.82*10**(-3), 'b':11.4, 'n':0.50, 'Ksat':0.04, 'sh':0.47, 'sfc':1.0, 'sst':0.78},
             'cover':       {'Ps': -1.345*10**(-2), 'b':3.42, 'n':0.47, 'Ksat':0.2, 'sh':0.12, 'sfc':0.65, 'sst':0.57},
             'cover2':      {'Ps': -1.095333, 'b':0.74689, 'n':0.47, 'Ksat':0.2, 'sh':0.2, 'sfc':0.6, 'sst':0.57},#using this soil Jan 2021
             'tuned':      {'Ps': -1.095333, 'b':0.74689, 'n':0.47, 'Ksat':0.2, 'sh':0.2, 'sfc':0.6, 'sst':0.57}, #overwrites values with parameter set
             'cover3':      {'Ps': -0.002, 'b':7.0, 'n':0.47, 'Ksat':0.2, 'sh':0.14, 'sfc':0.56, 'sst':0.46}}

