#%%
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
from scipy.signal import argrelextrema
from scipy.integrate import trapz
from scipy.optimize import Bounds, minimize

import meteostat

from geopy.geocoders import Nominatim

import requests
import json

import weatherapi

import plotly.io as pio


import warnings
pio.templates.default = "plotly"
warnings.filterwarnings('ignore')



#%%
kelvin = lambda x: x+273.15 # convert Celsius to Kelvin
celsius = lambda x: x-273.15 # convert Kelvin to Celsius
kWh = lambda x: x / 3600000 # convert Joule to kWh
Joule = lambda x: x * 3600000 # convert kWh to Joule

#freezing_zone = pd.Interval(left=-6, right=4) # temp between hp starts freezing and has to be thawed
freezing_zone = 4
freezing_factor = 1.2# factor which the efficiency coefficient is divided by in case of being in freezing zone
flow_max = 50
aussen_min = -10
room_temp = 21 # temperature at which room should be heated
ruecklauf = 25
ny = 0.4 #Gütegrad

temp_start_heat = 15 # temperature under which house has to be heated
house_energy = 200 # Watt / Kelvin

timeframe_h = 24 #hours
timeframe_s = timeframe_h * 3600 #seconds

#take Newton cooling data

#20 and 60 after some DIN-Norms
t_env = kelvin(20) #environmental temperature where water storage is located
t_0 = kelvin(60) #temperature at which water storage power is measured
#coef_h = 0.000625

#specifications of water storage (taken from real example)
storage_power = 200 #Watts
h2o_energy = 4200 #Joule / (kilogram * Kelvin)
storage_volume = 5000 #Liter = kilogram -> to be defined later on

#specifications of heatpump
hp_power = 7000 #Watts
hp_power_opt = 5000 #Watts -> Power in which heatpump works most efficiently (can be taken into calculation later)
hp_power_el_opt = 15000 #Watts -> Power in which heatpump works most efficiently (can be taken into calculation later)
max_heat_cap = 75 # temperature, up to which the heat pump can maximally heat

# other measures
groundwater = celsius(t_env)

flowtemp = lambda t: - (flow_max - room_temp) / (room_temp-aussen_min) *t + flow_max + aussen_min * (flow_max - room_temp) / (room_temp-aussen_min)
needed_house_energy = lambda t: np.maximum(0, room_temp - t) * house_energy # returns needed energy to heat house in Watt

#eta = lambda x: np.maximum(0, np.minimum(300, ny * kelvin(flowtemp(x)) / (kelvin(flowtemp(x))-kelvin(x))))
h = True
def eta(zu, ab = None): #theoretischer Wirkungsgrad
    if ab is None:
        ab = flowtemp(zu)
    return np.maximum(0, np.minimum(300, ny * kelvin(ab) / (kelvin(ab)-kelvin(zu))))

def eta_heating(zu, max_ab = None, min_ab = None ): #"gleitender" Wirkungsgrad -> da sich über temperaturerhöhung ändert
    try:
        if type(min_ab) == type(None):
            min_ab = ruecklauf
        if type(max_ab) == type(None):
            max_ab = flowtemp(zu)
    except:
        print("Error: function 'eta_heating'")
    return ny * (kelvin(zu) * np.log((kelvin(max_ab) - kelvin(zu)) / (kelvin(min_ab) - kelvin(zu))) + kelvin(max_ab) - kelvin(min_ab)) / (kelvin(max_ab) - kelvin(min_ab))

def eta_heating2(zu, max_ab = None, min_ab=None):
    try:
        if type(min_ab) == type(None):
            min_ab = ruecklauf
        if type(max_ab) == type(None):
            max_ab = flowtemp(zu)
    except:
        print("Error: function 'eta_heating'")

    return ny * ((kelvin(max_ab) - kelvin(min_ab)) / (kelvin(max_ab) - kelvin(min_ab) - kelvin(zu) * np.log(kelvin(max_ab) / kelvin(min_ab))))

def eta_freezing_adj(zu, max_ab = None, min_ab = None, heating =  False): # if heating = True, then use eta_heating
    try:
        if type(min_ab) == type(None):
            min_ab = ruecklauf
        if type(max_ab) == type(None):
            max_ab = flowtemp(zu)
    except:
        print("Error: function 'eta_heating'")

    if heating == True:
        e = eta_heating2(zu, max_ab, min_ab)
    else:
        e = eta(zu, max_ab)
    #print(type(eta))
    if type(e) == np.float64:
        if zu < 4:
            e = e / freezing_factor
    else:
        e[zu < 4] = e[zu < 4] / 1.2

    return e

def newton_cooling(t):
    return t_env + (t_0 - t_env) * np.exp(-coef_h * t) # = temperature at time t

def newton_cooling_diff(t):
    return coef_h * (t_env - newton_cooling(t)) # = temperature loss rate at time t

def newton_cooling_t_0(temp_t, t=7*3600): #enter in Kelvin
    #return ((temp_t - t_env) * np.exp((coef_h * t).T)).T + t_env #returns t_0 for fix temperature after time t
    return (temp_t - t_env) * np.exp((coef_h * t)) + t_env #returns t_0 for fix temperature after time t

def newton_cooling_t(temp_t):#, coef_h = coef_h):
    return -np.log((temp_t - t_env) / (t_0 - t_env)) / coef_h #returns how much time it takes for t_0 to cool down to wished temperature temp_t

def storage_power_temp(temp):
    return  -h2o_energy * storage_volume * newton_cooling_diff(newton_cooling_t(temp))#returns how much power in Watt is needed to keep temperature at temp level

def get_d_t():
    return storage_power / (h2o_energy * storage_volume) # returns heat loss in kelvin per second

def get_coef_h(d_t):
    return d_t / np.abs(t_env - t_0) #returns coefficient at which heat decreases in storage

def get_coef_h_adj(temp, flow, t):
    return - np.log((flow - t_env) / (temp - t_env)) / t # returns adjusted coef_h for given flow, overheat_temp and time after which flow should be reached

def storage_power_temp_adj(temp, flow, t):
    return h2o_energy * storage_volume * (coef_h - get_coef_h_adj(temp, flow, t)) * (temp - t_env) # returns needed power to keep temperature in storage and decrease at given rate

def get_needed_storage_vol(temp_diff, kWh):
    return (kWh * 3600000 )/ (temp_diff * h2o_energy) #returns waterstorage volume to store #kWh with given temp difference

d_t = get_d_t()
coef_h = get_coef_h(d_t)

integral = lambda x, y: trapz(y, x) # returns the estimated integral for data points x and y
integral_avg = lambda x, y: integral(x, y) / (x.shape[0] - 1) # returns the average of the integral for points x and y

heat = lambda ground_temp, flowtemp, eta_t_0: (flowtemp - ground_temp) * h2o_energy * storage_volume / eta_t_0 # energy needed to heat water from ground temp to flow with given eta (in Joule)
overheat = lambda time, flowtemp, eta_t_0: (newton_cooling_t_0(kelvin(flowtemp), t = time*3600) - 273.15 - flowtemp) * h2o_energy * storage_volume / eta_t_0 # energy needed to heat water storage from flowtemp to t_0
keepheat = lambda time, flowtemp, eta_avg: storage_power_temp(kelvin(flowtemp))*3600*time / eta_avg # energy needed to keep heat in water storage on flowtemp with changing eta
def over_keep_mix(time, flowtemp, eta_t_0, eta_avg): # get the minimum of a mix of overheating and keeping the heat

    try:
        if time == 0: #if time equals 0, then there is no optimum that can be computed
            return flowtemp, heat(groundwater, flowtemp, eta_t_0), 0

        def target_function(temp): #returns the energy needed for heating, overheating and keeping the heat level in Joule
            over = (temp - flowtemp) * h2o_energy * storage_volume / eta_t_0
            keep = storage_power_temp_adj(kelvin(temp), kelvin(flowtemp), time * 3600) * 3600 * time / eta_avg
            h = heat(groundwater, flowtemp, eta_t_0)

            return over + keep + h

        lb = flowtemp
        ub = max(flowtemp, (celsius(newton_cooling_t_0(kelvin(flowtemp), time*3600))))

        cons = Bounds(lb = lb, ub= ub)
        #print(cons)
        min = minimize(fun = target_function, x0=(flowtemp+0.0001, ), bounds=cons)

        return min.x[0], min.fun[0], storage_power_temp_adj(kelvin(min.x[0]), kelvin(flowtemp), time * 3600) #TODO insert what is returned
    except:
        return np.nan, np.nan, np.nan

over_keep_mix_vectorized = np.vectorize(over_keep_mix)

energy_total_overheat = lambda ground_temp, flowtemp, eta_t_0, time: heat(ground_temp, flowtemp, eta_t_0) + overheat(time, flowtemp, eta_t_0) # just consider overheating
energy_total_keepheat = lambda ground_temp, flowtemp, eta_t_0, eta_avg, time: heat(ground_temp, flowtemp, eta_t_0) + keepheat(time, flowtemp, eta_avg) # just consider keeping heat



#%%
def add_data(df):
    df['efficiency'] = eta_freezing_adj(df.temperature, heating=h)
    df['flow_temp'] = flowtemp(df.temperature)
    df['energy_needed_water'] = needed_house_energy(df.temperature)
    df['energy_needed_elec'] = needed_house_energy(df.temperature) / df.efficiency
    df['time'] = df.date.dt.time
    df['freezing_zone'] = df['temperature'] < freezing_zone

    return df

def get_local_extrema(df): #get local minima and maxima in function
    check_interval = 10

    df['local_min'] = df.iloc[argrelextrema(df.temperature.values, np.less_equal, order=check_interval)[0]]['temperature']
    df['local_max'] = df.iloc[argrelextrema(df.temperature.values, np.greater_equal, order=check_interval)[0]]['temperature']
    df['local_max_flow'] = df.iloc[argrelextrema(df.flow_temp.values, np.greater_equal, order=check_interval)[0]]['flow_temp']

    #delete doubled high and low points -> take only the later ones
    df.local_min[(pd.notna(df.local_min - df.local_min.shift(-1))) | (pd.notna(df.local_min - df.local_min.shift(-2))) | (pd.notna(df.local_min - df.local_min.shift(-3)))] = np.nan
    df.local_max[(pd.notna(df.local_max - df.local_max.shift(1))) | (pd.notna(df.local_max - df.local_max.shift(2)))| (pd.notna(df.local_max - df.local_max.shift(3)))] = np.nan

    # high and low points have to alternate
    local_max_index = df[pd.notna(df.local_max)].index # get indexes of local max
    local_min_index = df[pd.notna(df.local_min)].index # get indexes of local min

    index = np.sort(np.concatenate((local_min_index.values, local_max_index.values))) # concat and sort both indexes
    index_pd = pd.Series(np.isin(index, local_min_index))

    index = index[(index_pd - index_pd.shift(1)) != 0] # when two max or two min indexes occur, shifting and taking the difference results in 0

    local_min_index = local_min_index[np.isin(local_min_index, index)] # take only remaining indexes
    local_max_index = local_max_index[np.isin(local_max_index, index)]

    if local_max_index[0] < local_min_index[0]: #reassure that period counting starts with period_min -> heating period comes before deheating
        local_max_index = local_max_index[1:]

    for i, index in enumerate(local_min_index):
        df.loc[df.index > index, 'period_min'] = i

    for i, index in enumerate(local_max_index):
        df.loc[df.index >= index, 'period_max'] = i

    return df


#%%
def read_dat(path):
    data = pd.read_table(path, skiprows=30)
    col = data.columns[0]
    data = data[col].str.split(expand=True)
    data.columns = col.split()
    data = data.loc[1:]

    data.MM = data.MM.str.zfill(2)
    data.DD = data.DD.str.zfill(2)
    data.HH = data.HH.str.zfill(2)
    index_24 = data.HH[data.HH == '24'].index
    data.HH[data.HH == '24'] = '23'
    data['date'] = pd.to_datetime('2022-' + data.MM.astype(str) + '-' + data.DD.astype(str) + ' ' + data.HH.astype(str) + ':00', format = '%Y-%m-%d %H:%M')

    data.loc[index_24, 'date'] = data.loc[index_24].date + pd.to_timedelta(1, 'h')
    data.rename({'t': 'temperature'}, axis=1, inplace=True)
    data = data[['date','temperature']]
    data.temperature = data.temperature.astype(float)
    data = data.sort_values(by='date')
    return data

data_temp = read_dat('data/TRY_data/TRY2015_490148084221_Jahr.dat')
data_temp = data_temp[(data_temp.date >= pd.to_datetime('2022-01-01')) & (data_temp.date <= pd.to_datetime('2022-04-30'))]
#data_temp = data_temp[(data_temp.date.dt.month >= 10) | (data_temp.date.dt.month <= 3)]
data_temp.reset_index(drop = True, inplace=True)


#%%


#data_temp = pd.read_table('data/stundenwerte_TU_VS_akt/produkt_tu_stunde_20201216_20220618_05229.txt', sep=';') #Villingen-Schwenningen
data_temp = pd.read_table('data/weatherstationdata/produkt_tu_stunde_19480101_20081031_02522.txt', sep=';') #Karlsruhe
#data_temp = pd.read_table('data/weatherstationdata/produkt_tu_stunde_19500101_20211231_05792.txt', sep=';') #Zugspitze
#data_temp = pd.read_table('data/stundenwerte_TU_05100_19410101_20211231_hist/produkt_tu_stunde_19410101_20211231_05100.txt', sep=';') #Trier?
data_temp = data_temp.drop(['eor', 'QN_9', 'RF_TU'], axis=1) # drop end-of-row indicator
data_temp.rename(columns = {'TT_TU':'temperature', 'MESS_DATUM': 'date', }, inplace = True)
data_temp.date = pd.to_datetime(data_temp.date.astype(str), format='%Y%m%d%H')
data_temp[data_temp == -999] = np.nan # rename missing values to numpy nan
data_temp = data_temp[(data_temp.date >= pd.to_datetime('2007-10-01')) & (data_temp.date <= pd.to_datetime('2008-04-30'))]
#data_temp = data_temp[(data_temp.date.dt.month >= 10) | (data_temp.date.dt.month <= 3)]
data_temp.reset_index(drop = True, inplace=True)


#%%
data_temp = add_data(data_temp)
data_temp = get_local_extrema(data_temp)