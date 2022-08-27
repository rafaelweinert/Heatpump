import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
from scipy.signal import argrelextrema
from scipy.integrate import trapz
from scipy.optimize import Bounds, minimize

import meteostat
import geopy
from geopy.geocoders import Nominatim

import requests
import json

import weatherapi

import plotly.io as pio

import warnings

pio.templates.default = "plotly"
warnings.filterwarnings('ignore')

"""# Define variables and formulas"""


class Heatpump():


    def get_heating_data(self, P_e_hp=3000, storage_vol=5000, P_e_storage=300, P_house=200, **kwargs):

        kelvin = lambda x: x + 273.15  # convert Celsius to Kelvin
        celsius = lambda x: x - 273.15  # convert Kelvin to Celsius
        kWh = lambda x: x / 3600000  # convert Joule to kWh
        Joule = lambda x: x * 3600000  # convert kWh to Joule

        # freezing_zone = pd.Interval(left=-6, right=4) # temp between hp starts freezing and has to be thawed
        freezing_zone = 4
        freezing_factor = 1.2  # factor which the efficiency coefficient is divided by in case of being in freezing zone
        vorlauf_max = 50
        aussen_min = -10
        raum = 21  # temperature at which room should be heated
        ruecklauf = 25
        ny = 0.4  # Gütegrad

        temp_start_heat = 15  # temperature under which house has to be heated
        house_energy = 200  # Watt / Kelvin

        timeframe_h = 24  # hours
        timeframe_s = timeframe_h * 3600  # seconds

        # take Newton cooling data

        # 20 and 60 after some DIN-Norms
        t_env = kelvin(20)  # environmental temperature where water storage is located
        t_0 = kelvin(60)  # temperature at which water storage power is measured
        # coef_h = 0.000625

        # specifications of water storage (taken from real example)
        storage_power = 137  # Watts
        h2o_energy = 4200  # Joule / (kilogram * Kelvin)
        storage_volume = 6000  # Liter = kilogram -> to be defined later on

        # specifications of heatpump
        hp_power = 7000  # Watts
        hp_power_opt = 5000  # Watts -> Power in which heatpump works most efficiently (can be taken into calculation later)
        hp_power_el_opt = 3000  # Watts -> Power in which heatpump works most efficiently (can be taken into calculation later)

        # other measures
        groundwater = celsius(t_env)

        vorlauftemp = lambda t: - (vorlauf_max - raum) / (raum - aussen_min) * t + vorlauf_max + aussen_min * (
                    vorlauf_max - raum) / (raum - aussen_min)
        needed_house_energy = lambda t: np.maximum(0,
                                                   raum - t) * house_energy  # returns needed energy to heat house in Watt

        # eta = lambda x: np.maximum(0, np.minimum(300, ny * kelvin(vorlauftemp(x)) / (kelvin(vorlauftemp(x))-kelvin(x))))
        def eta(zu, ab=None):  # theoretischer Wirkungsgrad
            if ab is None:
                ab = vorlauftemp(zu)
            return np.maximum(0, np.minimum(300, ny * kelvin(ab) / (kelvin(ab) - kelvin(zu))))

        def eta_heating(zu, max_ab=None,
                        min_ab=None):  # "gleitender" Wirkungsgrad -> da sich über temperaturerhöhung ändert
            try:
                if type(min_ab) == type(None):
                    min_ab = ruecklauf
                if type(max_ab) == type(None):
                    max_ab = vorlauftemp(zu)
            except:
                print("Error: function 'eta_heating'")
            return ny * (kelvin(zu) * np.log((kelvin(max_ab) - kelvin(zu)) / (kelvin(min_ab) - kelvin(zu))) + kelvin(
                max_ab) - kelvin(min_ab)) / (kelvin(max_ab) - kelvin(min_ab))

        def eta_freezing_adj(zu, max_ab=None, min_ab=None, heating=False):  # if heating = True, then use eta_heating
            try:
                if type(min_ab) == type(None):
                    min_ab = ruecklauf
                if type(max_ab) == type(None):
                    max_ab = vorlauftemp(zu)
            except:
                print("Error: function 'eta_heating'")
            if heating == True:
                e = eta_heating(zu, max_ab, min_ab)
            else:
                e = eta(zu, max_ab)
            # print(type(eta))
            if type(e) == np.float64:
                if zu < 4:
                    e = e / freezing_factor
            else:
                e[zu < 4] = e[zu < 4] / 1.2

            return e

        def newton_cooling(t):
            return t_env + (t_0 - t_env) * np.exp(-coef_h * t)  # = temperature at time t

        def newton_cooling_diff(t):
            return coef_h * (t_env - newton_cooling(t))  # = temperature loss rate at time t

        def newton_cooling_t_0(temp_t, t=7 * 3600):  # enter in Kelvin
            # return ((temp_t - t_env) * np.exp((coef_h * t).T)).T + t_env #returns t_0 for fix temperature after time t
            return (temp_t - t_env) * np.exp((coef_h * t)) + t_env  # returns t_0 for fix temperature after time t

        def newton_cooling_t(temp_t):  # , coef_h = coef_h):
            return -np.log((temp_t - t_env) / (
                        t_0 - t_env)) / coef_h  # returns how much time it takes for t_0 to cool down to wished temperature temp_t

        def storage_power_temp(temp):
            return -h2o_energy * storage_volume * newton_cooling_diff(
                newton_cooling_t(temp))  # returns how much power in Watt is needed to keep temperature at temp level

        def get_d_t():
            return storage_power / (h2o_energy * storage_volume)  # returns heat loss in kelvin per second

        def get_coef_h(d_t):
            return d_t / np.abs(t_env - t_0)  # returns coefficient at which heat decreases in storage

        def get_coef_h_adj(temp, vorlauf, t):
            return - np.log((vorlauf - t_env) / (
                        temp - t_env)) / t  # returns adjusted coef_h for given vorlauf, overheat_temp and time after which vorlauf should be reached

        def storage_power_temp_adj(temp, vorlauf, t):
            return h2o_energy * storage_volume * (coef_h - get_coef_h_adj(temp, vorlauf, t)) * (
                        temp - t_env)  # returns needed power to keep temperature in storage and decrease at given rate

        def get_needed_storage_vol(temp_diff, kWh):
            return (kWh * 3600000) / (
                        temp_diff * h2o_energy)  # returns waterstorage volume to store #kWh with given temp difference

        d_t = get_d_t()
        coef_h = get_coef_h(d_t)

        integral = lambda x, y: trapz(y, x)  # returns the estimated integral for data points x and y
        integral_avg = lambda x, y: integral(x, y) / (
                    x.shape[0] - 1)  # returns the average of the integral for points x and y

        heat = lambda ground_temp, vorlauftemp, eta_t_0: (
                                                                     vorlauftemp - ground_temp) * h2o_energy * storage_volume / eta_t_0  # energy needed to heat water from ground temp to vorlauf with given eta (in Joule)
        overheat = lambda time, vorlauftemp, eta_t_0: (newton_cooling_t_0(kelvin(vorlauftemp),
                                                                          t=time * 3600) - 273.15 - vorlauftemp) * h2o_energy * storage_volume / eta_t_0  # energy needed to heat water storage from vorlauftemp to t_0
        keepheat = lambda time, vorlauftemp, eta_avg: storage_power_temp(kelvin(
            vorlauftemp)) * 3600 * time / eta_avg  # energy needed to keep heat in water storage on vorlauftemp with changing eta

        def over_keep_mix(time, vorlauftemp, eta_t_0,
                          eta_avg):  # get the minimum of a mix of overheating and keeping the heat

            try:
                if time == 0:  # if time equals 0, then there is no optimum that can be computed
                    return vorlauftemp, heat(groundwater, vorlauftemp, eta_t_0), 0

                def target_function(
                        temp):  # returns the energy needed for heating, overheating and keeping the heat level in Joule
                    over = (temp - vorlauftemp) * h2o_energy * storage_volume / eta_t_0
                    keep = storage_power_temp_adj(kelvin(temp), kelvin(vorlauftemp),
                                                  time * 3600) * 3600 * time / eta_avg
                    h = heat(groundwater, vorlauftemp, eta_t_0)

                    return over + keep + h

                lb = vorlauftemp
                ub = max(vorlauftemp, (celsius(newton_cooling_t_0(kelvin(vorlauftemp), time * 3600))))

                cons = Bounds(lb=lb, ub=ub)
                # print(cons)
                min = minimize(fun=target_function, x0=(vorlauftemp + 0.0001,), bounds=cons)

                return min.x[0], min.fun[0], storage_power_temp_adj(kelvin(min.x[0]), kelvin(vorlauftemp),
                                                                    time * 3600)  # TODO insert what is returned
            except:
                return np.nan, np.nan, np.nan

        over_keep_mix_vectorized = np.vectorize(over_keep_mix)

        energy_total_overheat = lambda ground_temp, vorlauftemp, eta_t_0, time: heat(ground_temp, vorlauftemp,
                                                                                     eta_t_0) + overheat(time,
                                                                                                         vorlauftemp,
                                                                                                         eta_t_0)  # just consider overheating
        energy_total_keepheat = lambda ground_temp, vorlauftemp, eta_t_0, eta_avg, time: heat(ground_temp, vorlauftemp,
                                                                                              eta_t_0) + keepheat(time,
                                                                                                                  vorlauftemp,
                                                                                                                  eta_avg)  # just consider keeping heat

        """# Add calculated Data to Dataframe"""

        def add_data(df):
            df['efficiency'] = eta_freezing_adj(df.temperature)
            df['vorlauf_temp'] = vorlauftemp(df.temperature)
            df['energy_needed_water'] = needed_house_energy(df.temperature)
            df['energy_needed_elec'] = needed_house_energy(df.temperature) / df.efficiency
            df['time'] = df.date.dt.time
            df['freezing_zone'] = df['temperature'] < freezing_zone

            return df

        def get_local_extrema(df):  # get local minima and maxima in function
            check_interval = 10

            df['local_min'] = df.iloc[argrelextrema(df.temperature.values, np.less_equal, order=check_interval)[0]][
                'temperature']
            df['local_max'] = df.iloc[argrelextrema(df.temperature.values, np.greater_equal, order=check_interval)[0]][
                'temperature']
            df['local_max_vorlauf'] = \
            df.iloc[argrelextrema(df.vorlauf_temp.values, np.greater_equal, order=check_interval)[0]]['vorlauf_temp']

            # delete doubled high and low points -> take only the later ones
            df.local_min[(pd.notna(df.local_min - df.local_min.shift(-1))) | (
                pd.notna(df.local_min - df.local_min.shift(-2))) | (
                             pd.notna(df.local_min - df.local_min.shift(-3)))] = np.nan
            df.local_max[
                (pd.notna(df.local_max - df.local_max.shift(1))) | (pd.notna(df.local_max - df.local_max.shift(2))) | (
                    pd.notna(df.local_max - df.local_max.shift(3)))] = np.nan

            # high and low points have to alternate
            local_max_index = df[pd.notna(df.local_max)].index  # get indexes of local max
            local_min_index = df[pd.notna(df.local_min)].index  # get indexes of local min

            index = np.sort(
                np.concatenate((local_min_index.values, local_max_index.values)))  # concat and sort both indexes
            index_pd = pd.Series(np.isin(index, local_min_index))

            index = index[(index_pd - index_pd.shift(
                1)) != 0]  # when two max or two min indexes occur, shifting and taking the difference results in 0

            local_min_index = local_min_index[np.isin(local_min_index, index)]  # take only remaining indexes
            local_max_index = local_max_index[np.isin(local_max_index, index)]
            for i, index in enumerate(local_min_index):
                df.loc[df.index > index, 'period_min'] = i

            for i, index in enumerate(local_max_index):
                df.loc[df.index >= index, 'period_max'] = i

            return df

        """# Import data from dwd

        (Deutscher Wetterdienst)
        """

        data_temp = pd.read_table(r'C:\Users\rafaelweinert\PycharmProjects\Heatpump\data\stundenwerte_TU_VS_akt\produkt_tu_stunde_20201216_20220618_05229.txt',
                                  sep=';')  # Villingen-Schwenningen
        # data_temp = pd.read_table('data/stundenwerte_TU_05100_19410101_20211231_hist/produkt_tu_stunde_19410101_20211231_05100.txt', sep=';') #Trier?
        data_temp = data_temp.drop(['eor', 'QN_9', 'RF_TU'], axis=1)  # drop end-of-row indicator
        data_temp.rename(columns={'TT_TU': 'temperature', 'MESS_DATUM': 'date', }, inplace=True)
        data_temp.date = pd.to_datetime(data_temp.date.astype(str), format='%Y%m%d%H')
        data_temp[data_temp == -999] = np.nan  # rename missing values to numpy nan
        data_temp = data_temp[
            (data_temp.date >= pd.to_datetime('2021-10-02')) & (data_temp.date <= pd.to_datetime('2022-03-31'))]
        # data_temp = data_temp[(data_temp.date.dt.month >= 10) | (data_temp.date.dt.month <= 3)]
        data_temp.reset_index(drop=True, inplace=True)

        data_temp = add_data(data_temp)
        data_temp = get_local_extrema(data_temp)

        px.line(data_temp, x='date', y='temperature', title='temperature')

        """# Visualizations for temperature"""

        px.histogram(data_temp.freezing_zone, title='Number of hours in Freezing zone (<4°C)')

        fig = make_subplots(specs=[[{"secondary_y": True}]])

        fig.add_trace(go.Scatter(x=data_temp.date, y=data_temp.temperature, name='temperature'), secondary_y=False)
        fig.add_trace(go.Scatter(x=data_temp.date, y=data_temp.local_min, mode='markers', name='local_minima'),
                      secondary_y=False)
        fig.add_trace(
            go.Scatter(x=data_temp.date, y=data_temp.local_max, mode='markers', hovertext=data_temp.period_max,
                       name='local_maxima'), secondary_y=False)
        fig.update_layout(dict1={'title': 'temperature with high and low points'})
        fig.show()

        data_temp

        """# Weather forecast

        -> the majority of the energy need for water is needed in profitable period (compare two lines)
        """

        def get_weather_forecast(address, days=4, key='676ea42d0beb4636b96160506221607'):

            # Get coordinates for adress
            geolocator = Nominatim(user_agent="user")
            location = geolocator.geocode(address)

            # Get nearby weather stations
            stations = meteostat.Stations()
            stations = stations.nearby(location.latitude, location.longitude)
            station = stations.fetch(1)

            # Print DataFrame
            print(station)
            weather_data = requests.get(
                f'http://api.weatherapi.com/v1/forecast.json?key={key}&q={address}&days={days}&aqi=no&alerts=no')

            f = pd.DataFrame.from_dict(json.loads(weather_data.content.decode())['forecast']['forecastday'])
            forecast = pd.json_normalize(pd.json_normalize(f['hour']).stack())

            forecast.time = pd.to_datetime(forecast.time, )
            forecast = forecast[['time', 'temp_c']]

            forecast.columns = ['date', 'temperature']
            return forecast

        forecast = get_weather_forecast('Karlsruhe', days=4)
        forecast.temperature = forecast.temperature - 17

        forecast = add_data(forecast)
        forecast = get_local_extrema(forecast)
        forecast = forecast[forecast.date > pd.to_datetime('now')]

        fig = make_subplots(specs=[[{"secondary_y": True}]])

        fig.add_trace(go.Scatter(x=forecast.date, y=forecast.temperature), secondary_y=False)
        fig.add_trace(go.Scatter(x=forecast.date, y=forecast.local_min, mode='markers', name='local_minima'),
                      secondary_y=False)
        fig.add_trace(go.Scatter(x=forecast.date, y=forecast.local_max, mode='markers', hovertext=forecast.period_max,
                                 name='local_maxima'), secondary_y=False)
        # fig.add_trace(go.Scatter(x = forecast.date, y = forecast.energy_needed_water, mode='lines', name='energy_needed_water'), secondary_y= True)
        # fig.add_trace(go.Scatter(x = forecast.date, y = forecast.energy_needed_water_integ, mode='lines', name='energy_needed_water_integ'), secondary_y= True)
        # fig.add_trace(go.Scatter(x = forecast.date, y = forecast.en_produced_total, mode='lines', name='en_produced_total'), secondary_y= True)
        # fig.add_trace(go.Scatter(x = forecast.date, y = forecast.en_produced_hourly, mode='lines', name='en_produced_hourly'), secondary_y= True)
        # fig.add_trace(go.Scatter(x = forecast.date, y = forecast.eta_hp_heating_nextmin, mode='lines', name='eta_nextmin'), secondary_y= False)
        # fig.add_trace(go.Scatter(x = forecast.date, y = forecast.eta_hp_heating_nextmin_freezing_adjusted, mode='lines', name='eta_nextmin_freezing_adj.'), secondary_y= False)
        fig.add_hline(y=4)
        fig.update_layout(title='temperature Simulation (From forecast)')
        fig.show()

        """# Optimization problem 1

        simplified version: period_lim = None
        assumptions:
        - Unlimited Storage
        - Heating temperature is equal for every temperature (max_ab_i = x for all i)
        - Efficiency only depends on environmental Temperature
        - ruecklauf is equal for each temperature outside -> whenever storage temperature is higher, it is sufficient for heating and can be adjusted by water flow



        extended version: period_lim = periods restrict optimization on
        """

        def heating_opt(df, period_lim=None, period=None, vorl=None):
            # df_per = df[(df.period_max == p) | (df.period_min == p)]
            if type(period) != type(None):
                df_per = df[(df.period_max == period) | (df.period_min == period)]
            else:
                df_per = df.copy()
            per_min = df[df.period_min == period]  # period from miminum to minimum
            per_max = df[df.period_max == period]  # period from max to max

            per = df_per.shape[0]  # number of hours in period
            M_timedelta = np.empty((per, per))

            for i in range(per):  # build timedelta matrix
                M_timedelta[i, range(0, i + 1)] = np.flip(np.linspace(0, i, i + 1))
            M_timedelta = M_timedelta.T
            vorlauf = df_per.vorlauf_temp.values

            M_vorlauf = celsius(newton_cooling_t_0(kelvin(vorlauf), M_timedelta * 3600))
            M_vorlauf = np.triu(
                M_vorlauf)  # needed heating temperature when heating in period i and using the heat in period j (i X j - Matrix)
            M_temperature = np.tile(df_per.temperature.values, (len(df_per.temperature.values), 1)).T

            M_eta = eta_freezing_adj(zu=M_temperature, max_ab=M_vorlauf)  # , min_ab=M_vorlauf)

            if type(period_lim) != type(None):

                M_vorlauf_max = celsius(newton_cooling_t_0(kelvin(vorl), t=M_timedelta * 3600))
                M_eta_heating = eta_freezing_adj(zu=M_temperature, max_ab=M_vorlauf_max, heating=True)
                M_eta_heating = np.triu(M_eta_heating, k=1)
                M_eta = np.tril(M_eta)

                M_eta = M_eta_heating + M_eta
                M_eta = np.triu(M_eta)  # upper diagonal
                M_eta = np.tril(M_eta, k=period_lim)

            else:
                M_eta = np.triu(M_eta)  # upper diagonal

            P_w_sum = df_per.energy_needed_water.values  # needed energy in water for each hour
            P_e_max = np.full((per,), hp_power_el_opt)  # maximal available energy (electricity) per hour -> e.g. 2500W
            P_e_available = P_e_max.copy()  # at each period still available 'capacity' in electricity for future hours

            P_e = np.zeros((per, per))  # Matrix with used energy (e) from each period for each period that is optimized
            P_e_a = np.empty((per, per))  # Matrix to later calculate P_e

            for j in range(per):
                P_w_j = P_w_sum[j]  # energy needed (w) in hour j

                M_eta_j = M_eta[:, j]  # efficiencies in each hour
                eta_sorted_j_index = np.flip(np.argsort(M_eta_j))  # get indexes of highes values in M_eta_j
                eta_sorted_j = M_eta_j[eta_sorted_j_index]  # get highest values of M_eta_j (descending order)

                P_w_available_j = P_e_available[
                                      eta_sorted_j_index] * eta_sorted_j  # available energy (w) for hour j in descending cost order

                P_w_available_j_cumsum = P_w_available_j.cumsum()  # get cumsum for available P_w (descending cost order)
                P_w_available_j_cumsum_shift = P_w_available_j_cumsum - P_w_available_j  # shifted cumsum

                # calculate which P_w_available_j values are used fully and which only partly (only one will be used partly)
                cumsum_bin = P_w_available_j_cumsum < P_w_j  # binary values if cumsum is smaller than needed energy
                cumsum_shift_bin = P_w_available_j_cumsum_shift < P_w_j  # binary values if cumsum_shift is smaller than needed energy

                idx = (cumsum_bin == False) & (cumsum_shift_bin == True)  # only True for partly used P_w_j
                percent = (P_w_available_j_cumsum[idx] - P_w_j) / P_w_available_j[
                    idx]  # get remaining share of partly used P_w_j
                P_e_percent_remaining = np.ones((per,))  # first, all remaining values are set to 1 (=100%)
                P_e_percent_remaining[cumsum_bin] = 0  # set all fully used P_w_j to 0
                P_e_percent_remaining[idx] = percent  # set partly used to calculated percentage

                P_e_percent_remaining_sorted_orig = P_e_percent_remaining[
                    np.argsort(eta_sorted_j_index)]  # set remaining shares to original order
                P_e_available = P_e_available * P_e_percent_remaining_sorted_orig  # get new P_e_available
                P_e_a[:, j] = P_e_available

            P_e_a = np.triu(P_e_a)  # take only upper diagonal
            P_e[:, 0] = P_e_a[:, 0]
            for j in range(1, per):
                P_e[:, j] = np.abs(P_e_a[:, j] - P_e_a[:, j - 1])
            diagonal = np.zeros((per, per))
            np.fill_diagonal(diagonal, hp_power_el_opt)

            P_e = np.abs(diagonal - P_e)
            P_w = (P_e * M_eta)

            df_per['P_e_opt'] = P_e.sum(1)

            return M_eta, P_e, P_w, df_per

        M_eta, P_e, P_w, df_per = heating_opt(data_temp)

        # maximally possible saving potential, if there was no storage limit
        df_per.P_e_opt.sum() / df_per.energy_needed_elec.sum()

        """## Storage volume"""

        # px.line(temperatures.max(axis=0))

        """## Visualizations for Optimized Heating

        ## Energy storage forecast
        """

        fig = make_subplots(specs=[[{"secondary_y": True}]])

        fig.add_trace(go.Scatter(x=df_per.date, y=df_per.temperature, name='temperature', mode='lines'),
                      secondary_y=False)
        fig.add_trace(go.Scatter(x=df_per.date, y=df_per.P_e_opt, name='Power (e) optimized', mode='markers'),
                      secondary_y=True)
        fig.add_trace(
            go.Scatter(x=df_per.date, y=df_per.energy_needed_elec, name='Power (e) not optimized', mode='markers'),
            secondary_y=True)
        fig.update_layout(title='Power (elec) for optimized and not optimized heating')
        fig.write_html('optimization1.html')

        # check if optimized power line is always above needed power line
        fig = make_subplots(specs=[[{"secondary_y": True}]])

        fig.add_trace(
            go.Scatter(x=df_per.date, y=df_per.energy_needed_water.cumsum(), name='power needed (water)', mode='lines'),
            secondary_y=False)
        fig.add_trace(go.Scatter(x=df_per.date, y=P_w.sum(1).cumsum(), name='optimized power (water)', mode='lines'),
                      secondary_y=False)
        fig.update_layout(title='Cumulated Energy (water) stored vs energy needed in Watthours')

        fig = make_subplots(specs=[[{"secondary_y": True}]])

        fig.add_trace(go.Scatter(x=df_per.date, y=df_per.P_e_opt.cumsum(), name='Power (e) optimized', mode='lines'),
                      secondary_y=False)
        fig.add_trace(go.Scatter(x=df_per.date, y=df_per.energy_needed_elec.cumsum(), name='Power (e) not optimized',
                                 mode='lines'), secondary_y=False)
        fig.update_layout(title='Cumulated Power (elec) for optimized and not optimized heating')

        # check visually if needed power equals optimized power supply for every period -> lines have to match exactly
        fig = make_subplots(specs=[[{"secondary_y": True}]])

        fig.add_trace(
            go.Scatter(x=df_per.date, y=df_per.energy_needed_water.cumsum(), name='power needed (water)', mode='lines'),
            secondary_y=False)
        fig.add_trace(go.Scatter(x=df_per.date, y=P_w.sum(0).cumsum(), name='optimized power (water)', mode='lines'),
                      secondary_y=False)

        """## Heat maps"""

        fig = px.imshow(P_e)
        fig.write_html('P_e' + '.html', auto_open=False)

        fig = px.imshow(M_eta)
        # fig.show()
        fig.write_html('M_eta' + '.html', auto_open=False)

        fig = px.imshow(P_w)
        fig.update_layout(height=1000)
        fig.write_html('P_w' + '.html', auto_open=False)

        px.line(data_temp.energy_needed_water.rolling(24).sum()).write_html('energy_needed_water' + '.html',
                                                                            auto_open=False)

        """# Optimization Problem 2

        Bei Betrachtung des Speichervolumens fällt mit einem Temperaturdelta im Speicher das gespeicherte Volumen selten über 6000 Liter -> künftige Annahme: 6000 Liter Tank und Vorlauftemperatur im Speicher konstantes Temperaturdelta von 10 Grad
        Speicher hat im Nullzustand 27 Grad, wird auf max_vorlauf Grad erhitzt. Vorlauftemp sonst ist immer 30 Grad
        """

        # optimize over vorlauftemp

        def heating_optimization(df, period_lim=None, period=None, **kwargs):
            if type(period) != type(None):
                df_per = df[(df.period_max == period) | (df.period_min == period)]
            else:
                df_per = df.copy()

            per = df_per.shape[0]  # number of hours in period
            M_timedelta = np.empty((per, per))

            for i in range(per):  # build timedelta matrix
                M_timedelta[i, range(0, i + 1)] = np.flip(np.linspace(0, i, i + 1))
            M_timedelta = M_timedelta.T
            vorlauf = df_per.vorlauf_temp.values

            M_vorlauf = celsius(newton_cooling_t_0(kelvin(vorlauf), M_timedelta * 3600))
            M_vorlauf = np.triu(
                M_vorlauf)  # needed heating temperature when heating in period i and using the heat in period j (i X j - Matrix)
            M_temperature = np.tile(df_per.temperature.values, (len(df_per.temperature.values), 1)).T

            def target_fun(vorl, optimizing=True, constraint=False):

                M_eta = eta_freezing_adj(zu=M_temperature, max_ab=M_vorlauf)  # , min_ab=M_vorlauf)

                if type(period_lim) != type(None):
                    M_vorlauf_max = celsius(newton_cooling_t_0(kelvin(vorl), t=M_timedelta * 3600))
                    M_eta_heating = eta_freezing_adj(zu=M_temperature, max_ab=M_vorlauf_max, heating=True)
                    M_eta_heating = np.triu(M_eta_heating, k=1)
                    M_eta = np.tril(M_eta)

                    M_eta = M_eta_heating + M_eta
                    M_eta = np.triu(M_eta)  # upper diagonal
                    M_eta = np.tril(M_eta, k=period_lim)

                else:
                    M_eta = np.triu(M_eta)  # upper diagonal

                P_w_sum = df_per.energy_needed_water.values  # needed energy in water for each hour
                P_e_max = np.full((per,),
                                  hp_power_el_opt)  # maximal available energy (electricity) per hour -> e.g. 2500W
                P_e_available = P_e_max.copy()  # at each period still available 'capacity' in electricity for future hours

                P_e = np.zeros(
                    (per, per))  # Matrix with used energy (e) from each period for each period that is optimized
                P_e_a = np.empty((per, per))  # Matrix to later calculate P_e

                for j in range(per):
                    P_w_j = P_w_sum[j]  # energy needed (w) in hour j

                    M_eta_j = M_eta[:, j]  # efficiencies in each hour
                    eta_sorted_j_index = np.flip(np.argsort(M_eta_j))  # get indexes of highes values in M_eta_j
                    eta_sorted_j = M_eta_j[eta_sorted_j_index]  # get highest values of M_eta_j (descending order)

                    P_w_available_j = P_e_available[
                                          eta_sorted_j_index] * eta_sorted_j  # available energy (w) for hour j in descending cost order

                    P_w_available_j_cumsum = P_w_available_j.cumsum()  # get cumsum for available P_w (descending cost order)
                    P_w_available_j_cumsum_shift = P_w_available_j_cumsum - P_w_available_j  # shifted cumsum

                    # calculate which P_w_available_j values are used fully and which only partly (only one will be used partly)
                    cumsum_bin = P_w_available_j_cumsum < P_w_j  # binary values if cumsum is smaller than needed energy
                    cumsum_shift_bin = P_w_available_j_cumsum_shift < P_w_j  # binary values if cumsum_shift is smaller than needed energy

                    idx = (cumsum_bin == False) & (cumsum_shift_bin == True)  # only True for partly used P_w_j
                    percent = (P_w_available_j_cumsum[idx] - P_w_j) / P_w_available_j[
                        idx]  # get remaining share of partly used P_w_j
                    P_e_percent_remaining = np.ones((per,))  # first, all remaining values are set to 1 (=100%)
                    P_e_percent_remaining[cumsum_bin] = 0  # set all fully used P_w_j to 0
                    P_e_percent_remaining[idx] = percent  # set partly used to calculated percentage

                    P_e_percent_remaining_sorted_orig = P_e_percent_remaining[
                        np.argsort(eta_sorted_j_index)]  # set remaining shares to original order
                    P_e_available = P_e_available * P_e_percent_remaining_sorted_orig  # get new P_e_available
                    P_e_a[:, j] = P_e_available

                P_e_a = np.triu(P_e_a)  # take only upper diagonal
                P_e[:, 0] = P_e_a[:, 0]
                for j in range(1, per):
                    P_e[:, j] = np.abs(P_e_a[:, j] - P_e_a[:, j - 1])
                diagonal = np.zeros((per, per))
                np.fill_diagonal(diagonal, hp_power_el_opt)

                P_e = np.abs(diagonal - P_e)
                P_w = (P_e * M_eta)

                if optimizing == False:
                    df_per['P_e_opt'] = P_e.sum(1)
                    return M_eta, P_e, P_w, df_per

                if constraint == True:
                    return P_w

                return P_e.sum()

            def constraint_fun(vorl):

                P_w = target_fun(vorl, constraint=True)
                stor_vol = np.zeros((P_w.shape[0],))

                for i in range(P_w.shape[0] - 1):
                    stor_vol[i] = get_needed_storage_vol(vorl - ruecklauf, P_w[:(i + 1), (i + 1):].sum() / 1000)
                max_stor = max(stor_vol)
                return storage_volume - max_stor

            if 'optimize' in kwargs:  # same as heating_opt
                if kwargs['optimize'] == False:
                    if 'vorlauf' in kwargs:
                        M_eta_opt, P_e_opt, P_w_opt, df_per_opt = target_fun(kwargs['vorlauf'], optimizing=False)
                        return M_eta_opt, P_e_opt, P_w_opt, df_per_opt
                    else:
                        return 'Must enter a temperature'

            # optimization problem
            lb = max(vorlauf) + 0.001
            ub = 100
            bounds = Bounds(lb=lb, ub=ub)

            min = minimize(target_fun, x0=(lb + 1,), bounds=bounds, constraints={'type': 'ineq',
                                                                                 'fun': constraint_fun})

            M_eta_opt, P_e_opt, P_w_opt, df_per_opt = target_fun(min.x[0],
                                                                 optimizing=False)  # get values for optimal heating temperature

            return min, M_eta_opt, P_e_opt, P_w_opt, df_per_opt

        minimum, M_eta_opt, P_e_opt, P_w_opt, df_per_opt = heating_optimization(data_temp, period=2, period_lim=21)

        """# Optimization with varying preheatings

        TODO
        add P_w to big P_w with all periods
        """

        max_periods = 36
        index = 0
        data_temp_new = pd.DataFrame()

        for p in data_temp.period_max.dropna().unique():
            print('period: ', p)

            min_x = []
            min_fun = []

            for i in range(max_periods):
                # print('period_lim: ', i)
                minimum, M_eta, P_e, P_w, df_per = heating_optimization(data_temp.loc[index:, :], period=p,
                                                                        period_lim=i)

                min_x.append(minimum.x[0])
                min_fun.append(minimum.fun)

            temp = min_x[np.argmin(min_fun)]
            # TODO change period_lim = i to actual period_lim
            # M_eta, P_e,  P_w, df_per = heating_opt(data_temp.loc[index:, :], period=p, period_lim=np.argmin(min_fun), vorl=temp)
            M_eta, P_e, P_w, df_per = heating_optimization(data_temp.loc[index:, :], period=p,
                                                           period_lim=np.argmin(min_fun), optimize=False, vorlauf=temp)

            stor_vol = np.zeros((P_w.shape[0],))

            for i in range(P_w.shape[0] - 1):
                stor_vol[i] = get_needed_storage_vol(temp - ruecklauf, P_w[:(i + 1), (i + 1):].sum() / 1000)
            df_per['Storage_vol'] = stor_vol
            df_per['Preheat_temp'] = temp
            df_per['Period_lim'] = np.argmin(min_fun)
            index = df_per[df_per.Storage_vol == 0].Storage_vol.index[-1] + 1

            data_temp_new = pd.concat([data_temp_new, df_per])

        fig = make_subplots(specs=[[{"secondary_y": True}]])
        fig.add_trace(go.Scatter(x=data_temp_new.date, y=data_temp_new.temperature, name='temperature', mode='lines'),
                      secondary_y=False)
        fig.add_trace(
            go.Scatter(x=data_temp_new.date, y=data_temp_new.P_e_opt, name='Power (e) optimized', mode='markers'),
            secondary_y=True)
        fig.add_trace(
            go.Scatter(x=data_temp_new.date, y=data_temp_new.energy_needed_elec, name='Power (e) not optimized',
                       mode='markers'),
            secondary_y=True)
        fig.add_trace(go.Scatter(x=data_temp_new.date, y=data_temp_new.Storage_vol, name='Storage Volume'),
                      secondary_y=True)
        fig.add_trace(go.Scatter(x=data_temp_new.date, y=data_temp_new.Preheat_temp, name='Preheat Temperature'),
                      secondary_y=False)
        fig.add_trace(go.Scatter(x=data_temp_new.date, y=data_temp_new.vorlauf_temp, name='Vorlauftemperatur'),
                      secondary_y=False)
        fig.update_layout(title='Power (elec) for optimized and not optimized heating')
        fig.write_html('optimization_data_temp.html')

        data_temp_new.P_e_opt.sum() / data_temp_new.energy_needed_elec.sum()

        stored_water_energy = np.zeros((P_w.shape[0],))
        stor_vol_fc = np.zeros((P_w.shape[0],))
        for i in range(P_w.shape[0] - 1):
            stored_water_energy[i] = P_w[:(i + 1), (i + 1):].sum()
            # print(temperatures.min(axis=0)[i] - ruecklauf)
            stor_vol_fc[i] = get_needed_storage_vol(minimum.x[0] - ruecklauf, P_w[:(i + 1), (i + 1):].sum() / 1000)

        fig = go.Figure()
        fig.add_trace(go.Scatter(y=stor_vol_fc))
        fig.add_trace(go.Scatter(x=data_temp.temperature.index, y=data_temp.temperature))
        fig.add_trace(
            go.Scatter(x=data_temp.local_min.dropna().index, y=np.full(data_temp.local_min.shape, 5), mode='markers'))
        fig.add_trace(go.Scatter(x=data_temp.index, y=df_per.P_e_opt, name='Power (e) optimized', mode='markers'))
        fig.add_trace(
            go.Scatter(x=data_temp.local_max.dropna().index, y=np.full(data_temp.local_max.shape, 5), mode='markers'))
        fig.write_html('water_storage.html')

        M_eta_fc, P_e_fc, P_w_fc, df_per_fc = heating_opt(forecast, period_lim=27, vorl=minimum.x[0])

        stored_water_energy = np.zeros((P_w_fc.shape[0],))
        stor_vol_fc = np.zeros((P_w_fc.shape[0],))
        for i in range(P_w_fc.shape[0] - 1):
            stored_water_energy[i] = P_w_fc[:(i + 1), (i + 1):].sum()
            # print(temperatures.min(axis=0)[i] - ruecklauf)
            stor_vol_fc[i] = get_needed_storage_vol(33.3428 - ruecklauf, P_w_fc[:(i + 1), (i + 1):].sum() / 1000)

        fig = make_subplots(specs=[[{"secondary_y": True}]])

        fig.add_trace(go.Scatter(x=df_per_fc.date, y=df_per_fc.temperature, name='temperature', mode='lines'),
                      secondary_y=False)
        fig.add_trace(go.Scatter(x=df_per_fc.date, y=df_per_fc.P_e_opt, name='Power (e) optimized', mode='markers'),
                      secondary_y=True)
        fig.add_trace(go.Scatter(x=df_per_fc.date, y=df_per_fc.energy_needed_elec, name='Power (e) not optimized',
                                 mode='markers'),
                      secondary_y=True)
        fig.update_layout(title='Power (elec) for optimized and not optimized heating')
        fig.write_html('optimization_fc.html')

        forecast.energy_needed_elec.sum()

        stored_water_energy = np.zeros((P_w_opt.shape[0],))
        stor_vol = np.zeros((P_w_opt.shape[0],))
        for i in range(P_w_opt.shape[0] - 1):
            stored_water_energy[i] = P_w_opt[:(i + 1), (i + 1):].sum()
            stor_vol[i] = get_needed_storage_vol(minimum.x[0] - ruecklauf, P_w_opt[:(i + 1), (i + 1):].sum() / 1000)
        max(get_needed_storage_vol(minimum.x[0] - ruecklauf, stored_water_energy / 1000))

        px.line(stor_vol).write_html('stor_vol3.html')



        print(data_temp_new.P_e_opt.sum() / data_temp_new.energy_needed_elec.sum())
        return data_temp_new.P_e_opt.sum() / data_temp_new.energy_needed_elec.sum()


hp = Heatpump()
hp.get_heating_data()

'''
        """# Optimization with varying period_lim"""

        # optimize over vorlauftemp

        def heating_optimization_period(df, period = None):

            def tf(period_lim=None):
                period_lim = np.array(period_lim, dtype=int)
                print(period_lim)
                if type(period) != type(None):
                    df_per = df[(df.period_max == period) | (df.period_min == period)]
                else:
                    df_per = df.copy()

                per = df_per.shape[0]  #number of hours in period
                M_timedelta = np.empty((per, per))

                for i in range(per):  # build timedelta matrix
                    M_timedelta[i, range(0, i + 1)] = np.flip(np.linspace(0, i, i + 1))
                M_timedelta = M_timedelta.T
                vorlauf = df_per.vorlauf_temp.values

                M_vorlauf = celsius(newton_cooling_t_0(kelvin(vorlauf), M_timedelta * 3600))
                M_vorlauf = np.triu(M_vorlauf)  # needed heating temperature when heating in period i and using the heat in period j (i X j - Matrix)
                M_temperature = np.tile(df_per.temperature.values, (len(df_per.temperature.values), 1)).T

                def target_fun(vorl, optimizing=True):

                    M_eta = eta_freezing_adj(zu=M_temperature, max_ab=M_vorlauf)  #, min_ab=M_vorlauf)

                    if type(period_lim) != type(None):
                        M_vorlauf_max = celsius(newton_cooling_t_0(kelvin(vorl), t=M_timedelta * 3600))
                        M_eta_heating = eta_freezing_adj(zu=M_temperature, max_ab=M_vorlauf_max, heating=False)
                        M_eta_heating = np.triu(M_eta_heating, k=1)
                        M_eta = np.tril(M_eta)

                        M_eta = M_eta_heating + M_eta
                        M_eta = np.triu(M_eta)  # upper diagonal
                        M_eta_final = np.zeros(M_eta.shape)
                        for i in range(M_eta.shape[0]):

                            M_eta_final[i, i:i +period_lim[i]+1] = M_eta[i, i: i + period_lim[i] +1]

                        M_eta = M_eta_final
                        #M_eta = np.tril(M_eta, k=period_lim)

                    else:
                        M_eta = np.triu(M_eta)  # upper diagonal

                    P_w_sum = df_per.energy_needed_water.values  # needed energy in water for each hour
                    P_e_max = np.full((per,), hp_power_el_opt)  # maximal available energy (electricity) per hour -> e.g. 2500W
                    P_e_available = P_e_max.copy()  # at each period still available 'capacity' in electricity for future hours

                    P_e = np.zeros((per, per))  # Matrix with used energy (e) from each period for each period that is optimized
                    P_e_a = np.empty((per, per))  # Matrix to later calculate P_e

                    for j in range(per):
                        P_w_j = P_w_sum[j]  # energy needed (w) in hour j

                        M_eta_j = M_eta[:, j]  # efficiencies in each hour
                        eta_sorted_j_index = np.flip(np.argsort(M_eta_j))  # get indexes of highes values in M_eta_j
                        eta_sorted_j = M_eta_j[eta_sorted_j_index]  # get highest values of M_eta_j (descending order)

                        P_w_available_j = P_e_available[
                                              eta_sorted_j_index] * eta_sorted_j  # available energy (w) for hour j in descending cost order

                        P_w_available_j_cumsum = P_w_available_j.cumsum()  # get cumsum for available P_w (descending cost order)
                        P_w_available_j_cumsum_shift = P_w_available_j_cumsum - P_w_available_j  # shifted cumsum

                        # calculate which P_w_available_j values are used fully and which only partly (only one will be used partly)
                        cumsum_bin = P_w_available_j_cumsum < P_w_j  # binary values if cumsum is smaller than needed energy
                        cumsum_shift_bin = P_w_available_j_cumsum_shift < P_w_j  # binary values if cumsum_shift is smaller than needed energy

                        idx = (cumsum_bin == False) & (cumsum_shift_bin == True)  # only True for partly used P_w_j
                        percent = (P_w_available_j_cumsum[idx] - P_w_j) / P_w_available_j[
                            idx]  # get remaining share of partly used P_w_j
                        P_e_percent_remaining = np.ones((per,))  # first, all remaining values are set to 1 (=100%)
                        P_e_percent_remaining[cumsum_bin] = 0  # set all fully used P_w_j to 0
                        P_e_percent_remaining[idx] = percent  # set partly used to calculated percentage

                        P_e_percent_remaining_sorted_orig = P_e_percent_remaining[
                            np.argsort(eta_sorted_j_index)]  # set remaining shares to original order
                        P_e_available = P_e_available * P_e_percent_remaining_sorted_orig  # get new P_e_available
                        P_e_a[:, j] = P_e_available

                    P_e_a = np.triu(P_e_a)  # take only upper diagonal
                    P_e[:, 0] = P_e_a[:, 0]
                    for j in range(1, per):
                        P_e[:, j] = np.abs(P_e_a[:, j] - P_e_a[:, j - 1])
                    diagonal = np.zeros((per, per))
                    np.fill_diagonal(diagonal, hp_power_el_opt)

                    P_e = np.abs(diagonal - P_e)

                    if optimizing == False:
                        P_w = (P_e * M_eta)
                        df_per['P_e_opt'] = P_e.sum(1)
                        return M_eta, P_e, P_w, df_per

                    return P_e.sum()

                def constraint_fun(vorl):

                    M_eta = eta_freezing_adj(zu=M_temperature, max_ab=M_vorlauf)  #, min_ab=M_vorlauf)

                    if type(period_lim) != type(None):
                        M_vorlauf_max = celsius(newton_cooling_t_0(kelvin(vorl), t=M_timedelta * 3600))
                        M_eta_heating = eta_freezing_adj(zu=M_temperature, max_ab=M_vorlauf_max, heating=False)
                        M_eta_heating = np.triu(M_eta_heating, k=1)
                        M_eta = np.tril(M_eta)

                        M_eta = M_eta_heating + M_eta
                        M_eta = np.triu(M_eta)  # upper diagonal

                        M_eta_final = np.zeros(M_eta.shape)
                        for i in range(M_eta.shape[0]):
                            M_eta_final[i, i:i +period_lim[i]+1] = M_eta[i, i: i + period_lim[i] +1]

                        M_eta = M_eta_final#np.tril(M_eta, k=period_lim)

                    else:
                        M_eta = np.triu(M_eta)  # upper diagonal

                    P_w_sum = df_per.energy_needed_water.values  # needed energy in water for each hour
                    P_e_max = np.full((per,), hp_power_el_opt)  # maximal available energy (electricity) per hour -> e.g. 2500W
                    P_e_available = P_e_max.copy()  # at each period still available 'capacity' in electricity for future hours

                    P_e = np.zeros((per, per))  # Matrix with used energy (e) from each period for each period that is optimized
                    P_e_a = np.empty((per, per))  # Matrix to later calculate P_e

                    for j in range(per):
                        P_w_j = P_w_sum[j]  # energy needed (w) in hour j

                        M_eta_j = M_eta[:, j]  # efficiencies in each hour
                        eta_sorted_j_index = np.flip(np.argsort(M_eta_j))  # get indexes of highes values in M_eta_j
                        eta_sorted_j = M_eta_j[eta_sorted_j_index]  # get highest values of M_eta_j (descending order)

                        P_w_available_j = P_e_available[
                                              eta_sorted_j_index] * eta_sorted_j  # available energy (w) for hour j in descending cost order

                        P_w_available_j_cumsum = P_w_available_j.cumsum()  # get cumsum for available P_w (descending cost order)
                        P_w_available_j_cumsum_shift = P_w_available_j_cumsum - P_w_available_j  # shifted cumsum

                        # calculate which P_w_available_j values are used fully and which only partly (only one will be used partly)
                        cumsum_bin = P_w_available_j_cumsum < P_w_j  # binary values if cumsum is smaller than needed energy
                        cumsum_shift_bin = P_w_available_j_cumsum_shift < P_w_j  # binary values if cumsum_shift is smaller than needed energy

                        idx = (cumsum_bin == False) & (cumsum_shift_bin == True)  # only True for partly used P_w_j
                        percent = (P_w_available_j_cumsum[idx] - P_w_j) / P_w_available_j[
                            idx]  # get remaining share of partly used P_w_j
                        P_e_percent_remaining = np.ones((per,))  # first, all remaining values are set to 1 (=100%)
                        P_e_percent_remaining[cumsum_bin] = 0  # set all fully used P_w_j to 0
                        P_e_percent_remaining[idx] = percent  # set partly used to calculated percentage

                        P_e_percent_remaining_sorted_orig = P_e_percent_remaining[
                            np.argsort(eta_sorted_j_index)]  # set remaining shares to original order
                        P_e_available = P_e_available * P_e_percent_remaining_sorted_orig  # get new P_e_available
                        P_e_a[:, j] = P_e_available

                    P_e_a = np.triu(P_e_a)  # take only upper diagonal
                    P_e[:, 0] = P_e_a[:, 0]
                    for j in range(1, per):
                        P_e[:, j] = np.abs(P_e_a[:, j] - P_e_a[:, j - 1])
                    diagonal = np.zeros((per, per))
                    np.fill_diagonal(diagonal, hp_power_el_opt)

                    P_e = np.abs(diagonal - P_e)
                    P_w = (P_e * M_eta)

                    stor_vol = np.zeros((P_w.shape[0],))

                    for i in range(P_w.shape[0] - 1):
                        stor_vol[i] = get_needed_storage_vol(vorl - ruecklauf, P_w[:(i + 1), (i + 1):].sum() / 1000)
                    max_stor = max(stor_vol)
                    return storage_volume - max_stor

                lb = ruecklauf + 0.01
                ub = 100
                bounds = Bounds(lb=lb, ub=ub)

                min = minimize(target_fun, x0=(lb + 1,), bounds=bounds, constraints={'type': 'ineq',
                                                                                     'fun': constraint_fun})

                M_eta_opt, P_e_opt, P_w_opt, df_per_opt = target_fun(min.x[0], optimizing=False)

                return min.fun#, M_eta_opt, P_e_opt, P_w_opt, df_per_opt

            lower_b = np.zeros((df.shape[0], ), dtype=int)
            upper_b = np.full((df.shape[0], ), 36, dtype=int)

            b = Bounds(lower_b, upper_b)
            mini = minimize(tf, x0=lower_b + 2, bounds=b, constraints = {'type':'eq','fun': lambda x : max([x[i]-int(x[i]) for i in range(len(x))])})

            return mini

        heating_optimization_period(forecast)

        '''
