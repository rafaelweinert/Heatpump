from Funs import Functions
import Cons as c

import numpy as np
import pandas as pd
from scipy.signal import argrelextrema


class Preprocessor:

    def __init__(self, flow_max, room_temp, out_min, return_flow, house_energy, storage_volume, storage_power, ny):
        self.f = Functions(flow_max, room_temp, out_min, return_flow, house_energy, storage_volume, storage_power, ny)

    def add_data(self, df):
        df['efficiency'] = self.f.eta_freezing_adj(df.temperature, heating=c.h)
        df['flow_temp'] = self.f.flowtemp(df.temperature)
        df['energy_needed_water'] = self.f.needed_house_energy(df.temperature)
        df['energy_needed_elec'] = self.f.needed_house_energy(df.temperature) / df.efficiency
        df['time'] = df.date.dt.time
        df['freezing_zone'] = df['temperature'] < c.freezing_zone

        return df

    def get_local_extrema(self, df):  # get local minima and maxima in function
        check_interval = 10

        df['local_min'] = df.iloc[argrelextrema(df.temperature.values, np.less_equal, order=check_interval)[0]][
            'temperature']
        df['local_max'] = df.iloc[argrelextrema(df.temperature.values, np.greater_equal, order=check_interval)[0]][
            'temperature']
        df['local_max_flow'] = df.iloc[argrelextrema(df.flow_temp.values, np.greater_equal, order=check_interval)[0]][
            'flow_temp']

        # delete doubled high and low points -> take only the later ones
        df.local_min[
            (pd.notna(df.local_min - df.local_min.shift(-1))) | (pd.notna(df.local_min - df.local_min.shift(-2))) | (
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

        if local_max_index[0] < local_min_index[
            0]:  # reassure that period counting starts with period_min -> heating period comes before deheating
            local_max_index = local_max_index[1:]

        for i, index in enumerate(local_min_index):
            df.loc[df.index > index, 'period_min'] = i

        for i, index in enumerate(local_max_index):
            df.loc[df.index >= index, 'period_max'] = i

        return df
