import pandas as pd
import numpy as np
from geopy.geocoders import Nominatim
import time
import os

import warnings
warnings.filterwarnings('ignore')


class StationFinder:

    def __init__(self, address):
        self.address = address
        self.app = Nominatim(user_agent="user")
    #station_data = pd.read_table(r'C:\Users\rafaelweinert\PycharmProjects\Wärmepumpe\data\TU_Stundenwerte_Beschreibung_Stationen.txt', sep='\t', encoding='windows-1252')
        self.station_data = pd.read_csv(r'C:\Users\rafaelweinert\PycharmProjects\Heatpump\data\Beschreibung_Stationen.txt', sep='\t', skiprows=[1])

        self.data = pd.DataFrame()

        for entry in self.station_data.iterrows():
            k = pd.Series(entry[1].values[0][:60].split(), dtype=float)
            self.data = self.data.append(k, ignore_index=True)

        self.data.columns = ['Stations_id', 'von_datum', 'bis_datum', 'Stationshoehe', 'geoBreite', 'geoLaenge']
        #print(self.data)

    def get_latlon(self, address):
        try:
            location = self.app.geocode(address).raw
            #print('location: ', location)
            return location
        except:
            print('Address not found! Retry with a different address.')

            exit()

    def get_closest_station(self, lat, long):
        dist = np.sqrt(np.square(self.data.geoBreite - lat) + np.square(self.data.geoLaenge - long))
        min_dist = np.argsort(dist)

        return min_dist

    def get_closest_try(self, lat, long):
        dist = np.sqrt(np.square(self.try_data.geoBreite - lat) + np.square(self.try_data.geoLaenge - long))
        min_dist = np.argsort(dist)

        return min_dist

    #assume data is downloaded already
    def get_station_data(self, station_ids):
        dir = os.listdir(r"C:\Users\rafaelweinert\PycharmProjects\Heatpump\data\weatherstationdata")
        for station_id in station_ids:
            for entry in dir:
                if str(station_id) in entry:
                    return pd.read_table(r"C:\Users\rafaelweinert\PycharmProjects\Heatpump\data\weatherstationdata" + '\\' + entry, sep=';')
        return 'No data found'

    def get_try_data(self, coordinates):

        lat = coordinates[0]
        lon = coordinates[1]

        dir = os.listdir(r"C:\Users\rafaelweinert\PycharmProjects\Heatpump\heatpumpcalculator\data")
        dir = pd.Series(dir)

        try_lat = dir.str.split('_', expand=True)[1].str[:6].astype(float) / 10000
        try_lon = dir.str.split('_', expand=True)[1].str[6:].astype(float) / 10000

        dist = np.sqrt(np.square(try_lat - lat) + np.square(try_lon - lon))
        min_dist = np.argsort(dist)[0]

        lat_min_dist = try_lat[min_dist]
        if lat_min_dist < 10:
            lat_min_dist = '0' + str(int(lat_min_dist * 10000))
        else:
            lat_min_dist = str(int(lat_min_dist * 10000))
        lon_min_dist = try_lon[min_dist]
        if lon_min_dist < 10:
            lon_min_dist = '0' + str(int(lon_min_dist * 10000))
        else:
            lon_min_dist = str(int(lon_min_dist * 10000))

        min_dist_folder = 'TRY_' + lat_min_dist + lon_min_dist
        min_dist_entry = 'TRY2015_' + lat_min_dist + lon_min_dist + '_Jahr.dat'

        #print('min_dist_entry:', min_dist_entry)

        data = pd.read_table(
            r"C:\Users\rafaelweinert\PycharmProjects\Heatpump\heatpumpcalculator\data" + '\\' + min_dist_folder + '\\' + min_dist_entry,
            skiprows=30)
        col = data.columns[0]
        data = data[col].str.split(expand=True)
        data.columns = col.split()
        data = data.loc[1:]

        data.MM = data.MM.str.zfill(2)
        data.DD = data.DD.str.zfill(2)
        data.HH = data.HH.str.zfill(2)
        index_24 = data.HH[data.HH == '24'].index
        data.HH[data.HH == '24'] = '23'
        data['date'] = pd.to_datetime(
            '2022-' + data.MM.astype(str) + '-' + data.DD.astype(str) + ' ' + data.HH.astype(str) + ':00',
            format='%Y-%m-%d %H:%M')

        data.loc[index_24, 'date'] = data.loc[index_24].date + pd.to_timedelta(1, 'h')
        data.loc[data.date < pd.to_datetime('2022-07-01'), 'date'] = data.loc[data.date < pd.to_datetime('2022-07-01')].date + pd.to_timedelta(365, 'D')
        data.rename({'t': 'temperature'}, axis=1, inplace=True)
        data = data[['date', 'temperature']]
        data.temperature = data.temperature.astype(float)
        data = data.sort_values(by='date')

        data = data[(data.date >= pd.to_datetime('2022-10-01')) & (data.date <= pd.to_datetime('2023-03-31'))]
        data.reset_index(inplace=True, drop=True)

        return data


    def main(self, try_data = True):
        location = self.get_latlon(self.address)

        lat, lon = float(location['lat']), float(location['lon'])

        #print('lat:', lat)
        #print('lon:', lon)

        if try_data:
            data_temp = self.get_try_data((lat, lon))

        else:
            closest_station_idx = self.get_closest_station(lat, lon)[0]
            closest_station = int(self.data.loc[closest_station_idx, 'Stations_id'])

            data_temp = self.get_station_data([closest_station])

        return data_temp







