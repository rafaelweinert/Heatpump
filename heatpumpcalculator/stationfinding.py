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
        self.station_data = pd.read_csv(r'C:\Users\rafaelweinert\PycharmProjects\Wärmepumpe\data\Beschreibung_Stationen.txt', sep='\t', skiprows=[1])

        self.data = pd.DataFrame()

        for entry in self.station_data.iterrows():
            k = pd.Series(entry[1].values[0][:60].split(), dtype=float)
            self.data = self.data.append(k, ignore_index=True)

        self.data.columns = ['Stations_id', 'von_datum', 'bis_datum', 'Stationshoehe', 'geoBreite', 'geoLaenge']
        #print(self.data)

    def get_latlon(self, address):
        try:
            print('addresse: ', address)
            location = self.app.geocode(address).raw
            return location
        except:
            print('Address not found! Retry with a different address.')

            exit()

    def get_closest_station(self, lat, long):
        dist = np.sqrt(np.square(self.data.geoBreite - lat) + np.square(self.data.geoLaenge - long))
        min_dist = np.argsort(dist)

        return min_dist

    #assume data is downloaded already
    def get_station_data(self, station_ids):
        dir = os.listdir(r"C:\Users\rafaelweinert\PycharmProjects\Wärmepumpe\data\weatherstationdata")
        for station_id in station_ids:
            for entry in dir:
                if str(station_id) in entry:
                    return pd.read_table(r"C:\Users\rafaelweinert\PycharmProjects\Wärmepumpe\data\weatherstationdata" + '\\' + entry, sep=';')
        return 'No data found'

    def main(self):
        location = self.get_latlon(self.address)

        lat, lon = float(location['lat']), float(location['lon'])
        closest_station_idx = self.get_closest_station(lat, lon)[0]
        closest_station = int(self.data.loc[closest_station_idx, 'Stations_id'])

        data_temp = self.get_station_data([closest_station])

        return data_temp




