import pandas as pd
import numpy as np
from geopy.geocoders import Nominatim
import time
import os

#station_data = pd.read_table(r'C:\Users\rafaelweinert\PycharmProjects\Wärmepumpe\data\TU_Stundenwerte_Beschreibung_Stationen.txt', sep='\t', encoding='windows-1252')
station_data = pd.read_csv(r'C:\Users\rafaelweinert\PycharmProjects\Wärmepumpe\data\Beschreibung_Stationen.txt', sep='\t', skiprows=[1])

data = pd.DataFrame()

for entry in station_data.iterrows():
    k = pd.Series(entry[1].values[0][:60].split(), dtype=float)
    data = data.append(k, ignore_index=True)

data.columns = ['Stations_id', 'von_datum', 'bis_datum', 'Stationshoehe', 'geoBreite', 'geoLaenge']
print(data)
#%%
type(data.geoLaenge[0])

#%%

address = 'Rintheimer Straße 84, 76131 Karlsruhe, Germany'
app = Nominatim(user_agent="user")
def get_latlon(address):
    try:
        location = app.geocode(address).raw
        return location
    except:
        print('Not found')
        time.sleep(1)
        return get_latlon(address)

location = get_latlon(address)
lat, lon = float(location['lat']), float(location['lon'])

def get_closest_station(lat, long):
    dist = np.sqrt(np.square(data.geoBreite - lat) + np.square(data.geoLaenge - long))
    min_dist = np.argsort(dist)

    return min_dist

closest_station_idx = get_closest_station(lat, lon)[0]
closest_station = int(data.loc[closest_station_idx, 'Stations_id'])

print('closest station: ', closest_station)


#%%
#assume data is downloaded already
def get_station_data(station_ids):
    dir = os.listdir(r"C:\Users\rafaelweinert\PycharmProjects\Wärmepumpe\data\weatherstationdata")
    for station_id in station_ids:
        for entry in dir:
            if str(station_id) in entry:
                return pd.read_table(r"C:\Users\rafaelweinert\PycharmProjects\Wärmepumpe\data\weatherstationdata" + '\\' + entry, sep=';')
    return 'No data found'

data_temp = get_station_data([5100])
