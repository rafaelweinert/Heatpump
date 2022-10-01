from geopy.geocoders import Nominatim
import meteostat
import requests
import pandas as pd
import json

def get_weather_forecast(address, days=4, key = '676ea42d0beb4636b96160506221607'):

    # Get coordinates for adress
    geolocator = Nominatim(user_agent="user")
    location = geolocator.geocode(address)

    # Get nearby weather stations
    stations = meteostat.Stations()
    stations = stations.nearby(location.latitude, location.longitude)
    station = stations.fetch(1)

    # Print DataFrame
    print(station)
    weather_data = requests.get(f'http://api.weatherapi.com/v1/forecast.json?key={key}&q={address}&days={days}&aqi=no&alerts=no')

    f = pd.DataFrame.from_dict(json.loads(weather_data.content.decode())['forecast']['forecastday'])
    forecast = pd.json_normalize(pd.json_normalize(f['hour']).stack())

    forecast.time = pd.to_datetime(forecast.time, )
    forecast = forecast[['time', 'temp_c']]

    forecast.columns = ['date', 'temperature']
    return forecast