import numpy as np


kelvin = lambda x: x+273.15 # convert Celsius to Kelvin
celsius = lambda x: x-273.15 # convert Kelvin to Celsius
kWh = lambda x: x / 3600000 # convert Joule to kWh
Joule = lambda x: x * 3600000 # convert kWh to Joule

freezing_zone = 4
freezing_factor = 1.2# factor which the efficiency coefficient is divided by in case of being in freezing zone

#20 and 60 after some DIN-Norms
t_env = kelvin(20) #environmental temperature where water storage is located
t_0 = kelvin(65) #temperature at which water storage power is measured

h2o_energy = 4200 #Joule / (kilogram * Kelvin)

# Loefflers heating formula
h = False
