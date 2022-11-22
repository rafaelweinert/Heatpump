import constants as c

import numpy as np


class Functions:

    def __init__(self, flow_max, room_temp, out_min, return_flow, house_energy, storage_volume, storage_power, ny):

        self.flow_max = flow_max
        self.room_temp = room_temp
        self.out_min = out_min
        self.return_flow = return_flow
        self.house_energy = house_energy
        self.storage_volume = storage_volume
        self.storage_power = storage_power
        self.ny = ny

        self.coef_h = (self.storage_power / (c.h2o_energy * self.storage_volume)) / np.abs(c.t_env - c.t_0)

    def flowtemp(self, t):
        flow =  - (self.flow_max - self.room_temp) / (self.room_temp-self.out_min) *t + self.flow_max + self.out_min * (self.flow_max - self.room_temp) / (self.room_temp-self.out_min)
        flow[flow <= self.return_flow] = self.return_flow + 0.01
        return flow
    def needed_house_energy(self, t):
        return np.maximum(0, self.room_temp - t) * self.house_energy # returns needed energy to heat house in Watt

    def eta(self, zu, ab = None): #theoretischer Wirkungsgrad
        if ab is None:
            ab = self.flowtemp(zu)
        return np.maximum(0, np.minimum(300, self.ny * c.kelvin(ab) / (c.kelvin(ab)-c.kelvin(zu))))

    def eta_heating(self, zu, max_ab = None, min_ab = None,  ): #"gleitender" Wirkungsgrad -> da sich über temperaturerhöhung ändert
        try:
            if type(min_ab) == type(None):
                min_ab = self.return_flow
            if type(max_ab) == type(None):
                max_ab = self.flowtemp(zu)
        except:
            print("Error: function 'eta_heating'")
        return self.ny * (c.kelvin(zu) * np.log((c.kelvin(max_ab) - c.kelvin(zu)) / (c.kelvin(min_ab) - c.kelvin(zu))) + c.kelvin(max_ab) - c.kelvin(min_ab)) / (c.kelvin(max_ab) - c.kelvin(min_ab))

    def eta_heating2(self, zu, max_ab = None, min_ab=None, ):
        try:
            if type(min_ab) == type(None):
                min_ab = self.return_flow
            if type(max_ab) == type(None):
                max_ab = self.flowtemp(zu)
        except:
            print("Error: function 'eta_heating'")

        return self.ny * ((c.kelvin(max_ab) - c.kelvin(min_ab)) / (c.kelvin(max_ab) - c.kelvin(min_ab) - c.kelvin(zu) * np.log(c.kelvin(max_ab) / c.kelvin(min_ab))))

    def eta_heating3(self, zu, max_ab = None, min_ab=None, ): # uses an average of maximum and minimum temperature
        try:
            if type(min_ab) == type(None):
                min_ab = self.return_flow
            if type(max_ab) == type(None):
                max_ab = self.flowtemp(zu)
        except:
            print("Error: function 'eta_heating'")

        ab_new = (max_ab + min_ab) / 2
        return np.maximum(0, self.ny * c.kelvin(ab_new) / (c.kelvin(ab_new)-c.kelvin(zu)))

    def eta_freezing_adj(self, zu, max_ab = None, min_ab = None, heating =  False): # if heating = True, then use eta_heating
        try:
            if type(min_ab) == type(None):
                min_ab = self.return_flow
            if type(max_ab) == type(None):
                max_ab = self.flowtemp(zu)
        except:
            print("Error: function 'eta_heating'")

        if heating == True:
            e = self.eta_heating3(zu, max_ab, min_ab)
        else:
            e = self.eta(zu, max_ab)
        #print(type(eta))
        if type(e) == np.float64:
            if zu < 4:
                e = e / c.freezing_factor
        else:
            e[zu < 4] = e[zu < 4] / 1.2

        return e

    def newton_cooling(self, t):
        return c.t_env + (c.t_0 - c.t_env) * np.exp(-self.coef_h * t)  # = temperature at time t

    def newton_cooling_t_0(self, temp_t, t=7 * 3600):  # enter in Kelvin
        # return ((temp_t - t_env) * np.exp((coef_h * t).T)).T + t_env #returns t_0 for fix temperature after time t
        return (temp_t - c.t_env) * np.exp((self.coef_h * t)) + c.t_env  # returns t_0 for fix temperature after time t

    def newton_cooling_t(self, temp_t):  # , coef_h = coef_h):
        return -np.log((temp_t - c.t_env) / (
                    c.t_0 - c.t_env)) / self.coef_h  # returns how much time it takes for t_0 to cool down to wished temperature temp_t

    def storage_power_temp(self, temp):
        return -c.h2o_energy * self.storage_volume * self.newton_cooling_diff(
            self.newton_cooling_t(temp))  # returns how much power in Watt is needed to keep temperature at temp level

    def get_d_t(self, ):
        return self.storage_power / (c.h2o_energy * self.storage_volume)  # returns heat loss in kelvin per second

    def get_coef_h(self, d_t):
        return d_t / np.abs(c.t_env - c.t_0)  # returns coefficient at which heat decreases in storage

    def get_coef_h_adj(self, temp, flow, t):
        return - np.log((flow - c.t_env) / (
                    temp - c.t_env)) / t  # returns adjusted coef_h for given flow, overheat_temp and time after which flow should be reached

    def storage_power_temp_adj(self, temp, flow, t):
        return c.h2o_energy * self.storage_volume * (self.coef_h - self.get_coef_h_adj(temp, flow, t)) * (
                    temp - c.t_env)  # returns needed power to keep temperature in storage and decrease at given rate

    def get_needed_storage_vol(self, temp_diff, kWh):
        return (kWh * 3600000) / (
                    temp_diff * c.h2o_energy)  # returns waterstorage volume to store #kWh with given temp difference

    #d_t = get_d_t()
    #coef_h = get_coef_h(d_t)