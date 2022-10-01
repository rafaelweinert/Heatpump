import constants as c
import functions

import numpy as np
import pandas as pd
from scipy.optimize import Bounds, minimize


class Optimizer:

    def __init__(self, functions, max_heat_cap, hp_power_el_opt, return_flow, storage_volume):
        self.f = functions
        self.max_heat_cap = max_heat_cap
        self.hp_power_el_opt = hp_power_el_opt
        self.return_flow = return_flow
        self.storage_volume = storage_volume


# optimize over flowtemp
# TODO divide by period-length (=1, when hourly temperatures)

    def heating_optimization(self, df, period_lim=None, period=None, **kwargs):
        if type(period) != type(None):
            df_per = df[(df.period_max == period) | (df.period_min == period)]
        else:
            df_per = df.copy()

        per = df_per.shape[0]  # number of hours in period
        M_timedelta = np.empty((per, per))

        for i in range(per):  # build timedelta matrix
            M_timedelta[i, range(0, i + 1)] = np.flip(np.linspace(0, i, i + 1))
        M_timedelta = M_timedelta.T
        flow = df_per.flow_temp.values

        M_flow = c.celsius(self.f.newton_cooling_t_0(c.kelvin(flow), M_timedelta * 3600))
        M_flow = np.triu(
            M_flow)  # needed heating temperature when heating in period i and using the heat in period j (i X j - Matrix)
        M_temperature = np.tile(df_per.temperature.values, (len(df_per.temperature.values), 1)).T

        def target_fun(fl, optimizing=True, constraint=False):

            M_eta = self.f.eta_freezing_adj(zu=M_temperature, max_ab=M_flow)  # , min_ab=M_flow)

            if type(period_lim) != type(None):  # period_lim also means storage limitation
                per_lim = np.minimum(period_lim, per)
                M_flow_max = c.celsius(self.f.newton_cooling_t_0(c.kelvin(fl), t=M_timedelta * 3600))
                M_eta_heating = self.f.eta_freezing_adj(zu=M_temperature, max_ab=M_flow_max, heating=c.h)
                M_eta_heating = np.triu(M_eta_heating, k=1)
                M_eta = np.tril(M_eta)

                M_eta = M_eta_heating + M_eta
                M_eta = np.triu(M_eta)  # upper diagonal
                M_eta = np.tril(M_eta, k=per_lim)

            else:
                M_eta = np.triu(M_eta)  # upper diagonal
            P_w_sum = df_per.energy_needed_water.values  # needed energy in water for each hour
            P_e_max = np.full((per,), self.hp_power_el_opt)  # maximal available energy (electricity) per hour -> e.g. 2500W
            P_e_available = P_e_max.copy()  # at each period still available 'capacity' in electricity for future hours

            P_e = np.zeros((per, per))  # Matrix with used energy (e) from each period for each period that is optimized
            P_e_a = np.empty((per, per))  # Matrix to later calculate P_e

            for j in range(per):
                P_w_j = P_w_sum[j]  # energy needed (w) in hour j

                M_eta_j = M_eta[:, j]  # efficiencies in each hour
                eta_sorted_j_index = np.flip(np.argsort(M_eta_j))  # get indexes of highest values in M_eta_j
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
            np.fill_diagonal(diagonal, self.hp_power_el_opt)

            P_e = np.abs(diagonal - P_e)
            P_w = (P_e * M_eta)

            if optimizing == False:
                df_per['P_e_opt'] = P_e.sum(1)
                return M_eta, P_e, P_w, df_per

            if constraint == True:
                return P_w

            return P_e.sum()

        def constraint_fun(fl):

            P_w = target_fun(fl, constraint=True)
            stor_vol = np.zeros((P_w.shape[0],))

            for i in range(P_w.shape[0] - 1):
                stor_vol[i] = self.f.get_needed_storage_vol(fl - self.return_flow, P_w[:(i + 1), (i + 1):].sum() / 1000)
            max_stor = max(stor_vol)

            return self.storage_volume - max_stor

        if 'optimize' in kwargs:  # same as heating_opt
            if kwargs['optimize'] == False:
                if 'flow' in kwargs:
                    M_eta_opt, P_e_opt, P_w_opt, df_per_opt = target_fun(kwargs['flow'], optimizing=False)
                    stor_vol = np.zeros((P_w_opt.shape[0],))

                    for i in range(P_w_opt.shape[0] - 1):
                        stor_vol[i] = self.f.get_needed_storage_vol(kwargs['flow'] - self.return_flow,
                                                             P_w_opt[:(i + 1), (i + 1):].sum() / 1000)
                    max_stor = max(stor_vol)

                    return M_eta_opt, P_e_opt, P_w_opt, df_per_opt
                else:
                    return 'Must enter a temperature'

        # optimization problem
        lb = max(flow) + 0.001
        ub = self.max_heat_cap
        bounds = Bounds(lb=lb, ub=ub)


        min = minimize(target_fun, x0=(lb + 1,), bounds=bounds, constraints={'type': 'ineq',
                                                                             'fun': constraint_fun})

        M_eta_opt, P_e_opt, P_w_opt, df_per_opt = target_fun(min.x[0],
                                                             optimizing=False)  # get values for optimal heating temperature

        return min, M_eta_opt, P_e_opt, P_w_opt, df_per_opt


    def optimize_periods(self, df, max_periods=36):
        index = 0
        df_new = pd.DataFrame()
        M_eta_new = np.zeros((df.shape[0], df.shape[0]))
        P_e_new = np.zeros((df.shape[0], df.shape[0]))
        P_w_new = np.zeros((df.shape[0], df.shape[0]))
        for p in df.period_max.dropna().unique():
            print('period: ', int(p), f' / {int(max(df.period_max.dropna()))}', end='\r')

            min_df = pd.DataFrame(columns=['period_lim', 'min_x', 'min_fun'])

            for i in range(max_periods):
                minimum, M_eta, P_e, P_w, df_per = self.heating_optimization(df.loc[index:, :], period=p, period_lim=i)
                if minimum.success == True:
                    min_df = min_df.append({'period_lim': i,
                                            'min_x': minimum.x[0],
                                            'min_fun': minimum.fun}, ignore_index=True)

            temp = min_df.loc[np.argmin(min_df.min_fun)].min_x
            per_lim = min_df.loc[np.argmin(min_df.min_fun)].period_lim

            M_eta_opt, P_e_opt, P_w_opt, df_per_opt = self.heating_optimization(df.loc[index:, :], period=p,
                                                                           period_lim=per_lim, optimize=False,
                                                                           flow=temp)

            stor_vol = np.zeros((P_w_opt.shape[0],))
            for i in range(P_w_opt.shape[0] - 1):
                stor_vol[i] = self.f.get_needed_storage_vol(temp - self.return_flow, P_w_opt[:(i + 1), (i + 1):].sum() / 1000)

            df_per_opt['Storage_vol'] = stor_vol
            df_per_opt['Preheat_temp'] = temp
            df_per_opt['Period_lim'] = per_lim

            df_new = pd.concat([df_new, df_per_opt])
            M_eta_new[index:index + M_eta_opt.shape[0], index:index + M_eta_opt.shape[1]] = M_eta_opt
            P_e_new[index:index + P_e_opt.shape[0], index:index + P_e_opt.shape[1]] = P_e_opt
            P_w_new[index:index + P_w_opt.shape[0], index:index + P_w_opt.shape[1]] = P_w_opt

            index = df_per_opt[df_per_opt.Storage_vol == 0].Storage_vol.index[-1] + 1
        print()
        return M_eta_new, P_e_new, P_w_new, df_new

