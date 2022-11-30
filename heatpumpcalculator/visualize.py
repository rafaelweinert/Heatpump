import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.io as pio
pio.templates.default = "plotly"

def plot_temperature(df, address, extrema = False):
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=df.date, y=df.temperature, name='Temperature in [°C]'))
    fig.update_layout(dict1={'title': f'Temperature of {address}'})
    if extrema:
        fig.add_trace(go.Scatter(x=df.date, y=df.local_min, mode='markers', name='local_minima'))
        fig.add_trace(go.Scatter(x=df.date, y=df.local_max, mode='markers', hovertext=df.period_max, name='local_maxima'))
        fig.add_trace(go.Scatter(x=df.date, y=df.flow_temp, name='Flow temperature in [°C]'))
    #print(type(fig))
    return fig

def plot_heatmap(matrix):
    fig = px.imshow(matrix)

    return fig


def plot_optimization(df):
    fig = make_subplots(specs=[[{"secondary_y": True}]])
    fig.add_trace(go.Scatter(x=df.date, y=df.temperature, name='Temperature in [°C]', mode='lines'),
                  secondary_y=False)
    fig.add_trace(go.Scatter(x=df.date, y=df.P_e_opt, name='Electric power optimized in [W]', mode='markers'),
                  secondary_y=True)
    fig.add_trace(go.Scatter(x=df.date, y=df.energy_needed_elec, name='Electric power not optimized in [W]',
                             mode='markers'),
                  secondary_y=True)
    fig.add_trace(go.Scatter(x=df.date, y=df.Storage_vol, name='Storage volume in [L]'),
                  secondary_y=True)
    fig.add_trace(go.Scatter(x=df.date, y=df.Preheat_temp, name='Preheating temperature in [°C]'),
                  secondary_y=False)
    fig.add_trace(go.Scatter(x=df.date, y=df.flow_temp, name='Flow temperature in [°C]'), secondary_y=False)
    #fig.update_layout(title='Power (elec) for optimized and not optimized heating')
    fig.update_layout(plot_bgcolor="white",
                      margin=dict(l=50, r=50, t=50, b=60), )
    fig.for_each_xaxis(lambda x: x.update(showgrid=True, gridwidth=0.5, gridcolor='rgb(242,242,242)'))
    fig.for_each_yaxis(lambda x: x.update(showgrid=True, gridwidth=0.5, gridcolor='rgb(242,242,242)'))
    return fig
