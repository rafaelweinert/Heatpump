import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.io as pio
pio.templates.default = "plotly"

def plot_temperature(df, extrema = False):
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=df.date, y=df.temperature, name='temperature'))
    fig.update_layout(dict1={'title': 'temperature'})
    if extrema:
        fig.add_trace(go.Scatter(x=df.date, y=df.local_min, mode='markers', name='local_minima'))
        fig.add_trace(go.Scatter(x=df.date, y=df.local_max, mode='markers', hovertext=df.period_max, name='local_maxima'))

    fig.show()

    #print(type(fig))
    return fig

def plot_heatmap(matrix):
    fig = px.imshow(matrix)

    return fig


def plot_optimization(df):
    fig = make_subplots(specs=[[{"secondary_y": True}]])
    fig.add_trace(go.Scatter(x=df.date, y=df.temperature, name='temperature', mode='lines'),
                  secondary_y=False)
    fig.add_trace(go.Scatter(x=df.date, y=df.P_e_opt, name='Power (e) optimized', mode='markers'),
                  secondary_y=True)
    fig.add_trace(go.Scatter(x=df.date, y=df.energy_needed_elec, name='Power (e) not optimized',
                             mode='markers'),
                  secondary_y=True)
    fig.add_trace(go.Scatter(x=df.date, y=df.Storage_vol, name='Storage Volume'),
                  secondary_y=True)
    fig.add_trace(go.Scatter(x=df.date, y=df.Preheat_temp, name='Preheat temperature'),
                  secondary_y=False)
    fig.add_trace(go.Scatter(x=df.date, y=df.flow_temp, name='Flow temperature'), secondary_y=False)
    fig.update_layout(title='Power (elec) for optimized and not optimized heating')
    fig.write_html('optimization_data_temp.html', auto_open=False)
