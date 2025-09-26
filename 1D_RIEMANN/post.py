import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.io as pio
pio.templates.default = "seaborn"
import h5py

def read_frame(filename, var1, var2=None):
    mesh         = h5py.File(filename, 'r')['mesh']
    points       = mesh['points']
    connectivity = mesh['connectivity']

    segments          = np.zeros((connectivity.shape[0], 2, 2))
    segments[:, :, 0] = points[:][connectivity[:]][:, :, 0]

    centers = 0.5*(segments[:, 0, 0] + segments[:, 1, 0])

    field1 = mesh['fields'][var1][:]
    if var2 is not None:
        field2 = mesh['fields'][var2][:]

    index = np.argsort(centers)

    x = centers[index]
    field1 = field1[index]
    if var2 is not None:
        field2 = field2[index]

    if var2 is not None:
        return x, field1, field2

    return x, field1

def plot_1D_Riemann_results(filename_base, filename_time, field_name_1, field_name_2=None):
    frames = []
    steps = []

    t = np.loadtxt(filename_time)

    for i, ti in enumerate(t):
        name = str(ti)
        filename = filename_base + f"{i:05d}" + ".h5"

        if field_name_2 is not None:
            x, field_1, field_2 = read_frame(filename, field_name_1, field_name_2)
            frames.append(go.Frame(data=[go.Scatter(x=x, y=field_1, name=field_name_1), \
                                         go.Scatter(x=x, y=field_2, name=field_name_2)], \
                                   name=name))
        else:
            x, field_1 = read_frame(filename, field_name_1, field_name_2)
            frames.append(go.Frame(data=[go.Scatter(x=x, y=field_1, name=field_name_1)], \
                                   name=name))

        step = dict(method="animate", label = f"{ti:.2e}", \
                    args=[[name], {"frame": {"redraw": False, "duration": 300}, \
                    "mode": "immediate", "transition": {"duration": 10}}])
        steps.append(step)

    sliders = [dict(currentvalue={'prefix': 't = '}, active=i, visible=True, steps=steps)]

    fig = make_subplots(specs=[[{"secondary_y": True}]])
    fig.add_trace(go.Scatter(x=frames[i].data[0].x, y=frames[i].data[0].y, name=field_name_1,line=dict(color="red")), secondary_y=False)
    if field_name_2 is not None:
        fig.add_trace(go.Scatter(x=frames[i].data[1].x, y=frames[i].data[1].y, name=field_name_2,line=dict(color="blue")), secondary_y=True)
    fig.update_layout(sliders=sliders, height=500,paper_bgcolor='rgba(0,0,0,0)')
    fig.update(frames=frames)
    fig.update_yaxes(title_text=field_name_1, secondary_y=False)
    if field_name_2 is not None:
        fig.update_yaxes(title_text=field_name_2, secondary_y=True)

    fig.show()
