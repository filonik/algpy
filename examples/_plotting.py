"""Small Plotly slider and vector helpers shared by the notebooks."""
import numpy as np
import plotly.graph_objects as go


def slider_controls(figure, labels, prefix="t = "):
    """Select precomputed frames immediately, with no playback or transitions.

    Plotly uses the animate method to select a named frame. All durations are
    zero; redraw is needed for 3D surfaces and meshes, and for images, to
    update reliably.
    """
    if len(labels) != len(figure.frames):
        raise ValueError("Each frame needs one slider label")
    three_d = any(trace.type in {"scatter3d", "surface", "mesh3d"}
                  for trace in figure.data)
    redraw = three_d or any(trace.type == "image" for trace in figure.data)
    options = dict(frame=dict(duration=0, redraw=redraw),
                   transition=dict(duration=0), mode="immediate", fromcurrent=True)
    if not three_d:
        figure.update_xaxes(autorange=False)
        figure.update_yaxes(autorange=False)
    if three_d:
        figure.update_scenes(xaxis_autorange=False, yaxis_autorange=False,
                             zaxis_autorange=False, aspectmode="cube")
    figure.update_layout(
        template="plotly_white",
        uirevision="examples",
        margin=dict(l=30, r=30, t=90, b=90),
        updatemenus=[],
        sliders=[dict(
            active=0, x=0, len=1, xanchor="left", y=-.12, yanchor="middle",
            pad=dict(t=0, b=0, l=0, r=0),
            currentvalue=dict(visible=False), ticklen=0, minorticklen=0,
            tickcolor="rgba(0,0,0,0)",
            steps=[dict(label="", method="animate", args=[[frame.name], options],
                        value=f"{prefix}{label}")
                   for frame, label in zip(figure.frames, labels)],
        )],
    )
    return figure


def vector2(point, direction, name, color, length=0.7):
    end = np.asarray(point) + length*np.asarray(direction)
    return go.Scatter(
        x=[point[0], end[0]], y=[point[1], end[1]], mode="lines+markers",
        name=name, line=dict(color=color, width=3),
        marker=dict(color=color, size=[0, 10], symbol="arrow", angleref="previous"),
    )


def vector3(point, direction, name, color, length=0.5):
    end = np.asarray(point) + length*np.asarray(direction)
    return go.Scatter3d(
        x=[point[0], end[0]], y=[point[1], end[1]], z=[point[2], end[2]],
        mode="lines+markers", name=name, line=dict(color=color, width=7),
        marker=dict(color=color, size=[0, 4]),
    )
