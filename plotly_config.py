import plotly.io as pio

# Define colours
DEFAULT_COLOURS = {
    "background": "#C9D0D6",
    "paper": "yellow",
    "line": "#1A2228",
    "accent": "#B03060"
}

# Layout defaults
DEFAULT_LAYOUT = dict(
    plot_bgcolor=DEFAULT_COLOURS["background"],
    paper_bgcolor=DEFAULT_COLOURS["paper"],
    font=dict(family="Arial", size=14, color="#1A2228"),
    margin=dict(l=80, r=50, t=60, b=60),
    showlegend=False
)

# Axis defaults
DEFAULT_AXES = dict(
    showline=True,
    linecolor="black",
    mirror=True,
    showgrid=False,
    zeroline=False
)

# Trace defaults
DEFAULT_TRACES = dict(
    scatter=[dict(line=dict(color=DEFAULT_COLOURS["line"], width=1.2))]
)


def init_plotly():
    """
    Registers and activates a consistent Plotly style for AIDDTools figures.
    Call this once before creating any figures.
    """
    pio.templates["aiddtools"] = dict(
        layout=DEFAULT_LAYOUT,
        data=DEFAULT_TRACES
    )
    pio.templates.default = "aiddtools"
    print("[plotly_config] Plotly settings applied (template: 'aiddtools')")


# Optional: automatically apply settings on import
init_plotly()
