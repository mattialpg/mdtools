# Colour palettes for scientific figures and data visualisation
# Sources:
# - https://www.molecularecologist.com/2020/04/23/simple-tools-for-mastering-color-in-scientific-figures/
# - https://github.com/OrdnanceSurvey/GeoDataViz-Toolkit/blob/master/Colours/GDV-colour-palettes.json


palette_1 = ["#28353E", "#44576C", "#768B96", "#AAC7D6", "#DFEBF7", "#E7E7E7"]
palette_2 = ["#2A1D1B", "#861D1C", "#768B96", "#F5B34F", "#ECCEB6", "#FFFFFF"]
palette_3 = ["#290909", "#57473A", "#6F8758", "#74A1B1", "#D0D5CF", "#F0E8E2"]
palette_4 = ["#392D26", "#A58376", "#CCAC8E", "#D1C8BD", "#EBE3DC", "#F2ECE7"]
palette_5 = ['#D3D3D3', '#9EC9E2', '#3082BD', '#0D4A70', '#002147']
palette_6 = ["#0D4A70", "#3082BD", "#9EC9E2", "#E4F1F7"]
palette_7 = ["#B3B3B3", "#E6B8B8", "#CB3658"]
palette_8 = ['#434e56', '#CF597E']

squirtle = ['#94D3CB', '#4C8984', '#1C4340', '#AB4E3D', '#702D24']
mewtwo =   ['#792670', '#BF7B89', '#3D1228', '#DB67B5', '#BB8CBF', '#974B59']
pikachu =  ['#F5D43F', '#EFA042', '#E46D36', '#B84244', '#5D1115']
chopper =  ['#6E0F12', '#E46650', '#9C5A4C', '#E6B892', '#8A4A6A']
mihawk =   ['#1F2430', '#3A3F4B', '#9E2F2F', '#D65A3A', '#F08A3C']

pokemon_pantone = [
    '#D3C4D7', '#6E4F3A', '#38509C', '#53295C', '#F2EADF', '#02B0D3',
    '#F5CAC8', '#FDD87A', '#DE8F65', '#8EC696', '#BE1B32', '#2E2B2E', '#8E908F']

geodataviz_diverging_1 = [
    "#045275", "#089099", "#39B185", "#7CCBA2", "#FCDE9C", "#EEB479",
    "#E88471", "#CF597E", "#A33A6F", "#7E2A5C", "#5F1848"]

blacks = [
    "#232325", "#1A2228", "#08100C", "#23262A", "#101923", "#161110",
    "#21242A", "#100C08", "#333333", "#2B2E27", "#2D383A", "#36454F",
    "#2A2A2A", "#343434", "#1F262A", "#16161D"]  # 16161D is the Eigengrau!

traffic_lights = {'red': '#d10808', 'yellow': '#ffc32d', 'green': '#009b73'}

all_palettes = {
    "palette_1": palette_1,
    "palette_2": palette_2,
    "palette_3": palette_3,
    "palette_4": palette_4,
    "geodataviz_diverging_1": geodataviz_diverging_1,
}


def discrete_colorscale(scale_name, tick_labels):
    n_levels = len(tick_labels)
    scale = all_palettes.get(scale_name)[:n_levels]
    if scale is None:
        raise KeyError(f"Palette '{scale_name}' not found.")

    # Define integer bin boundaries for normalisation (0 to n_levels)
    boundaries = list(range(n_levels + 1))
    normalised = [(v - boundaries[0]) / (boundaries[-1] - boundaries[0]) for v in boundaries]

    # Build discrete colourscale
    colorscale = []
    for i in range(n_levels):
        colorscale.extend([[normalised[i], scale[i]], [normalised[i + 1], scale[i]]])

    # Tick positions: centre of each bin (in data units, not normalised)
    tickvals = [i + 0.5 for i in range(n_levels)]

    # Tick labels: use provided names or fallback to numbers
    if len(tick_labels) != n_levels:
        raise ValueError("tick_labels must have the same length as n_levels.")
    ticktext = tick_labels

    colorbar = dict(
        thickness=25,
        tickvals=tickvals,
        ticktext=ticktext,
        len=0.757,
        x=0.94,
        y=1.027,
        yanchor="top")

    return colorscale, colorbar
