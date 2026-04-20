import warnings
warnings.filterwarnings('ignore')

import os, sys, glob
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
from tqdm import tqdm
# import mdtraj as md
import MDAnalysis as mda

import plotly.graph_objects as go
import plotly.io as pio
from plotly.subplots import make_subplots

import tools.utils as utils
import visualisation.palettes as palettes

sys.path.append('/home/mattia/aiddtools')
import interaction_tools as inttools


class OccupancyAnalysis:
    def __init__(self, trj_name, configs):
        self.configs = configs
        self.trj_name = trj_name
        self.workdir = Path(self.configs['workdir'])
        self.frame_dir = self.workdir / 'complex.frames'
        self.int_types = ['hydrophobic', 'hbond', 'waterbridge',
            'saltbridge', 'pistacking', 'pication', 'halogen', 'metal']

        cwd = Path.cwd()
        if cwd != self.workdir:
            print('\033[1;31m\nWorkdir in config.yaml is different from current directory!\n\033[0m')
            input("Press Enter to continue...\n")

        print("Loading trajectory...")
        self.u = mda.Universe(f"{self.trj_name}.gro", f"{self.trj_name}.xtc")


    def extract_frames(self):
        non_water = self.u.select_atoms("not water and not resname K and not resname CL")

        # Build ligand selection from config
        lig_resnames = [d.get("md_id", d.get("resname", "LIG")) for d in self.configs.get("ligands", [])]
        lig_sel = " or ".join([f"resname {r}" for r in lig_resnames])

        # Select water shell
        water = self.u.select_atoms("resname SOL or resname HOH or resname WAT")
        cutoff_A = float(self.configs.get("hydration_shell_cutoff_A", 5))
        shell = water.select_atoms(f"byres around {cutoff_A} global ({lig_sel})",
            updating=True) if lig_sel else self.u.atoms[[]]

        total_frames = min(len(self.u.trajectory), 1001)
        frame_indices = np.linspace(0, len(self.u.trajectory) - 1, total_frames, dtype=int)
        for i, fi in enumerate(tqdm(frame_indices, desc="Exporting MD frames", unit="frame")):
            self.u.trajectory[fi]
            export_ag = (non_water + shell).unique
            export_ag.write(f"{self.frame_dir}/frame_{i+1:04d}.pdb")


    def pivot_df(self, df_lig):
        # Preserve frame order from the source table
        int_priority = {name: idx for idx, name in enumerate(self.int_types)}
        all_pdb = pd.Index(pd.unique(df_lig['PDB']), name='PDB')
        df_lig['RESID'] = df_lig['RESID'].replace({'NA': np.nan, 'NANA': np.nan})
        df_lig['INT_TYPE'] = df_lig['INT_TYPE'].where(
            df_lig['INT_TYPE'].isin(self.int_types), np.nan)

        # # Flag cells with more than one interaction type in the same frame/residue.
        # valid_rows = df_lig['INT_TYPE'].notna()
        # if valid_rows.any():
        #     ntypes = df_lig.loc[valid_rows].groupby(['PDB', 'RESID'])['INT_TYPE'].transform('nunique')
        #     multi_idx = ntypes[ntypes > 1].index
        #     df_lig.loc[multi_idx, 'INT_TYPE'] = 'multiple'

        df_lig['INT_PRIORITY'] = df_lig['INT_TYPE'].map(int_priority)
        df_lig = df_lig.sort_values('INT_PRIORITY')

        # Keep one row per frame even when there is no valid interaction
        df_pivot = df_lig.pivot_table(index='PDB', columns='RESID',
            values='INT_TYPE', aggfunc='first')
        df_pivot = df_pivot.reindex(all_pdb)
        df_pivot = df_pivot.replace(int_priority)

        if df_pivot.shape[1] > 0:
            occupancy = df_pivot.notna().sum() / df_pivot.shape[0]
            kept = occupancy[occupancy >= 0.05].index
            if len(kept) == 0:
                kept = occupancy.sort_values(ascending=False).index
            df_pivot = df_pivot[kept]

            corr = df_pivot.fillna(0).corr(method='pearson')
            ordered = corr.mean().sort_values(ascending=False).index
            df_pivot = df_pivot[ordered]

            # Keep 12 most populated residues
            top_populated = df_pivot.notna().sum().sort_values(
                ascending=False).head(12).index
            df_pivot = df_pivot[top_populated]
        
        return df_pivot


    def make_occupancy_trace(self, df_pivot):
        dcolorsc, colorbar_cfg = palettes.discrete_colorscale(
            'geodataviz_diverging_1', self.int_types)
        colorbar_cfg = dict(colorbar_cfg, x=1.0, y=0.5,
            xanchor='left', yanchor='middle', len=1.048)
        
        # Render missing values as grey
        left = abs(-0.01) / (len(self.int_types) + abs(-0.01))
        dcolorsc = [[left + x * (1.0 - left), c] for x, c in dcolorsc]
        dcolorsc = [[0.0, '#C9D0D6'], [left, '#C9D0D6']] + dcolorsc

        # Build time vector scaled to number of displayed frames
        n_frames = len(df_pivot.index)
        total_time_ns = self.u.trajectory[-1].time / 1000
        time_ns = np.linspace(0, total_time_ns, n_frames)
        dt = time_ns[1] - time_ns[0] if len(time_ns) > 1 else 1.0
        x_range = [time_ns[0] - dt/2, time_ns[-1] + dt/2]

        if df_pivot.shape[1] == 0:
            y_vals = np.array(['No interactions'])
            z_vals = np.full((1, n_frames), np.nan)
        else:
            y_vals = df_pivot.columns.to_numpy()
            z_vals = df_pivot.to_numpy(dtype=float).T

        customdata = np.full(z_vals.shape, 'NA', dtype=object)
        valid = ~np.isnan(z_vals)
        if valid.any():
            int_labels = np.array(self.int_types, dtype=object)
            customdata[valid] = int_labels[z_vals[valid].astype(int)]
        z_vals = np.where(valid, z_vals + 1e-6, -0.01)

        trace = go.Heatmap(
            x=time_ns, 
            y=y_vals,
            z=z_vals,
            colorscale=dcolorsc,
            zsmooth=False,
            zmin=-0.01, zmax=len(self.int_types),
            colorbar=colorbar_cfg, name='Occupancy',
            customdata=customdata,
            hovertemplate="Time=%{x:.2f} ns<br>Residue=%{y}<br>Interaction=%{customdata}<extra></extra>")
        
        trace.meta = {
            'layout': dict(
                height=500, width=1100,
                margin=dict(l=100, r=50, t=65, b=0),
                font=dict(color=palettes.blacks[0]),
                title=dict(
                    text=f"Interaction occupancy ({self.pdb_id}-{self.ligand_resname})",
                    x=0.5, y=0.95, xanchor='center', yanchor='top',
                    font=dict(size=20, weight='bold')),
                plot_bgcolor='#C9D0D6',
                hovermode='closest'),
            'xaxis_params': dict(showgrid=False, showticklabels=False, range=x_range),
            'yaxis_params': dict(
                title_text='<b>Residue</b>', tickmode='linear',
                showgrid=False),}
        return trace


    def _read_xvg_series(self, xvg_path, y_col):
        xvg_path = Path(xvg_path)
        rows = []
        with xvg_path.open() as handle:
            for line in handle:
                stripped = line.strip()
                if not stripped or stripped.startswith(('@', '#')):
                    continue
                fields = stripped.split()
                if len(fields) < 2:
                    continue
                try:
                    rows.append((float(fields[0]), float(fields[1])))
                except ValueError:
                    continue

        series_df = pd.DataFrame(rows, columns=['time_ns', y_col])
        if series_df.empty:
            raise ValueError(f'No data found in {xvg_path}')

        max_frames = 1001
        indices = np.linspace(0, len(series_df) - 1, max_frames, dtype=int) \
            if len(series_df) >= max_frames else np.arange(len(series_df))
        return series_df.iloc[indices]


    def _make_line_trace(self, series_df, y_col, color, width, name, hovertemplate):
        return go.Scatter(
            x=series_df['time_ns'],
            y=series_df[y_col],
            mode='lines',
            line=dict(color=color, width=width),
            name=name,
            hovertemplate=hovertemplate)


    def make_rmsd_trace(self):
        rmsd_df = self._read_xvg_series('trj_rmsd_lig.xvg', 'rmsd_nm')
        trace = self._make_line_trace(rmsd_df, 'rmsd_nm', palettes.blacks[0],
            1.3, 'RMSD (nm)', 'RMSD=%{y:.3f} nm<extra></extra>')

        trace.meta = {
            'layout': dict(
                height=250, width=1040,
                margin=dict(l=100, r=50, t=10, b=65),
                font=dict(color=palettes.blacks[0]),
                xaxis_title='<b>Time (ns)</b>',
                plot_bgcolor='#FFFFFF',
                showlegend=False,
                hovermode='x unified',
                hoverlabel=dict(
                    bgcolor='#444444',
                    font=dict(color='#FFFFFF', size=13),
                    bordercolor='#FFFFFF')),
            'xaxis_params': dict(
                showgrid=False,
                showline=True,
                linewidth=2,
                linecolor=palettes.blacks[0],
                ticks='outside',
                ticklen=4,
                tickwidth=1,
                tickcolor=palettes.blacks[0],
                unifiedhovertitle=dict(text='Time=%{x:.2f} ns')),
            'yaxis_params': dict(
                title=dict(text='<b>RMSD (nm)</b>'),
                gridcolor='lightgrey', showgrid=True,
                zeroline=True, zerolinecolor='lightgrey',
                showticklabels=True)}
        return trace


    def make_distance_trace(self):
        dist_df = self._read_xvg_series('trj_dist_lig.xvg', 'dist_nm')
        trace = self._make_line_trace(dist_df, 'dist_nm', '#B03060',
            1.1, 'Distance (nm)', 'Distance=%{y:.3f} nm<extra></extra>')
        
        trace.meta = {
            'axis_id': 'y2',
            'yaxis_params': dict(
                title=dict(text='<b>Distance (nm)</b>',
                    font=dict(color='#B03060')),
                tickfont=dict(color='#B03060'),
                showgrid=False, showline=True, zeroline=False,
                showticklabels=True),}
        return trace


    def plot_occupancy(self, occupancy_trace, rmsd_trace, distance_trace):
        fig1 = go.Figure(occupancy_trace)
        fig1.update_layout(**occupancy_trace.meta['layout'])
        fig1.update_xaxes(**occupancy_trace.meta['xaxis_params'])
        fig1.update_yaxes(**occupancy_trace.meta['yaxis_params'])

        fig2 = go.Figure()
        fig2.add_trace(rmsd_trace)
        fig2.add_trace(distance_trace)
        fig2.update_layout(**rmsd_trace.meta['layout'])
        fig2.update_xaxes(**rmsd_trace.meta['xaxis_params'])

        # Compute tickvals and ranges
        rmsd_upper, rmsd_tickvals, dist_upper, dist_tickvals = utils.nice_ticks(
            np.nanmax(rmsd_trace.y), np.nanmax(distance_trace.y))
        rmsd_upper += rmsd_upper * 0.1
        dist_upper += dist_upper * 0.1

        fig2.update_layout(
            yaxis=dict(**rmsd_trace.meta['yaxis_params'],
                tickmode='array', tickvals=rmsd_tickvals,
                range=[0, rmsd_upper]),
            yaxis2=dict(**distance_trace.meta['yaxis_params'],
                tickmode='array', tickvals=dist_tickvals,
                range=[0, dist_upper], overlaying='y',
                side='right'))

        # Attach distance trace to right axis
        fig2.data[1].update(yaxis=distance_trace.meta['axis_id'])

        html = (pio.to_html(fig1, include_plotlyjs='cdn', full_html=False, div_id='occ') +
            pio.to_html(fig2, include_plotlyjs=False, full_html=False, div_id='lines'))

        # Generate PNG layout
        png_fig = make_subplots(
            rows=2, cols=1,
            shared_xaxes=True,
            vertical_spacing=0.04,
            row_heights=[0.67, 0.33],
            specs=[[{}], [{'secondary_y': True}]])
        png_fig.add_trace(go.Heatmap(occupancy_trace.to_plotly_json()), row=1, col=1)
        png_fig.add_trace(go.Scatter(rmsd_trace.to_plotly_json()), row=2, col=1, secondary_y=False)
        png_fig.add_trace(go.Scatter(distance_trace.to_plotly_json()), row=2, col=1, secondary_y=True)

        png_fig.update_layout(
            height=650, width=1100,
            margin=dict(l=100, r=70, t=65, b=65),
            title=occupancy_trace.meta['layout']['title'],
            plot_bgcolor=occupancy_trace.meta['layout']['plot_bgcolor'],
            paper_bgcolor='#FFFFFF',
            showlegend=False,
            hovermode='x unified',
            hoverlabel=rmsd_trace.meta['layout']['hoverlabel'])
        png_fig.update_xaxes(row=1, col=1, **occupancy_trace.meta['xaxis_params'])
        png_fig.update_yaxes(row=1, col=1, **occupancy_trace.meta['yaxis_params'])
        png_fig.update_xaxes(row=2, col=1,
            title_text=rmsd_trace.meta['layout']['xaxis_title'],
            **rmsd_trace.meta['xaxis_params'])
        png_fig.update_yaxes(row=2, col=1, secondary_y=False,
            tickmode='array', tickvals=rmsd_tickvals, range=[0, rmsd_upper],
            **rmsd_trace.meta['yaxis_params'])
        png_fig.update_yaxes(row=2, col=1, secondary_y=True,
            tickmode='array', tickvals=dist_tickvals, range=[0, dist_upper],
            **distance_trace.meta['yaxis_params'])

        # Restore white background for subplot 2
        row2_xdomain = png_fig.layout.xaxis2.domain
        row2_ydomain = png_fig.layout.yaxis2.domain
        png_fig.add_shape(
            type='rect',
            xref='paper', yref='paper',
            x0=row2_xdomain[0], x1=row2_xdomain[1],
            y0=row2_ydomain[0], y1=row2_ydomain[1],
            fillcolor='#FFFFFF',
            line=dict(width=0),
            layer='below')

        # Fit heatmap colorbar to the occupancy subplot
        occ_domain = png_fig.layout.yaxis.domain
        occ_center = (occ_domain[0] + occ_domain[1]) / 2
        occ_height = (occ_domain[1] - occ_domain[0]) * 1.055
        occ_xmax = png_fig.layout.xaxis.domain[1]
        png_fig.update_traces(selector=dict(type='heatmap'),
            colorbar=dict(x=occ_xmax, y=occ_center, len=occ_height,
            yanchor='middle', thickness=18))

        return html, png_fig


    def run_analysis(self):
        int_csv = self.frame_dir / 'interactions.csv'
        if not self.frame_dir.exists():
            self.frame_dir.mkdir(parents=True, exist_ok=True)
            self.extract_frames()
        if not int_csv.exists():
            pdb_files = glob.glob(f"{self.frame_dir}/*.pdb")
            df_int = inttools.analyse_pdb_files(pdb_files)
            df_int.to_csv(int_csv, index=False)
        else:
            print(f"Using file {int_csv}...")
            df_int = pd.read_csv(int_csv)
        return df_int
        

    def plot_analysis(self, df_pivot):
        occupancy_trace = self.make_occupancy_trace(df_pivot)
        rmsd_trace = self.make_rmsd_trace()
        distance_trace = self.make_distance_trace()
        html, png = self.plot_occupancy(occupancy_trace, rmsd_trace, distance_trace)
        return html, png


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Protein-ligand time-series analysis")
    parser.add_argument('--trj_name', required=True,
        help="Trajectory file name for time series analysis")
    args = parser.parse_args()

    configs = utils.read_config_file('config.yaml')
    workdir = Path(configs['workdir'])

    occ = OccupancyAnalysis(args.trj_name, configs)
    df_int = occ.run_analysis()

    for dct in configs['ligands']:
        lig_name = dct.get('name', 'ligand')
        lig_resname = dct.get('md_id', dct.get('resname', 'LIG'))
        int_dir = workdir / f'{lig_name}.interactions'
        int_dir.mkdir(parents=True, exist_ok=True)

        df_lig = df_int[df_int['RESNAME'] == lig_resname].copy()
        aa_keys = set(utils.standard_AAs.keys())
        resid_name = df_lig['RESID'].astype(str).str.extract(R'^([A-Za-z]+)', expand=False)
        df_lig = df_lig[resid_name.isin(aa_keys)].reset_index(drop=True)
        if df_lig.empty:
            continue
        df_lig.to_csv(int_dir / 'interactions.csv', index=False)
        df_pivot = occ.pivot_df(df_lig)

        occ.pdb_id = configs['receptor_id']
        occ.ligand_resname = lig_resname
        html, png = occ.plot_analysis(df_pivot)

        html_path = workdir / f"trj_occupancy_{lig_resname}.html"
        with html_path.open('w') as f:
            f.write(html)
        print(f"Saved {html_path.name}")

        png_path = workdir / f"trj_occupancy_{lig_resname}.png"
        png.write_image(png_path, scale=3)
        print(f"Saved {png_path.name}")
