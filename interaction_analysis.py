import warnings
warnings.filterwarnings('ignore')

import os, sys, glob
import re, yaml
from pathlib import Path
import asyncio, argparse
import subprocess
import numpy as np
import pandas as pd
from tqdm import tqdm
from itertools import combinations
from collections import Counter
import mdtraj as md

import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
from plotly.subplots import make_subplots

import utils, palettes

# Import AIDDTools utils
sys.path.append('/home/mattia/aiddtools')
import fetching_tools as fetchtools
import interaction_tools as inttools


class TimeSeriesAnalysis:
    def __init__(self, trj_name, config_path='config.yaml'):
        with open(config_path) as f:
            self.configs = yaml.safe_load(f)
        self.ligand_name = self.configs['ligand_name']
        self.ligand_md_id = self.configs['ligand_md_id']
        self.pdb_id = self.configs['receptor_id']
        self.trj_name = trj_name

        self.workdir = Path(self.configs['workdir'])
        self.int_dir = self.workdir / f'{self.ligand_name}.interactions'
  
        self.int_types = ['hydrophobic', 'hbond', 'waterbridge',
             'saltbridge', 'pistacking', 'pication', 'halogen', 'metal']

        self._ligand_atoms = None
        self._ligand_masses = None

        print("Loading trajectory...")
        self.traj = md.load_xtc(f'{self.trj_name}.xtc', top=f'{self.trj_name}.gro')


    def _get_ligand_atoms(self):
        if self._ligand_atoms is None:
            self._ligand_atoms = self.traj.topology.select(f'resname {self.ligand_md_id}')
            if not len(self._ligand_atoms):
                raise ValueError(f"Ligand '{self.ligand_md_id}' not found")
        return self._ligand_atoms


    def _get_ligand_masses(self, ligand_atoms):
        if self._ligand_masses is None:
            masses = np.array([atom.element.mass for atom in\
                self.traj.topology.atoms if atom.index in ligand_atoms])
            masses /= masses.sum()
            self._ligand_masses = masses
        return self._ligand_masses


    def extract_frames(self):
        total_frames = self.traj.n_frames
        n_extract = 1001
        frame_indices = np.linspace(0, total_frames - 1, n_extract, dtype=int)\
            if total_frames >= n_extract else np.arange(total_frames)
        for i, idx in enumerate(tqdm(frame_indices, desc="Exporting MD frames", unit='frame')):
            self.traj[idx].save_pdb(f'{self.int_dir}/frame_{i+1:04d}.pdb')

        # subprocess.run(['gmx', 'trjconv', '-f', f'{self.trj_name}.xtc',
        #     '-s', f'{self.trj_name}.gro', '-o', f'{self.int_dir}/frame_.pdb',
        #     '-sep', '-nzero', '4'], input="System\n", text=True, check=True)


    def make_occupancy_trace(self, df_int):
        dcolorsc, colorbar_cfg = palettes.discrete_colorscale(
            'geodataviz_diverging_1', self.int_types)
        colorbar_cfg = dict(colorbar_cfg, x=1.0, y=0.5,
            xanchor='left', yanchor='middle', len=1.048)

        int_priority = {name: idx for idx, name in enumerate(self.int_types)}
        df_lig = df_int[df_int['LIGNAME'] == self.ligand_md_id].copy()
        df_lig['INT_PRIORITY'] = df_lig['INT_TYPE'].map(int_priority)
        df_lig = df_lig.sort_values('INT_PRIORITY')

        df_pivot = df_lig.pivot_table(index='PDB', columns='RESID', values='INT_TYPE', aggfunc='first')
        df_pivot = df_pivot.drop(columns=['NANA'], errors='ignore')
        df_pivot = df_pivot.replace(int_priority)

        occupancy = df_pivot.notna().sum() / df_pivot.shape[0]
        kept = occupancy[occupancy >= 0.05].index
        df_pivot = df_pivot[kept]

        corr = df_pivot.fillna(0).corr(method='pearson')
        ordered = corr.mean().sort_values(ascending=False).index
        df_pivot = df_pivot[ordered]

        # Build time vector scaled to number of frames
        total_time_ns = self.traj.time[-1] / 1000
        n_frames = len(df_pivot.index)
        time_ns = np.linspace(0, total_time_ns, n_frames)

        trace = go.Heatmap(
            x=time_ns, 
            y=df_pivot.columns,
            z=df_pivot.values.T,
            colorscale=dcolorsc,
            zmin=0, zmax=len(self.int_types),
            colorbar=colorbar_cfg, name='Occupancy',
            customdata=np.vectorize(lambda i: self.int_types[int(i)]
                if i == i else 'NA')(df_pivot.values.T),
            hovertemplate="Time=%{x:.2f} ns<br>Residue=%{y}<br>Interaction=%{customdata}<extra></extra>")
        
        trace.meta = {
            'layout': dict(
                height=500, width=1100,
                margin=dict(l=100, r=50, t=65, b=0),
                title=dict(
                    text=f"Interaction occupancy ({self.pdb_id}·{self.ligand_md_id})",
                    x=0.5, y=0.95, xanchor='center', yanchor='top',
                    font=dict(size=20, weight='bold')),
                plot_bgcolor='#C9D0D6',
                hovermode='closest'),
            'xaxis_params': dict(showgrid=False, showticklabels=False),
            'yaxis_params': dict(
                title_text='<b>Residue</b>', tickmode='linear',
                showgrid=False),}
        return trace


    def make_rmsd_trace(self):
        ligand_atoms = self._get_ligand_atoms()
        ligand_traj = self.traj.atom_slice(ligand_atoms)
        rmsd = md.rmsd(ligand_traj, ligand_traj, 0)

        trace = go.Scatter(
            x=self.traj.time / 1000,
            y=rmsd,
            mode='lines',
            line=dict(color='#1A2228', width=1.2),
            name='RMSD (nm)',
            hovertemplate='RMSD=%{y:.3f} nm<extra></extra>')

        trace.meta = {
            'layout': dict(
                height=250, width=1040,
                margin=dict(l=100, r=50, t=10, b=65),
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
                unifiedhovertitle=dict(text='Time=%{x:.2f} ns')),
            'yaxis_params': dict(
                title=dict(text='<b>RMSD (nm)</b>', font=dict(size=12)),
                gridcolor='lightgrey',
                showgrid=True,
                zeroline=True, zerolinecolor='lightgrey',
                showticklabels=True)}
        return trace


    def make_distance_trace(self):
        ligand_atoms = self._get_ligand_atoms()
        masses = self._get_ligand_masses(ligand_atoms)
        com = np.average(self.traj.xyz[:, ligand_atoms, :], axis=1, weights=masses)
        com_ref = com[0]
        dist = np.linalg.norm(com - com_ref, axis=1)

        trace = go.Scatter(
            x=self.traj.time / 1000,
            y=dist,
            mode='lines',
            line=dict(color='#B03060', width=1.2),
            name='Distance (Å)',
            hovertemplate='Distance=%{y:.2f} Å<extra></extra>')
        
        trace.meta = {
            'axis_id': 'y2',
            'yaxis_params': dict(
                title=dict(text='<b>Distance (Å)</b>',
                    font=dict(color='#B03060', size=12)),
                # tickvals=tickvals, ticktext=ticktext, range=[0, upper],
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
        with open(f'{self.int_dir}/../trj_occupancy.html', 'w') as f:
            f.write(html)
        print("Saved trj_occupancy.html")


    def run(self):
        if not os.path.exists(self.int_dir):
            os.makedirs(self.int_dir)
            self.extract_frames()
            pdb_files = glob.glob(f'{self.int_dir}/*.pdb')
            df_int = inttools.analyse_pdb_files(pdb_files)
            df_int.to_csv(f'{self.int_dir}/interactions.csv', index=False)
        else:
            print(f"Using file {self.int_dir}/interactions.csv...")
            df_int = pd.read_csv(f'{self.int_dir}/interactions.csv')

        occupancy_trace = self.make_occupancy_trace(df_int)
        rmsd_trace = self.make_rmsd_trace()
        distance_trace = self.make_distance_trace()
        self.plot_occupancy(occupancy_trace, rmsd_trace, distance_trace)


class CrossTargetAnalysis:
    def __init__(self, pdb_files):
        self.pdb_files = pdb_files
        self.pair_counts = Counter()
        self.matrix = None

    def analyse_receptors(self):
        # Analise multiple receptors
        df_int = inttools.analyse_pdb_files(self.pdb_files,
            xml_outdir='.', pdb_outdir='.')

        # Identify ligands of interest (LOI)
        lig_cci = set()
        for file in self.pdb_files:
            ligs = inttools.get_ligands(file, select_loi=True)
            lig_cci.update([x.split(':')[0] for x in ligs.keys()])
        df_int = df_int[df_int['LIGNAME'].isin(lig_cci)]
        df_int.to_csv('df_int.csv', index=False)

        # Count co-occurrences across receptors
        res_per_pdb = df_int.groupby('PDB')['RESID'].apply(set)
        self.pair_counts = Counter()
        for residues in res_per_pdb:
            for r1, r2 in combinations(sorted(residues), 2):
                self.pair_counts[(r1, r2)] += 1

        # Build co-occurrence matrix
        all_residues = sorted({r for pair in self.pair_counts for r in pair})
        self.matrix = pd.DataFrame(0, index=all_residues, columns=all_residues)
        for (r1, r2), count in self.pair_counts.items():
            self.matrix.loc[r1, r2] = count
            self.matrix.loc[r2, r1] = count

        # Reorder matrix using weighted order
        values = sorted(self.matrix.stack().unique(), reverse=True)
        counts = pd.DataFrame({v: (self.matrix == v).sum(axis=1) for v in values})
        counts['sort_key'] = list(counts.itertuples(index=False, name=None))
        sorted_residues = counts.sort_values(by=values, ascending=False).index.tolist()
        self.matrix = self.matrix.loc[sorted_residues, sorted_residues]

        # svae matrix to CSV
        self.matrix.to_csv('matrix.csv')


    def plot_heatmap(self):

        colorscale = [[float(p), c] for p, c in zip(
            np.arange(0, 1.01, 0.25), palettes.palette_5)]
        
        fig = px.imshow(
            self.matrix,
            x=self.matrix.columns,
            y=self.matrix.index,
            color_continuous_scale=colorscale,
            zmin=0,
            zmax=self.matrix.to_numpy().max(),
            aspect='auto')

        fig.update_layout(
            showlegend=False,
            coloraxis_showscale=True,
            xaxis=dict(
                title='Residue',
                tickangle=-90,
                title_font=dict(size=16, weight='bold'),
                showgrid=False, zeroline=False,
                linecolor='black', mirror=True,
                scaleanchor='y', constrain='domain'),
            yaxis=dict(
                title='Residue',
                title_font=dict(size=16, weight='bold'),
                showgrid=False, zeroline=False,
                linecolor='black', mirror=True,),
            # title=dict(
            #     text="Antibiogram",
                # x=0.5,  y=0.97, xanchor='center', yanchor='top',
                # font=dict(size=22, weight='bold')),
            coloraxis_colorbar=dict(
                tickmode='array',
                tickvals=list(range(int(self.matrix.to_numpy().max()) + 1)),
                ticktext=[str(i) for i in range(int(self.matrix.to_numpy().max()) + 1)]),
            height=650, width=750,
            margin=dict(l=60, r=50, t=50, b=105))

        fig.write_image('matrix.png', scale=5)
        print("Saved matrix.png")


    def find_top_residue_triplets(self):
        # Compute residue totals
        # residue_totals = self.matrix.sum(axis=1).astype(int).to_dict()

        # Select all pairs with maximum co-occurrence frequency
        max_score = max(self.pair_counts.values())
        top_pairs = [pair for pair, score in self.pair_counts.items() if score == max_score]

        # Combine those top pairs into unique triplets
        triplets = set()
        for p1, p2 in combinations(top_pairs, 2):
            combo = set(p1 + p2)
            if len(combo) == 3:  # valid triplet if exactly three unique residues
                combo = sorted(combo, key=lambda x: int(x[3:]))  # sort numerically by residue number
                # score = sum(residue_totals.get(x, 0) for x in combo)
                triplets.add(tuple(combo))
        
        print("Top residue triplets:")
        for triplet in triplets:
            print('   ' + '-'.join(triplet))
        return triplets
    

    def run(self):
        self.analyse_receptors()
        self.plot_heatmap()
        self.find_top_residue_triplets()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Protein-ligand analysis")
    parser.add_argument('--analysis', choices=['time', 'targets', 'both'], default='time',
        help="Select which analysis to run: time = time series only, targets = cross-target only, both = run both")
    parser.add_argument('--trj_name', help="Trajectory file name for time series analysis")
    parser.add_argument('--pdb_files', help="Path or pattern to PDB files for cross-target analysis (e.g. \"./*.pdb\")")
    parser.add_argument('--cif_files', help="Path or pattern to CIF files for cross-target analysis (e.g. \"./*.cif\")")
    args = parser.parse_args()

    if args.analysis in ['time', 'both']:
        if not args.trj_name:
            sys.exit("Provide --trj_name for time series analysis")
        time_analysis = TimeSeriesAnalysis(trj_name=args.trj_name)
        time_analysis.run()

    if args.analysis in ['targets', 'both']:
        pdb_files = glob.glob(args.pdb_files if args.pdb_files else '*.pdb')
        if not pdb_files:
            try:
                cif_files = glob.glob(args.cif_files if args.cif_files else '*.cif')
                for cif in cif_files:
                    asyncio.run(fetchtools.download_pdbs_async(cif.split('.')[0]))
            except:
                sys.exit("No PDB files found for cross-target analysis")
                
        target_analysis = CrossTargetAnalysis(pdb_files)
        target_analysis.run()
