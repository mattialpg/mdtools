


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

        # Save matrix to CSV
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