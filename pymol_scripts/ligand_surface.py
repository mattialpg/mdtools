# -*- coding: utf-8 -*-
from pymol import cmd

def ligand_surface(selection='org', type='gaussian_max'):
    
    cmd.set('surface_carve_cutoff', '4.5')
    cmd.set('surface_carve_selection', selection)
    cmd.set('surface_carve_normal_cutoff', -0.1)
    cmd.set('surface_quality', 2)

    cmd.set('ray_transparency_contrast', 0.20)
    cmd.set('ray_transparency_oblique', 1.0)
    cmd.set('ray_transparency_oblique_power', 20)

    cmd.set('two_sided_lighting')
    cmd.set('transparency', 0.5)
    cmd.set('ray_shadows', 0)

    cmd.set('gaussian_resolution', 2)
    cmd.ramp_new('pRamp', selection, selection=selection, range=[1.5,2.5], color='rainbow')


    if type=='gaussian_max':
        ## Gaussian surface approximating vdW surface ##
        cmd.map_new('gauss_map', 'gaussian_max', selection=selection, grid=.3, buffer=5)
        cmd.isosurface('gauss_surf', 'gauss_map', level=0.2)
        cmd.set('surface_color', 'pRamp', 'gauss_surf')
        cmd.delete('gauss_map')
    else:
        ## True vdW surface ##
        cmd.map_new('vdw_map', 'vdw', grid=.2, selection=selection, buffer=3)
        cmd.isosurface('vdw_surf', 'vdw_map', level=0.01, buffer=3)
        cmd.set('surface_color', 'pRamp', 'vdw_surf')
        cmd.delete('vdw_map')

    # color protein residues by distance from ligand
    #> show within 4
    #color pRamp, prot
    
    # cmd.show('mesh', 'prot within 8 of lig')
#    cmd.set('surface_type', 2)
#    cmd.set('surface_color', 'white')
cmd.extend ("ligand_surface", ligand_surface)
