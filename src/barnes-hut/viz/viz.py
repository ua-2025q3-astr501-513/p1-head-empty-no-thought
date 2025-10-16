import glob
import numpy as np
import pyvista as pv

import os
import shutil


def make_viz( datadir: str, rad_multiplier = 2) :

    OUTPATH = os.path.join('frames', datadir)
    DATPATH = os.path.join('..', 'data', datadir, f'globr_{datadir}_*.dat')

    shutil.rmtree( OUTPATH, ignore_errors=True) # removes directory if exists already before rendering
    os.mkdir(  OUTPATH  )
    directory = f'{datadir}'
    files = sorted(glob.glob(DATPATH))

    pv.OFF_SCREEN = True 

    plotter = pv.Plotter(off_screen=True)
    plotter.set_background("black")

    for i, fname in enumerate(files):
        plotter.clear()

        data = np.loadtxt(fname, comments='#', skiprows=11)
        mass = data[:,1]
        xyz = data[:,2:5]

        # create point cloud object
        cloud = pv.PolyData(xyz)

        # scale point size by normalized mass
        sizes = (mass / mass.max()) * 5  # tweak value

        plotter.add_points(cloud, render_points_as_spheres=True, point_size=2, color="white")

        # fixed camera
        # plotter.camera_position = "xz"  # you can also set ((x,y,z), focus, up)

        # Let PyVista auto-center the view\n,
        center = xyz.mean(axis=0)  # dynamic center of mass
        radius = np.linalg.norm(xyz - center, axis=1).max()  # size of system
        
        plotter.camera_position = [
            center + [0, 0, radius * rad_multiplier],  # camera location (adjust Z multiplier for zoom)
            center,                       # focal point
            [0, 1, 0]                     # view up direction
        ]

        plotter.screenshot(f"frames/{directory}/frame_{i:05d}.png")
