# exomoons

These are Python tools for creating and displaying exorings and exomoon artefacts (rings) in circumplanetary meterial,
with the main goal to generate ligh curves.
This code is based on and adapted from Matthew Kenworthy's GitHub repository "exorings", which in turn is based on
the paper by Kenworthy and Mamajek (2015).

Full licence can be found in "LICENCE.txt".

Currently I am bringing the code up to speed with Python 3 and trying to restructure it so that it is easier to use
and understand.

#######################################
Current instructions to make stuff run:

--- CREATING DATA (you can skip steps 1-4 if data already exists) ---

1) create the data folder in which you want to work (e.g. test)

2) open bring_disk_sim_data.py (You never actually run this code. It gets used by other codes.)
    a) set masses m
    b) set moons' semi-major axes a
    c) adjust any other parameters, like ring radii or tau values

3) open data_creation.py
    a) set tergetdir to folder to which data should be saved (e.g. test)
    b) set impact - it's the impact parameter in terms of parts of the Hill radius  =>  b = impact * R_Hill

4) run data_creation.py

--- CREATING PLOT ---

5) open grid_lightcurve_plots.py
    a) set target folder "data" (to e.g. test)
    b) set figure suptitle
    c) maybe need to adjust days_side and y_range

6) run grid_lightcurve_plots.py
    => it uses the data in folder "data" to make a grid in the same folder


y-axis size of the ring plots can be set in draw_rings.py, #new_y_rings; but technically I set them in grid_lightcurve_plots.py (parameter "y_range")