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

1) Open config.ini or config_local.ini and set your system pathd and system parameters.

2) Run data_creation.py

--- CREATING PLOT ---

3) In grid_lightcurve_plots.py
    a) can change figure suptitle
    b) maybe need to adjust days_side and y_range

4) run grid_lightcurve_plots.py
    => it uses the data in folder "data" to make a grid in the same folder


y-axis size of the ring plots can be set in draw_rings.py, #new_y_rings; but technically I set them in grid_lightcurve_plots.py (parameter "y_range")