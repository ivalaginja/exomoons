# exomoons

These are Python tools for creating and displaying exorings and exomoon artefacts (rings) in circumplanetary meterial,
with the main goal to generate light curves.
This code is based on and adapted from Matthew Kenworthy's GitHub repository "exorings", which in turn is based on
the paper by Kenworthy and Mamajek (2015):
https://github.com/mkenworthy/exorings

Full licence can be found in "LICENCE.txt".

Currently I am not actively working on this code but offering support to people who (want to) use it.
There is still a lot of instances in the code where things are hard coded and need to be moved to the configfile (as of 10/08/2018).

The code is written in Python 3 but is Python 2 compatible. I would encourage anyone to use Python 3.x though, since Python 2 is just not a thing anymore. Really. Just write Python 3 code.


Current instructions to make stuff run:


## Quick run for direct example results from template:

- cd into the directory you want to have your repo in
- $ git clone git@github.com:ivalaginja/exomoons.git
- $ cd -repo-/exomoons_py
- open config.ini with a text editor and change the entry for "local_data_path" to the directory you want to save the data to and "cur_data_path" to the name you want the experiment folder to have
- $ python data_creation.py
- $ python grid_lightcurve_plots.py
- look at the result in the directory you defined for local_data_path/cur_data_path !


## Going beyond just running the default code:

--- FIRST: SET UP THE CONFIGURATION FILE

The main configuration file is config.ini, which holds all of your simulation paramers. This file,
however, is version controlled, and the paths to local directories will get messed up if you push this
file. This is why config.ini is supposed to be a TEMPLATE. In order to make it work for you,
copy config.ini and rename the copy to config_local.ini. This new file will override anything present
in config.ini. Make sure you tell your version control system to ignore config_local.ini.


--- CREATING DATA (you can skip steps 1-3 if data already exists) ---

1) Open config_local.ini and set your system path and system parameters.
    - "local_data_path" sets your overal directory you want to have data from these simulations saved to
    - "curr_data_path" sets the current specific subdirectory you want an experiment saved to (need to change this for each new run/experiment)

2) Set the physical parameters of your exomoon system in the different sections of the config_local.ini

3) Run data_creation.py


--- CREATING PLOT ---

4) In grid_lightcurve_plots.py

    a) can change figure suptitle
    
    b) maybe need to adjust days_side and y_range

5) run grid_lightcurve_plots.py
    => it uses the data in folder "data" to make a grid in the same folder


y-axis size of the ring plots can be set in draw_rings.py, #new_y_rings; but technically I set them in grid_lightcurve_plots.py (parameter "y_range")


## Contributing

If you want to suggest changes to this content do the following:

1. Fork it
2. Create your feature branch (`git checkout -b my-new-feature`)
3. Commit your changes (`git commit -am 'Added some feature'`)
4. Push to the branch (`git push origin my-new-feature`)
5. Create new Pull Request
