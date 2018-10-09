###########################
#    Iva Laginja, 2016    #
###########################

"""
Creating the data for an image grid of ring systems and light curves, with different i and phi.
Created data are a fits, dat and txt file. They will be saved to a subfolder specifed by "targetdir" and can be further used with the ring drawing code.
This code uses 'bring_disk_sim_data.py' and the here created output data is used by 'grid_lightcurve_plots.py'.

Takes currently about 3 minutes to run.
"""

import sys
import os
from shutil import copy
import time
import exomoons_py.main.bring_disk_sim_data as disksim
from exomoons_py.config import CONFIG_INI


if __name__ == '__main__':

    local_path = CONFIG_INI.get('data_paths', 'local_data_path')
    current_exp = CONFIG_INI.get('data_paths', 'curr_data_path')

    # Define directory to work in and create it
    targetdir = os.path.join(local_path, current_exp)
    if os.path.isdir(targetdir):
        sys.exit('ERROR: The directory ' + targetdir + ' already exists, please delete or change folder to save data to in config.ini.')
    else:
        os.mkdir(targetdir)

    # Keeping track of code runtime
    start_time = time.time()

    impact = CONFIG_INI.getfloat('beta_Pic', 'impact')  # impact parameter in terms of parts of the Hill radius => b = impact * R_Hill

    # Copy the configfile to the experiment folder.
    #copy(os.path.join('..', 'config_local.ini'), targetdir)
    copy(os.path.join('..', 'config.ini'), targetdir)

    # Template: bring_disk_data(i_in, phi_in, output, targetdir, modelnumber)

    # Create data

    start_time_single = time.time()   # Keeping track of code runtime for single model
    print("\nMAKING MODEL 0 NOW:")
    disksim.bring_disk_data(impact, 30., 0., "model0", targetdir, '0', paramfile=True)
    end_time_single = time.time()


    print("\nMAKING MODEL 1 NOW:")
    disksim.bring_disk_data(impact, 30., 30., "model1", targetdir, '1')

    print("\nMAKING MODEL 2 NOW:")
    disksim.bring_disk_data(impact, 30., 60., "model2", targetdir, '2')

    print("\nMAKING MODEL 3 NOW:")
    disksim.bring_disk_data(impact, 45., 0., "model3", targetdir, '3')

    print("\nMAKING MODEL 4 NOW:")
    disksim.bring_disk_data(impact, 45., 30., "model4", targetdir, '4')

    print("\nMAKING MODEL 5 NOW:")
    disksim.bring_disk_data(impact, 45., 60., "model5", targetdir, '5')

    print("\nMAKING MODEL 6 NOW:")
    disksim.bring_disk_data(impact, 60., 0., "model6", targetdir, '6')

    print("\nMAKING MODEL 7 NOW:")
    disksim.bring_disk_data(impact, 60., 30., "model7", targetdir, '7')

    print("\nMAKING MODEL 8 NOW:")
    disksim.bring_disk_data(impact, 60., 60., "model8", targetdir, '8')

    # append stuff to parameter file
    params = open(os.path.join(targetdir, "system_parameters.txt"), "a")
    params.write("Folder name: %s \n" % targetdir)
    params.write("Impact parameter: %s of R_Hill\n" % impact)
    params.close()

    print('')
    print('---------------------------')
    print('impact = ', impact)
    print('DATA SAVED TO ', targetdir)
    print('---------------------------')

    #### Ending program
    end_time = time.time()

    print('Runtime:', str((end_time-start_time)/60) + " min")
    print('(Runtime of one model:', str((end_time_single-start_time_single)/60) + ' min)')
