###########################
#    Iva Laginja, 2016    #
# last change, 16-02-2017 #
###########################

"""
Creating the data for an image grid of ring systems and light curves, with different i and phi.
Created data are a fits, dat and txt file. They will be saved to a subfolder specifed by "targetdir" and can be further used with the ring drawing code.
This code uses 'bring_disk_sim_data.py' and the here created output data is used by 'grid_lightcurve_plots.py'.
"""
import os
import exomoons_py.main.bring_disk_sim_data as disksim
from exomoons_py.config import CONFIG_INI


# constants (taken from disk_sim.py)
G = 6.6738480e-11  # m3 kg-1 s-2
yr = 365.242189669 * 86400  # sec
msol = 1.98855e30  # kg
rsol = 6.5500e8  # m
mjup = 1.8986e27  # kg
rjup = 6.9911e7  # m
mearth = 5.97219e24  # kg
rearth = 6.371e3  # kg (Wikipedia)
mmoon = 7.3476e22  # kg
au = 1.49597870700e11  # m
pc = 3.0856e16  # m

# --------------------------------------------#
local_path = CONFIG_INI.get('data_paths', 'local_data_path')
current_exp = CONFIG_INI.get('data_paths', 'curr_data_path')

targetdir = os.path.join(local_path, current_exp)

impact = CONFIG_INI.getfloat('beta_Pic', 'impact')  # impact parameter in terms of parts of the Hill radius => b = impact * R_Hill

# Template: bring_disk_data(i_in, phi_in, output, targetdir, modelnumber)

# Create data

print("\nMAKING MODEL 0 NOW:")
disksim.bring_disk_data(impact, 30., 0., "model0", targetdir, '0', paramfile=True)

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
