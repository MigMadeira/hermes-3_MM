#!/usr/bin/env python3

import xhermes

# Following two lines needed so all variables are shown when printing the Dataset
import xarray as xr
xr.set_options(display_max_rows=1000)

# Set better figure size
from matplotlib import pyplot as plt
plt.rcParams["figure.figsize"] = (16,8)

ds = xhermes.open(".")
#print(ds)

# Get rid of size-1 radial and toroidal directions
ds = ds.squeeze()

import numpy as np

Tnorm = 50

qe = 1.602176634e-19           # Electron charge
Mp = 1.672621898e-27   # Proton mass
ion_sound_speed = np.sqrt(2*Tnorm * qe / Mp)


# In case you want to check the source.
# plt.figure()
# plt.plot(ds.SNi[-1,:])
# plt.plot(ds.SNi[0,:])
# plt.show()

# Make an animation
# Note: saving a gif can be slow. Comment out `save_as` argument to disable.
ds.bout.animate_list(
    ["Ni", "Vi"],
    show=True,
    save_as="hermes_animation",
)

parallel_velocity = ds.Vi[-1,:] #cutting first and last point as they have some weird values.
density = ds.Ni[-1,:]


#these parameters need to be set according to the input file. 
max_density = max(density.values)
length = 640
y_xpoint = 80

fig, ax1 = plt.subplots(constrained_layout=True, figsize=(7,5))
ax2 = ax1.twinx()
ax1.set_xticks([0, y_xpoint, length/2, length-y_xpoint, length])
ax1.set_xticklabels(['lower divertor', 'lower X-point', '0', 'upper X-point', 'upper divertor'], minor=False, rotation=45)
ax1.set_title('Profiles in absence of drifts')

ax1.plot(parallel_velocity, color='tab:blue')
# ax1.axhline(0, linestyle='--')
ax1.set_yticks([-ion_sound_speed, 0, ion_sound_speed])
ax1.set_yticklabels(['$-c_s$', '0', '$c_s$'], minor=False)
ax1.tick_params(axis='y', labelcolor='tab:blue')
ax1.set_ylabel(r'Parallel velocity', color='tab:blue')

ax2.plot(density, color='tab:red')
ax2.set_ylim(bottom=0)
ax2.set_yticks([0, max_density / 2, max_density])
ax2.set_yticklabels(['0', '$n_d$', '$2 n_d$'], minor=False)
ax2.tick_params(axis='y', labelcolor='tab:red')
ax2.set_ylabel('Density', color='tab:red')

density_kriete = np.load("../plot_db/Kriete_density.npy")
parallel_velocity_kriete = np.load("../plot_db/Kriete_parallel_velocity.npy")#*np.sqrt(2)
y_kriete = np.load("../plot_db/Kriete_y.npy")
ax1.plot((y_kriete+1)*320, parallel_velocity_kriete, linestyle='--', color='k', label=r'$v_\parallel$_Kriete', alpha = 0.5)
ax2.plot((y_kriete+1)*320, density_kriete, color='k', alpha = 0.5, label = "density_Kriete")

ax1.grid("both")
ax2.grid()
plt.show()

print("nd = ",max_density/2)
print("Cs divided by the maximum parallel velocity (should be 1): ",ion_sound_speed/max(parallel_velocity.values))
print("max density over min density (should give 2.0): ",max_density/min(density.values))