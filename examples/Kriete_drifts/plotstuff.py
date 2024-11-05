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

Tnorm = 80

qe = 1.602176634e-19           # Electron charge
Mp = 1.672621898e-27   # Proton mass
ion_sound_speed = np.sqrt(2*Tnorm * qe / Mp)


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
field_line_pitch_angle = 0.001 
drift_velocity = -62
poloidal_velocity = parallel_velocity * np.sin(field_line_pitch_angle) + drift_velocity * np.cos(field_line_pitch_angle)
max_density = max(density.values)
length = 640
y_xpoint = 80
y = ysim = ds.y


fig, ax1 = plt.subplots(constrained_layout=True,figsize=(7,5))
ax2 = ax1.twinx()
ax1.set_xticks([0, y_xpoint-1, length/2, length-y_xpoint, length])
ax1.set_xticklabels(['lower divertor', 'lower X-point', '0', 'upper X-point', 'upper divertor'], minor=False, rotation=60)
ax1.plot(y, parallel_velocity, linestyle='--', color='tab:blue', label=r'$v_\parallel$')
scale_factor = 1 / np.sin(field_line_pitch_angle)
ax1.plot(y, scale_factor * poloidal_velocity, linestyle=':', color='tab:blue', label=fr'${scale_factor:.0f} * v_\theta$')
ax1.set_yticks([-ion_sound_speed, 0, ion_sound_speed])
ax1.set_yticklabels(['$-c_s$', '0', '$c_s$'], minor=False)
ax1.tick_params(axis='y', labelcolor='tab:blue')
ax1.set_ylabel(r'Velocity', color='tab:blue')
ax1.legend()
ax1.grid()

ax2.plot(y, density, color='tab:red')
ax2.set_ylim(bottom=0)
ax2.set_yticks([0, max_density])
ax2.set_yticklabels(['0', '$n_\mathrm{{max}}$'], minor=False)
ax2.tick_params(axis='y', labelcolor='tab:red')
ax2.set_ylabel('Density', color='tab:red')
plt.show()

min_density = min(density.values)
print("the max density is: ", max_density)
print("the upstream density is: ", density.values[0])
print("the min density is: ", min_density)
gamma = drift_velocity/np.tan(field_line_pitch_angle)/2/ion_sound_speed
print("The upstream density over the downstream density must be (1+\gamma)/(1-\gamma), where \gamma is vD/(2Cs), so this should be 1: ",
      (density.values[-1]/density.values[0])/((1+gamma)/(1-gamma)))
