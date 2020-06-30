import espressomd
from espressomd import interactions
from espressomd import polymer
from espressomd import electrostatics
from espressomd.io.writer import vtf


import numpy as np
import sys
import argparse


######################################################################################
#
# Cylindrical cell model : MD mimic Deserno et al Macromolecules 2000,33, 199-206
#
# Generate integrated charge distribution around the rod
#
######################################################################################


print( " " )
print( "===================================================" )
print( "=      MD simulation of a charged rod             =" )
print( "===================================================" )
print( " " )

print( "Program Information: \n" )
print( espressomd.code_info.features() )
system = espressomd.System(box_l=[1,1,1])
#system.random

#############################################################
#  Parameters
#############################################################

# line charge density of the rod 
line_dens =  1.0

# Output parameters
energy_filename = "rod-energy_" + str(line_dens) + ".dat"
dist_filename = "rod-dist_" + str(line_dens) + ".dat"
vsf_filename = "rod_" + str(line_dens) + ".vsf"
vcf_filename = "rod_" + str(line_dens) + ".vcf"
pos_filename = "last_" + str(line_dens) + ".npy"
positions_filename = "positions_" + str(line_dens) + ".dat"

#############################################################
# System parameters
#############################################################

# System parameters
# Rod length
L = np.sqrt(np.pi) * 28.2

# rod radius
r0 = 1.0

# Bjerrum length
bjerrum_length = 1.0

# ion diameter
ion_diameter = 1.0

# valency of the counterions
valency_ci = 1

# number of beads of the rod
num_rod_beads = int(round(line_dens * L))


# Run parameters
#############################################################

# MD frames to go		
max_frames = 10000

# number of timesteps per frame
steps_per_frame = 100

# accuracy of p3m algorithm 
accuracy =  1e-3;


####################
# CONSTANTS
####################
# particle types 
rod_type = 1
ci_type = 2

####################
# DERIVED PARAMETERS
####################
# total charge of the rod
total_rod_charge = line_dens*L

# charge per rod bead
rod_charge = -1.0
#rod_charge = -total_rod_charge/num_rod_beads

# distance of beads on the rod
rod_distance = L/num_rod_beads

# number of counterions
num_ci = int(round(total_rod_charge/valency_ci))


####################
# ESPResSo system parameters
####################

# simulation box size
system.box_l = [L, L, L]

# coupling to heat bath via dissipation-fluctutation theory
system.thermostat.set_langevin(kT=1.0, gamma=0.5, seed=3)

# timestep in simulation units	
system.time_step = 0.01

# Verlet skin
system.cell_system.skin = 0.4


#
# SETTING UP THE INTERACTIONS
# 
# ion-ion interaction
lj0_eps = 1.0
lj0_sig = ion_diameter
lj0_cut =  2**(1./6.) * ion_diameter
lj0_shift = 0.25
lj0_off = 0

system.non_bonded_inter[ci_type,ci_type].lennard_jones.set_params(
      epsilon=lj0_eps, sigma=lj0_sig,
      cutoff=lj0_cut, shift=lj0_shift, offset=lj0_off)

# ion-rod interaction
lj1_eps = 1.0
lj1_sig = ion_diameter
lj1_cut =  2**(1./6.) * ion_diameter
lj1_shift = 0.25
lj1_off = r0-ion_diameter
print ("Offset for the ion-rod interaction = {}".format(lj1_off))

system.non_bonded_inter[ci_type,rod_type].lennard_jones.set_params(
      epsilon=lj0_eps, sigma=lj0_sig,
      cutoff=lj0_cut, shift=lj1_shift, offset=lj1_off)


####################
# Add beads to the system
####################

    
# CREATE ROD: set rod along the center of the box
px = L/2.
py = L/2.
pz = 0.0

print( "Placing rod beads... ")

for p in range(num_rod_beads):
    # the fix=[1,1,1] option fixes the particles in x,y,z.
    # in other words, these will NOT move from these positions!
    system.part.add(id=p, pos=[px,py,pz], type=rod_type, q=rod_charge, fix=[1,1,1])
    pz = pz + rod_distance
print( "Done." )

# CREATE COUNTERIONS: randomly placed beads in the box
start_pid = num_rod_beads
print( "Placing counter ions... ")

for p in range(num_ci):
    system.part.add(id=start_pid+p, pos=np.random.random(3) * system.box_l, type=ci_type, q=valency_ci)
print( "Done." )

    



#############################################################
#  Warm-up Integration (with capped LJ-potential)           #
#  NOTE: this just gets rid of bad overlaps                 #
#        you still need a proper warmup afterwards          #
#############################################################
 

dist=system.analysis.min_dist(p1=[ci_type],p2=[ci_type])
cap=1
system.force_cap = cap
print( "Warming up... ",)
sys.stdout.flush()
for t in range(4000):
    print( "Warming step: {} min_dist={} cap={}\r".format(t, dist, cap))
    system.integrator.run(200)
    dist=system.analysis.min_dist(p1=[ci_type],p2=[ci_type])
    cap = cap + 1.
    system.force_cap = cap
    if (dist >= ion_diameter):
        break
print( "Done.")
sys.stdout.flush()
print( "Remove capping of LJ-interactions." )
system.force_cap = 0
print( "Warmup finished.")
sys.stdout.flush()



####################
# SETTING UP THE COULOMB INTERACTION
####################

print("Starting P3M. Need to tune...")
sys.stdout.flush()
p3m = electrostatics.P3M(prefactor=bjerrum_length, accuracy=accuracy)
system.actors.add(p3m)
# run for a while to test stability
system.integrator.run(10)
sys.stdout.flush()




# Prepare output of VTF file (for VMD)
#############################################################
vsf_file = open(vsf_filename, 'w')
vcf_file = open(vcf_filename, 'w')

vtf.writevsf(system, vsf_file, types='all')
#vtf.writevsf(system, vsf_file, types=ci_type)
vsf_file.close()

energy_file = open(energy_filename, 'w')
energy_file.write("#time \t coulomb \n")

positions_file = open(positions_filename, 'w')
positions_file.write("#x \t y\n")

#############################################################
#      Integration                                          #
#############################################################



print("Starting simulation...")
system.time=0

for t in range(max_frames):

    print("frame: {}/{}\r".format(t, max_frames),)
    sys.stdout.flush()

    # run run the simulation for a few steps
    system.integrator.run(steps_per_frame)
    
    energy = system.analysis.energy()
    energy_file.write("{}  \t {} \n".format(system.time, energy['coulomb']))
    energy_file.flush()

    # off-line visualisation
    vtf.writevcf(system, vcf_file, types='all')
    vcf_file.flush()
    # note that you can also load the trajectroy file for analysis using:
    # data=np.loadtxt(vcf_filename, comments="t")

    # you can open the particle positions from the VCF file,
    # or alternatively, you can save information here
    # for later analysis
    for p in system.part:
        if p.type==ci_type:
            positions_file.write("{}  \t {} \n".format(p.pos_folded[0], p.pos_folded[1]))

vcf_file.close()
energy_file.close()
positions_file.close()

print( "\nFinished with simulation:" )
