import espressomd
from espressomd import interactions
from espressomd import polymer
import numpy as np
import sys
import argparse


#############################################################
#                                                           #
#                    Ideal Chain Polymer                    #
#                                                           #
#############################################################

print( " " )
print( "===================================================" )
print( "=              Ideal Chain Polymer                =" )
print( "===================================================" )
print( " " )

print( "Program Information: \n" )
print( espressomd.code_info.features() )
system = espressomd.System(box_l = 3*[1])

parser = argparse.ArgumentParser(description='Read simulation parameters')
parser.add_argument('-name', type=str,  help='Base name of output files', default="output")
parser.add_argument('-MPC', type=int,  help='Monomers per chain', default=10)
parser.add_argument('--with_LJ', action = 'store_true')
args = parser.parse_args()


#############################################################
#  Parameters                                               #
#############################################################

# Start of all files written by this script
name = args.name


# System parameters
#############################################################

# number of polymer chains
N_polymers = 1
# number of monomers (beads) in the chain
beads_per_chain = args.MPC
# initial polymer bond length
bond_l =        1.0
# density
density =       0.0005


# Interaction parameters
#############################################################

# repulsive Lennard Jones
lj1_eps =          1.0
lj1_sig =          1.0
lj1_cut =          2**(1./6.)

# attractive harmonic spring
k =             100.0
      
# Integration parameters
#############################################################

system.time_step = 0.01

# tuning parameter for frequency of Verlet rebuilds
system.cell_system.skin = 0.4

# coupling to heat bath via dissipation-fluctutation theory
kT = 1.
gamma = 1.
system.thermostat.set_langevin(kT=kT, gamma=gamma, seed=42)

# warmup integration (with capped LJ potential) until particles are at least min_dist apart
# see below for what the values do
warm_step =   6000
warm_loop =   5
warm_cap1 =   10
warm_incr =   10
min_dist =    0.90

# integration (with full LJ potential) for int_time
int_time = 10000
loop_time = 10



#############################################################
#  Setup System                                             #
#############################################################

n_part = N_polymers*beads_per_chain
box_l = N_polymers * beads_per_chain*bond_l
system.box_l = [box_l, box_l, box_l]

# Interaction setup
#############################################################



if args.with_LJ:
    raise NotImplementedError('fill in your code for interaction setup')

    # YOUR CODE HERE
    # look up in the espresso documentation how to set up a LJ interaction between
    # all particles of type 0 using the lj1_* variables defined above

harmonic_bond = interactions.HarmonicBond(k=k, r_0=bond_l)
system.bonded_inter.add(harmonic_bond)

print( f"\nSimulate {N_polymers} polymer chains with {beads_per_chain} monomers of bond-length {bond_l} each")
print( f"in a cubic simulation box of length {box_l} at density {density} with gamma {gamma} and temp {kT}.")
print( f"LJ-interaction: {args.with_LJ}")

# Particle setup
#############################################################

print( "Creating polymers... ",)
sys.stdout.flush()

polymer_positions = polymer.positions(n_polymers = N_polymers,
                                      beads_per_chain = beads_per_chain,
                                      bond_length=bond_l,
                                      seed = 42)
#polymer_positions: list of length N_polymers, with each elements containing the list of beads_per_chain monomer positions
for p in polymer_positions:
  for i,pos in enumerate(p):
    id = len(system.part)
    system.part.add(id=id, pos=pos)
    if i > 0:
      system.part[id].add_bond((harmonic_bond, id - 1))

#############################################################
#  Warm-up Integration (with capped LJ-potential)           #
#############################################################
if args.with_LJ:
    print( f"\nRemoving overlap by successively removing force capping for {warm_step*warm_loop} timesteps in {warm_loop} loops; ")
    system.time = 0
    cap = warm_cap1

    for t in range(warm_loop):
        system.force_cap = cap
        system.integrator.run(warm_step)
        dist = system.analysis.min_dist()
        print(f'min dist between monomers: {dist}')
    
        cap = cap + warm_incr
        
    print( "Remove capping of LJ-interactions." )
    system.force_cap = 0
    
    
# overlaps are taken care of, now run for a while to properly warm-up
warm_steps = int(0.05*int_time/system.time_step)
print(f'warming up for {warm_steps} steps')
system.integrator.run(warm_steps)


#############################################################
#      Integration                                          #
#############################################################

system.time = 0.
step = 0

loop_steps = int(loop_time/system.time_step)

# arrays to hold values as the simulation progresses
re2_ar = []
rg2_ar = []
t_ar  = []

print( f"\nStart integration with timestep {system.time_step} until time t>={int_time}",flush=True) 

obs_file = open("{}-obs.dat".format(name), "w")
obs_file.write("#t\tmin_dist\tre\trg\tT\n" )
obs_file.flush()
 
while system.time<int_time:
    system.integrator.run(loop_steps)

    #print( f" Step {system.time/system.time_step}/{int_step*int_loop} (t={system.time}): ", flush= True)
    
    dist = system.analysis.min_dist()
    temp = system.analysis.energy()['kinetic']/(n_part/(2./3.))
    re = system.analysis.calc_re( chain_start=0, number_of_chains=N_polymers, chain_length=beads_per_chain)[0]
    rg = system.analysis.calc_rg( chain_start=0, number_of_chains=N_polymers, chain_length=beads_per_chain)[0]
    re2_ar.append(re**2)
    rg2_ar.append(rg**2)
    t_ar.append(system.time)
    
    obs_file.write(f'{system.time}\t{dist}\t{re}\t{rg}\t{temp}\n')
    obs_file.flush()

    #print(F"rg={rg}, re={re}, T={temp} ..."),

obs_file.close()

print( "\nFinished with simulation:" )

# derive time averages
re2_av = np.mean(re2_ar)
rg2_av = np.mean(rg2_ar)

print( f"<re^2> = {re2_av}")
print( f"<rg^2> = {rg2_av}")
print( f"<re^2>/<rg^2> = {re2_av/rg2_av}")

print( "Done." )
