import espressomd
from espressomd import shapes
from espressomd import lb
from espressomd import lbboundaries

import numpy as np


print( " " )
print( "===================================================" )
print( "=            Lattice-Boltzmann Fluid              =" )
print( "===================================================" )
print( " " )

print( "Program Information: \n" )
print( espressomd.code_info.features() )

# geometry
box_l = 32.
padding = 1.

# fluid parameters
LB_params = {'agrid':1.,
             'dens':1.,
             'visc':1.,
             'tau':0.01,
             'ext_force_density':[0.01, 0., 0.],
             'kT':0.}
            
system = espressomd.System(box_l = 3*[box_l])
system.time_step = LB_params['tau']
system.cell_system.skin = 0.2



# choose between these two: GPU or CPU (depending on compilation features)
if espressomd.espressomd.cuda_init.gpu_available():
    lbf = lb.LBFluidGPU(**LB_params)
else:
    lbf = lb.LBFluid(**LB_params)


system.actors.add(lbf)

# create the boundary "shape"
upper_wall=shapes.Wall(normal=[0,1,0], dist=padding)
lower_wall=shapes.Wall(normal=[0,-1,0], dist=-(box_l-padding))

# from these shapes, define the LB boundary
upper_bound=lbboundaries.LBBoundary(shape=upper_wall)
lower_bound=lbboundaries.LBBoundary(shape=lower_wall)

system.lbboundaries.add(upper_bound)
system.lbboundaries.add(lower_bound)

#system.part.add(pos=0.5*system.box_l, type=0)

probe_ys = np.linspace(padding, box_l-padding, num = 200)

max_time=500
for t in range(max_time):
    system.integrator.run(int(1./system.time_step))
    
    pos = [0,system.box_l[1]/2.,0]
    vel = lbf.get_interpolated_velocity(pos)

    print("time: {} velocity:{}".format(system.time, vel))

outdir = ("./")
#lbf.print_vtk_velocity("{}/velocity.vtk".format(outdir))
print("**Simulation Completed** ")
