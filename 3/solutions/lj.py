#!/usr/bin/env python3


class Simulation:
    def __init__(self, dt, x, v, box, r_cut, shift, f_max):
        self.dt = dt
        self.x = x.copy()
        self.v = v.copy()
        self.box = box.copy()
        self.r_cut = r_cut
        self.shift = shift
        self.f_max = f_max

        self.n_dims = self.x.shape[0]
        self.n = self.x.shape[1]
        self.f = np.zeros_like(x)

        # both r_ij_matrix and f_ij_matrix are computed in self.forces()
        self.r_ij_matrix = np.zeros((self.n, self.n, self.n_dims))
        self.f_ij_matrix = np.zeros((self.n, self.n, self.n_dims))
        # computed in e_pot_ij_matrix
        self.e_pot_ij_matrix = np.zeros((self.n, self.n))

    def distances(self):
        self.r_ij_matrix = np.repeat([self.x.transpose()], self.n, axis=0)
        self.r_ij_matrix -= np.transpose(self.r_ij_matrix, axes=[1, 0, 2])
        # minimum image convention
        image_offsets = self.r_ij_matrix.copy()
        for nth_box_component, box_component in enumerate(self.box):
            image_offsets[:, :, nth_box_component] = \
                np.rint(image_offsets[:, :, nth_box_component] / box_component) * box_component
        self.r_ij_matrix -= image_offsets

    def energies(self):
        r = np.linalg.norm(self.r_ij_matrix, axis=2)
        with np.errstate(all='ignore'):
            self.e_pot_ij_matrix = np.where((r != 0.0) & (r < self.r_cut),
                                            4.0 * (np.power(r, -12.) - np.power(r, -6.)) + self.shift, 0.0)

    def forces(self):
        # first update the distance vector matrix, obeying minimum image convention
        self.distances()
        self.f_ij_matrix = self.r_ij_matrix.copy()
        r = np.linalg.norm(self.r_ij_matrix, axis=2)
        with np.errstate(all='ignore'):
            fac = np.where((r != 0.0) & (r < self.r_cut),
                           4.0 * (12.0 * np.power(r, -13.) - 6.0 * np.power(r, -7.)), 0.0)
        for dim in range(self.n_dims):
            with np.errstate(invalid='ignore'):
                self.f_ij_matrix[:, :, dim] *= np.where(r != 0.0, fac / r, 0.0)
        f = np.sum(self.f_ij_matrix, axis=0).transpose()
        self.f = np.clip(f, - self.f_max, self.f_max)

    def e_pot(self):
        return np.sum(self.e_pot_ij_matrix) / 2

    def e_kin(self):
        return np.sum(0.5 * np.power(self.v, 2))  # mass = 1

    def energy(self):
        """Compute and return the energy components of the system."""
        # compute energy matrix
        self.energies()
        # TODO compute interaction energy from self.e_pot_ij_matrix
        # TODO calculate kinetic energy from the velocities self.v and return both energy components

        return np.array((self.e_pot(), self.e_kin()))

    def temperature(self):
        # we use kB=1 for the scaling of our units
        return 2 * self.e_kin() / (self.n_dims * self.n)

    def pressure(self):
        def ideal():
            return self.e_kin() * 2 / (self.n_dims * np.product(self.box))
        # s = 0
        # for i in range(1, self.n):
        #     for j in range(i):
        #         s += np.dot(self.f_ij_matrix[i][j], self.r_ij_matrix[i][j])
        # f = s * 2

        f = np.multiply(self.f_ij_matrix, self.r_ij_matrix)
        f = np.sum(f)

        area = np.product(self.box)
        ret = 1 / area * (self.e_kin() + f * 0.25)
        # print("Pressure: calculated=%4.6f, ideal=%4.6f" % (ret, ideal()))
        return ret

    def rdf(self):
        self.distances()
        r = np.linalg.norm(self.r_ij_matrix, axis=2)
        hist, bins = np.histogram(r, bins=100, range=(0.8, 5))
        return hist

    def rescale(self, T0):
        temp = self.temperature()
        factor = np.sqrt(T0 / temp)
        self.v *= factor

    def propagate(self):
        # update positions
        self.x += self.v * self.dt + 0.5 * self.f * self.dt * self.dt

        # half update of the velocity
        self.v += 0.5 * self.f * self.dt

        # compute new forces
        self.forces()
        # we assume that all particles have a mass of unity

        # second half update of the velocity
        self.v += 0.5 * self.f * self.dt


def write_checkpoint(state, path, overwrite=False):
    if os.path.exists(path) and not overwrite:
        raise RuntimeError("Checkpoint file already exists")
    with open(path, 'wb') as fp:
        pickle.dump(state, fp)


if __name__ == "__main__":
    import argparse
    import pickle
    import itertools
    import logging

    import os.path

    import numpy as np
    import scipy.spatial  # todo: probably remove in template
    import tqdm

    logging.basicConfig(level=logging.INFO)

    parser = argparse.ArgumentParser()
    parser.add_argument(
        'N_per_side',
        type=int,
        help='Number of particles per lattice side.')
    parser.add_argument(
        '--cpt',
        type=str,
        help='Path to checkpoint.')
    parser.add_argument(
        '-t',
        '--temperature',
        type=float,
        help='Target temperature.')
    parser.add_argument(
        '-f',
        '--fmax',
        type=float,
        help='Enable force capping.')
    parser.add_argument(
        '-r',
        '--random_initial_positions',
        help='Start with random initial positions.')
    args = parser.parse_args()

    np.random.seed(2)

    DT = 0.01
    T_MAX = 1000.0
    N_TIME_STEPS = int(T_MAX / DT)

    R_CUT = 2.5
    SHIFT = 0.016316891136

    DIM = 2
    DENSITY = 0.316
    N_PER_SIDE = args.N_per_side
    N_PART = N_PER_SIDE**DIM
    VOLUME = N_PART / DENSITY
    BOX = np.ones(DIM) * VOLUME**(1. / DIM)

    SAMPLING_STRIDE = 3

    if not args.cpt or not os.path.exists(args.cpt):
        logging.info("Starting from scratch.")
        # particle positions
        if args.random_initial_positions:
            x = np.random.random((DIM, N_PART))
            for i in range(DIM):
                x[i] *= BOX[i]
        else:
            x = np.array(list(itertools.product(np.linspace(0, BOX[0], N_PER_SIDE, endpoint=False),
                                             np.linspace(0, BOX[1], N_PER_SIDE, endpoint=False)))).T

        #

        # random particle velocities
        v = 0.5*(2.0 * np.random.random((DIM, N_PART)) - 1.0)

        positions = []
        energies = []
        pressures = []
        temperatures = []
        rdfs = []
        potential_energies = []
        kinetic_energies = []
        if args.fmax:
            f_max = args.fmax
        else:
            f_max = np.inf
    #elif args.cpt and os.path.exists(args.cpt):
    else:
        logging.info("Reading state from checkpoint.")
        with open(args.cpt, 'rb') as fp:
            state = pickle.load(fp)
        positions = state['positions']
        energies = state['energies']
        pressures = state['pressures']
        temperatures = state['temperatures']
        rdfs = state['rdfs']
        potential_energies = state['potential_energies']
        kinetic_energies = state['kinetic_energies']

        x = state['x']
        v = state['v']
        f = state['f']
        f_max = state['f_max']

    sim = Simulation(DT, x, v, BOX, R_CUT, SHIFT, f_max)

    # If checkpoint is used, also the forces have to be reloaded!
    if args.cpt and os.path.exists(args.cpt):
        sim.f = f

    for i in tqdm.tqdm(range(N_TIME_STEPS)):
        sim.propagate()

        if i % SAMPLING_STRIDE == 0:
            positions.append(sim.x.copy())
            pressures.append(sim.pressure())
            energies.append(np.sum(sim.energy()))
            temperatures.append(sim.temperature())
            rdfs.append(sim.rdf())
            potential_energies.append(sim.e_pot())
            kinetic_energies.append(sim.e_kin())

            if args.temperature:
                sim.rescale(args.temperature)
                
            if args.fmax:
                sim.f_max *= 1.1
                if sim.f_max >= 1e9:
                    sim.f_max = np.inf

    if args.cpt:
        state = {'positions': positions,
                 'energies': energies,
                 'pressures': pressures,
                 'temperatures': temperatures,
                 'rdfs': rdfs,
                 'potential_energies': potential_energies,
                 'kinetic_energies': kinetic_energies,
                 'x':sim.x,
                 'v': sim.v,
                 'f': sim.f,
                 'f_max': sim.f_max}
        write_checkpoint(state, args.cpt, overwrite=True)
