import numpy as np
import matplotlib.pyplot as plt
import itertools
from scipy import signal

N_TEMPERATURE = 41

def calculate_energy(spin_config, Nx, Ny):
    ret = 0.0
    for nx in range(Nx):
        for ny in range(Ny):
            ret += -0.5 * spin_config[nx * Ny + ny] * (spin_config[((nx - 1) % Nx) * Ny + ny]\
                + spin_config[((nx + 1) % Nx) * Ny + ny] + spin_config[nx * Ny + (ny - 1) % Ny]\
                + spin_config[nx * Ny + (ny + 1) % Ny])
    return ret

def calculate_magnetization(spin_config):
    return np.sum(spin_config)

def ising_mean(Nx, Ny, T):
    Z = 0.0
    mean_energy = 0.0
    mean_magnetization = 0.0
    iterator = itertools.product([-1, +1], repeat = Nx*Ny)
    for i in iterator:
        m = calculate_magnetization(list(i))
        E = calculate_energy(list(i), Nx, Ny)
        boltzmann = np.exp( -E / T)
        Z += boltzmann
        mean_energy += E * boltzmann
        mean_magnetization += np.abs(calculate_magnetization(list(i))) * boltzmann
    mean_energy /= (Z * Nx * Ny)
    mean_magnetization /= (Z * Nx * Ny)
    return mean_energy, mean_magnetization

def ising_monte_carlo(N_MC_STEPS, Nx, Ny, T):
    N_SPINS = Nx * Ny
    spin_config = [(-1) ** np.random.randint(0, high=2) for i in range(N_SPINS)]
    energy = calculate_energy(spin_config, Nx, Ny)
    magnetization = calculate_magnetization(spin_config)
    mean_energy = 0.0
    mean_magnetization = 0.0
    
    energy_time_series = []
    magnetization_time_series = []
    
    for i in range(N_MC_STEPS):
        #Choose a random spin to flip
        random_spin = np.random.randint(0, high=N_SPINS)
        
        #Calculate the coordinates of the spin
        ny = random_spin % Ny
        nx = int((random_spin - ny) / Ny)
        
        #Calculate the new energy for this configuration
        energy_trial = energy + 2 * spin_config[nx * Ny + ny] * (spin_config[((nx - 1) % Nx) * Ny + ny]\
                + spin_config[((nx + 1) % Nx) * Ny + ny] + spin_config[nx * Ny + (ny - 1) % Ny]\
                + spin_config[nx * Ny + (ny + 1) % Ny])
        
        #Decide if new spin configuration is accepted
        if np.random.random_sample() < min(1, np.exp( (energy - energy_trial)/ T)):
            spin_config[random_spin] *= -1
            energy = energy_trial
            magnetization += 2 * spin_config[random_spin]
        
        energy_time_series.append(energy)
        magnetization_time_series.append(magnetization)
        
        mean_energy += energy
        mean_magnetization += np.abs(magnetization)
        
    mean_energy /= (N_MC_STEPS * N_SPINS)
    mean_magnetization /= (N_MC_STEPS * N_SPINS)
    return mean_energy, mean_magnetization, np.array(energy_time_series), np.array(magnetization_time_series)

#def error_analysis(O):
#    #calculate the mean value of the observables
#    mean_value = np.mean(O)
#    variance = np.mean((O - mean_value) ** 2)
#    
#    N = len(O)
#    
#    ac = signal.correlate((O - mean_value), (O - mean_value), mode='same')
#    
#    plt.plot(ac)
#    plt.show()
#    
#    tau_int = 0.5 * np.sum(ac) / max(ac)
#    
#    #calculate the effective statistics
#    N_eff = N / (2 * tau_int)
#    
#    #calculate the error of the observables mean value
#    error = np.sqrt(variance / N_eff)
#    
#    return mean_value, tau_int, N_eff, error

#def bivariate_gaussian(N, rho):
#    e = np.random.normal(size=N)
#    for i in range(N):
#        e[i] = rho * e[i-1] + np.sqrt(1 - rho ** 2) * e[i]
#    return e

#Test the if error_analysis returns the correct autocorrelation time
#N_TEST = 1000000
#rho = 0.99005
#e = bivariate_gaussian(N_TEST, rho)
#tau_exact = 0.5 * (1 + rho) / (1 - rho)
#tau_error_analysis = error_analysis(e)
#print(tau_exact)
#print(tau_error_analysis)


E_exact = np.zeros(N_TEMPERATURE)
M_exact = np.zeros(N_TEMPERATURE)

E_MC = np.zeros(N_TEMPERATURE)
M_MC = np.zeros(N_TEMPERATURE)

#energy_time_series = np.zeros([N_TEMPERATURE, 10000])
#magnetization_time_series = np.zeros([N_TEMPERATURE, 10000])

for i in range(N_TEMPERATURE):
    T = 1.0 + i*0.1
#    E_exact[i], M_exact[i] = ising_mean(4, 4, T)
    E_MC[i], M_MC[i], energy_time_series[i,:], magnetization_time_series[i,:] = ising_monte_carlo(10000, 4, 4, T)
#    print(E_exact[i])
    

#Analyse the time series from the Monte Carlo simulation
#for i in range(N_TEMPERATURE):
#    statistics_energy[i] = error_analysis(energy_time_series[i,:])
#    statistics_magnetization[i] = error_analysis(magnetization_time_series[i,:])


#Create Plots
T_range = np.linspace(1.0, 5.0, num = 41, endpoint=True)

width = 5.787
height = width*0.6
plt.rc('figure', figsize=(width,height))
plt.rc('text', usetex=True)

#plt.plot(T_range, E_exact, label=r'exact summation')
plt.plot(T_range, E_MC, label=r'Monte Carlo')
plt.xlabel(r'$\mathrm{dimensionless}$ $\mathrm{temperature}$ $\frac{k_\mathrm{B}T}{J}$')
plt.ylabel(r'$\mathrm{mean}$ $\mathrm{energy}$ $\frac{E}{J}$')
plt.tight_layout()
plt.legend()
plt.show()

#plt.plot(T_range, M_exact, label=r'exact summation')
plt.plot(T_range, M_MC, label=r'Monte Carlo')
plt.xlabel(r'$\mathrm{dimensionless}$ $\mathrm{temperature}$ $\frac{k_\mathrm{B}T}{J}$')
plt.ylabel(r'$\mathrm{mean}$ $\mathrm{magnetization}$ $\frac{\left|\mu\right|}{m}$')
plt.tight_layout()
plt.legend()
plt.show()
