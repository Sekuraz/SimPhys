import numpy as np
import matplotlib.pyplot as plt

Nsamples = 400000

def metropolis(N, P, trial_move, phi0, dx):
    ret = []
    ret.append(phi0)
    phi = phi0
    phi_next = phi0
    counter = 0
    for i in range(N):
        #Perform trial move
        phi_next = trial_move(phi, dx)
        
        #Generate a random number
        r = np.random.random_sample()
        
        #Decide if new state is accepted
        if r < min(1, P(phi_next)/P(phi)):
            phi = phi_next
            counter += 1
        
        #Add new state to the list
        ret.append(phi)
    #Calculate the acceptance rate
    acceptance_rate = counter / N
    return ret, acceptance_rate
    
def gaussian_distribution(x):
    ret = np.exp(-x ** 2) / np.sqrt(np.pi)
    return ret

def trial_move(x, dx):
    r = 2 * dx * np.random.random_sample() - dx
    ret = x + r
    return ret

width = 5.787
height = width*0.6
plt.rc('figure', figsize=(width,height))
plt.rc('text', usetex=True)

dx = np.logspace(-1.0, 2.0, num = 4)
for i in range(len(dx)):
    samples, acceptance_rate = metropolis(Nsamples, gaussian_distribution, trial_move, 0.0, dx[i])
    x = np.linspace(-5, 5, 1000)
    plt.hist(samples, bins=100, density = True, label=r'$\mathrm{histogram}$')
    plt.plot(x, gaussian_distribution(x), label=r'$p(x) = \frac{\exp(-x^2)}{\sqrt{2}}$')
    plt.legend()
    plt.xlabel(r'$x$')
    plt.xlim(-5.0, 5.0)
    plt.tight_layout()
    plt.show()
    print(dx[i])
    print(acceptance_rate)
