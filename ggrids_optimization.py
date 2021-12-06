import numpy as np
import scipy as sp
import itertools as it
import multiprocessing as mp
from simplified_model import expected_infected

def allocation_brute_force(n_pops,sim_params,n_sims,time_rescaled=False):
	Sinit = sim_params[3].copy()
	n = Sinit.shape[0]
	R0 = sim_params[1].copy()
	if len(R0)==1:
		R0 = R0*np.ones((n,1))
	herd_threshold = np.ceil(Sinit*(1-1/R0))
	initialState = sim_params[4].copy()

	allocations = list(it.combinations(np.arange(n)[np.where(initialState==0)[0]],n_pops))

	sizes = np.zeros((len(allocations),n))

	for i in range(len(allocations)):
		strategy_Sinit = Sinit.copy()
		for patch in allocations[i]:
			strategy_Sinit[patch] -= herd_threshold[patch]
		strategy_R0 = R0*strategy_Sinit/Sinit
		strategy_params = (sim_params[0],strategy_R0,
			sim_params[2],strategy_Sinit,sim_params[4])
		sizes[i,:] = expected_infected(strategy_params,n_sims,time_rescaled)[:,0]

	return(allocations,sizes)

def allocation_brute_force_optimize(n_pops,sim_params,n_sims,time_rescaled=False):
	allocations, sizes = allocation_brute_force(n_pops,sim_params,n_sims,time_rescaled)
	total_size = np.sum(sizes,axis=1)
	optimal_index = np.argmin(total_size)
	return(allocations[optimal_index], sizes[optimal_index])