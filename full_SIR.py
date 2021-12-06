import numpy as np
import multiprocessing as mp
import random

def SIR_simulation(sim_params):
	A = sim_params[0]
	R0 = sim_params[1]
	MU = sim_params[2]
	N = sim_params[3]
	Iinit = sim_params[4]
	n = A.shape[1]
	BETA = (R0*MU)/N
	if len(BETA)==1:
		BETA = BETA*np.ones((n,1))
	if len(MU)==1:
		MU = MU*np.ones((n,1))

	I = np.array(Iinit)
	S = np.array(N)-I

	while sum(I)!=0:
		cross_infection = np.sum(A*np.tile(I.T,(n,1)),axis=1)
		cross_infection.shape = (n,1)
		RATE_INFECT = BETA*S*(I+cross_infection)
		RATE_RECOVER = MU*I

		event_number = random.choices(np.arange(2*n), weights = \
			np.concatenate((RATE_INFECT.flatten(),RATE_RECOVER.flatten())))
		event_number = event_number[0]

		if event_number < n:
			S[event_number] -= 1
			I[event_number] += 1
		elif event_number>=n:
			I[event_number-n] -= 1

	total_infections = np.sum(N-S)

	return(total_infections)

def outbreaks_summary(sim_params,n_sims,time_rescaled=None):
	n=(sim_params[0]).shape[1]
	N = sim_params[3]
	Iinit = sim_params[4]*N/40
	sim_params_with_Iinit = (sim_params[0],sim_params[1],sim_params[2],sim_params[3]
		,Iinit)
	cpus = mp.cpu_count()

	pool=mp.Pool(cpus)
	y=pool.imap_unordered(SIR_simulation,
		[sim_params for i in range(n_sims)],chunksize=int(n_sims/cpus))
	pool.close()

	x = list(y)
	epidemic_sizes = [size for size in x if size>0.01*np.sum(N*sim_params[4])]
	if len(epidemic_sizes)==0:
		epidemic_sizes=[0]

	mean_infections = np.mean(epidemic_sizes)
	median_infections = np.median(epidemic_sizes)
	upper_quartile = np.percentile(epidemic_sizes,75)
	lower_quartile = np.percentile(epidemic_sizes,25)

	return(mean_infections,median_infections,upper_quartile,lower_quartile)

def outbreaks_sizes(sim_params,n_sims,time_rescaled=None):
	n=(sim_params[0]).shape[1]
	N = sim_params[3]
	Iinit = sim_params[4]*N/40
	sim_params_with_Iinit = (sim_params[0],sim_params[1],sim_params[2],sim_params[3]
		,Iinit)
	cpus = mp.cpu_count()

	pool=mp.Pool(cpus)
	y=pool.imap_unordered(SIR_simulation,
		[sim_params for i in range(n_sims)],chunksize=int(n_sims/cpus))
	pool.close()

	x = list(y)
	epidemic_sizes = [size for size in x if size>0.01*np.sum(N*sim_params[4])]
	if len(epidemic_sizes)==0:
		epidemic_sizes=[0]

	return(epidemic_sizes)

def SIR_importation_simulation(sim_params):
	A = sim_params[0]
	R0 = sim_params[1]
	MU = sim_params[2]
	N = sim_params[3]
	Iinit = sim_params[4]
	n = A.shape[1]
	BETA = (R0*MU)/N
	if len(BETA)==1:
		BETA = BETA*np.ones((n,1))
	if len(MU)==1:
		MU = MU*np.ones((n,1))

	I = np.array(Iinit)
	S = np.array(N)-I

	while sum(I)!=0:
		cross_infection = np.sum(A*np.tile(I.T,(n,1)),axis=1)
		cross_infection.shape = (n,1)
		RATE_INFECT = BETA*S*(I+cross_infection)
		RATE_RECOVER = MU*I

		event_number = random.choices(np.arange(2*n), weights = \
			np.concatenate((RATE_INFECT.flatten(),RATE_RECOVER.flatten())))
		event_number = event_number[0]

		if event_number < n:
			S[event_number] -= 1
			I[event_number] += 1
		elif event_number>=n:
			I[event_number-n] -= 1

	imporation = [N[i]!=S[i] for i in range(n)]

	return(imporation)

def importation_probabilities(sim_params,n_sims,time_rescaled=None):
	n=(sim_params[0]).shape[1]
	N = sim_params[3]
	Iinit = sim_params[4]*N/40
	sim_params_with_Iinit = (sim_params[0],sim_params[1],sim_params[2],sim_params[3]
		,Iinit)
	cpus = mp.cpu_count()

	pool=mp.Pool(cpus)
	y=pool.imap_unordered(SIR_importation_simulation,
		[sim_params_with_Iinit for i in range(n_sims)],chunksize=int(n_sims/cpus))
	pool.close()

	x = np.array(list(y))
	importation = np.sum(x,axis=0)/n_sims
	importation.shape = (n,1)

	return importation