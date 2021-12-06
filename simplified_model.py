import numpy as np
from scipy import optimize
from itertools import chain
import random
import multiprocessing as mp
import random
import itertools as it
import functools as ft

def calculate_Rinf(L,R0,MU,Sinit,initialState,n=None):
	if n is None:
		n = L.shape[1]
	if len(R0)==1:
		R0 = R0*np.ones((n,1))
	def Rinf_function(x):
		x.shape = (n,1)
		y=x - Sinit + Sinit*np.exp(-x*R0/Sinit)
		y.shape = (n,)
		return(y)
	Rinf = optimize.fsolve(Rinf_function,(Sinit-1))
	Rinf = np.array(Rinf)
	Rinf.shape = (n,1)
	return(Rinf)

def simulate(L,R0,MU,Sinit,initialState,Rinf=None):
	n = L.shape[1]
	if len(R0)==1:
		R0 = R0*np.ones((n,1))
	if len(MU)==1:
		MU = MU*np.ones((n,1))

	if Rinf is None:
		Rinf = calculate_Rinf(L,R0,MU,Sinit,initialState,n=n)
	rinf = np.broadcast_to(Rinf,(n,n))

	jump_rate = L*rinf*np.broadcast_to(np.transpose(np.maximum(0,(1-1/R0))),(n,n))

	state = initialState.copy()
	state.shape = (n,1)
	jumps = list()
	count = 0

	while (np.any(state==1)):
		RATE_EXPORT = jump_rate*np.broadcast_to(state%2,(n,n))\
		*np.broadcast_to(np.transpose(state)==0,(n,n))
		RATE_RECOVER = MU*(state%2)

		event_number = random.choices(np.arange(n+n**2), weights = 
			np.concatenate((RATE_RECOVER.flatten(),RATE_EXPORT.flatten())))
		event_number = event_number[0]

		if event_number<n:
			state[event_number] = 2
		else:
			state[(event_number-n)%n] = 1
			jumps.append([(event_number-n)//n,(event_number-n)%n])

		count+=1

		if count==1000000:
			print('while loop terminated after 1,000,000 iterations')

	return(np.array(jumps))

def simulate_time_rescaled(L,R0,MU,Sinit,initialState,Rinf=None):
	n = L.shape[1]
	if len(R0)==1:
		R0 = R0*np.ones((n,1))
	if len(MU)==1:
		MU = MU*np.ones((n,1))

	if Rinf is None:
		Rinf = calculate_Rinf(L,R0,MU,Sinit,initialState,n=n)
	rinf = np.broadcast_to(Rinf,(n,n))

	jump_rate = L*(1/MU)*rinf*np.broadcast_to(np.transpose(np.maximum(0,(1-1/R0))),(n,n))

	state = initialState.copy()
	state.shape = (n,1)
	jumps = list()
	count = 0

	while (np.any(state==1)) and count<1000000:
		RATE_EXPORT = jump_rate*np.broadcast_to(state%2,(n,n))\
		*np.broadcast_to(np.transpose(state)==0,(n,n))
		RATE_RECOVER = (state%2)

		event_number = random.choices(np.arange(n+n**2), weights = 
			np.concatenate((RATE_RECOVER.flatten(),RATE_EXPORT.flatten())))
		event_number = event_number[0]

		if event_number<n:
			state[event_number] = 2
		else:
			state[(event_number-n)%n] = 1
			jumps.append([(event_number-n)//n,(event_number-n)%n])

		count+=1

		if count==1000000:
			print('while loop terminated after 1,000,000 iterations')

	return(np.array(jumps))

def importation_function(sim_params):
	jumps=simulate(*sim_params)
	n=(sim_params[0]).shape[1]
	if np.any(jumps):
		importation = np.array([np.any(jumps[:,1]==j) for j in range(n)])
	else:
		importation = np.zeros(n)
	return importation

def importation_function_time_rescaled(sim_params):
	jumps=simulate_time_rescaled(*sim_params)
	n=(sim_params[0]).shape[1]
	if np.any(jumps):
		importation = np.array([np.any(jumps[:,1]==j) for j in range(n)])
	else:
		importation = np.zeros(n)
	return importation

def avg_importation(sim_params,n_sims,time_rescaled=False):
	n=(sim_params[0]).shape[1]

	if time_rescaled:
		impfunc = importation_function_time_rescaled
	else:
		impfunc = importation_function

	pool=mp.Pool(mp.cpu_count())
	y=pool.imap_unordered(impfunc,
		[sim_params for i in range(n_sims)])
	pool.close()

	importation = np.sum(np.array(list(y)),axis=0)/n_sims
	importation[sim_params[4][:,0]==1] = 1
	importation.shape = (n,1)

	return importation

def expected_infected(sim_params,n_sims,time_rescaled=False):
	n=(sim_params[0]).shape[1]
	Rinf=calculate_Rinf(*sim_params)
	cpus = mp.cpu_count()

	if time_rescaled:
		impfunc = importation_function_time_rescaled
	else:
		impfunc = importation_function

	pool=mp.Pool(cpus)
	y=pool.imap_unordered(impfunc,
		[sim_params for i in range(n_sims)],chunksize=int(n_sims/cpus))
	pool.close()

	importation = np.sum(np.array(list(y)),axis=0)/n_sims
	importation[sim_params[4][:,0]==1] = 1
	importation.shape = (n,1)
	E_of_i = importation*Rinf

	return E_of_i

def infected_summary(sim_params,n_sims,time_rescaled=False):
	n=(sim_params[0]).shape[1]
	Rinf=calculate_Rinf(*sim_params)
	cpus = mp.cpu_count()

	if time_rescaled:
		impfunc = importation_function_time_rescaled
	else:
		impfunc = importation_function

	pool=mp.Pool(cpus)
	y=pool.imap_unordered(impfunc,
		[sim_params for i in range(n_sims)],chunksize=int(n_sims/cpus))
	pool.close()

	x = np.array(list(y))
	x[:,sim_params[4][:,0]==1] = 1
	outbreak_sizes = np.sum(x*Rinf.T,axis=1)

	mean_infections = np.mean(outbreak_sizes)
	median_infections = np.median(outbreak_sizes)
	upper_quartile = np.percentile(outbreak_sizes,75)
	lower_quartile = np.percentile(outbreak_sizes,25)

	return(mean_infections,median_infections,upper_quartile,lower_quartile)

def infected_sizes(sim_params,n_sims,time_rescaled=False):
	n=(sim_params[0]).shape[1]
	Rinf=calculate_Rinf(*sim_params)
	cpus = mp.cpu_count()

	if time_rescaled:
		impfunc = importation_function_time_rescaled
	else:
		impfunc = importation_function

	pool=mp.Pool(cpus)
	y=pool.imap_unordered(impfunc,
		[sim_params for i in range(n_sims)],chunksize=int(n_sims/cpus))
	pool.close()

	x = np.array(list(y))
	x[:,sim_params[4][:,0]==1] = 1
	outbreak_sizes = np.sum(x*Rinf.T,axis=1)

	return(outbreak_sizes)