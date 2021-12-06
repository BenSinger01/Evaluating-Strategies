import numpy as np
from simplified_model import expected_infected
from simplified_model import infected_summary
from simplified_model import avg_importation
from simplified_model import infected_sizes
from random_walk_percolation_centrality import random_walk_percolation_centrality as rwpc
from full_SIR import importation_probabilities
import networkx as nx

def pro_rata_priority():
	print('pro_rata_priority function should not actually be called')
	return(None)

def size_given_strategy(priority_of_patches,doses,sim_params,n_sims
	,infections_function=expected_infected,time_rescaled=False):
	Sinit = sim_params[3].copy()
	n = Sinit.shape[0]
	R0 = sim_params[1].copy()
	initialState = sim_params[4].copy()
	if priority_of_patches==pro_rata_priority:
		Sinit_adj = Sinit*(1-initialState)
		allocation = Sinit_adj*doses/np.sum(Sinit_adj)
	else:
		if len(R0)==1:
			R0 = R0*np.ones((n,1))
		herd_threshold = np.ceil(Sinit*(1-1/R0))

		if callable(priority_of_patches):
			priority_order = priority_of_patches(sim_params,infections_function,n_sims
				,time_rescaled=time_rescaled)
		else:
			priority_order = priority_of_patches

		allocation = np.zeros((n,1))
		doses_temp = doses
		for patch in priority_order:
			if initialState[patch]==1:
				continue
			elif doses_temp >= herd_threshold[patch]:
				allocation[patch] = herd_threshold[patch]
				doses_temp -= herd_threshold[patch]
			elif doses_temp < herd_threshold[patch]:
				allocation[patch] = doses_temp
				doses_temp = 0

	strategy_R0 = R0*(Sinit-allocation)/Sinit
	strategy_Sinit = Sinit-allocation
	strategy_params = (sim_params[0],strategy_R0,
		sim_params[2],strategy_Sinit,sim_params[4])
	size = infections_function(strategy_params,n_sims,time_rescaled=time_rescaled)
	
	return(allocation,size)

def summary_given_strategy(priority_of_patches,doses,sim_params,n_sims
	,sizes_function=expected_infected,time_rescaled=False):
	Sinit = sim_params[3].copy()
	n = Sinit.shape[0]
	R0 = sim_params[1].copy()
	if sizes_function==infected_sizes:
		outbreak_sizes = np.zeros(n_sims*n)
	else:
		outbreak_sizes = np.array([])
	for start_patch in range(n):
		initialState = np.zeros((n,1))
		initialState[start_patch] = 1
		if priority_of_patches==pro_rata_priority:
			Sinit_adj = Sinit*(1-initialState)
			allocation = Sinit_adj*doses/np.sum(Sinit_adj)
		else:
			if len(R0)==1:
				R0 = R0*np.ones((n,1))
			herd_threshold = np.ceil(Sinit*(1-1/R0))

			priority_params = (sim_params[0],sim_params[1],sim_params[2],sim_params[3]
				,initialState)
			if callable(priority_of_patches):
				priority_order = priority_of_patches(priority_params,sizes_function
					,n_sims,time_rescaled=time_rescaled)
			else:
				priority_order = priority_of_patches

			allocation = np.zeros((n,1))
			doses_temp = doses
			for patch in priority_order:
				if initialState[patch]==1:
					continue
				elif doses_temp >= herd_threshold[patch]:
					allocation[patch] = herd_threshold[patch]
					doses_temp -= herd_threshold[patch]
				elif doses_temp < herd_threshold[patch]:
					allocation[patch] = doses_temp
					doses_temp = 0

		strategy_R0 = R0*(Sinit-allocation)/Sinit
		strategy_Sinit = Sinit-allocation
		strategy_params = (sim_params[0],strategy_R0,
			sim_params[2],strategy_Sinit,initialState)
		sizes = sizes_function(
				strategy_params,n_sims,time_rescaled=time_rescaled)
		if sizes_function==infected_sizes:
			outbreak_sizes[start_patch*n_sims:(start_patch+1)*n_sims] = sizes
		else:
			outbreak_sizes = np.concatenate((outbreak_sizes,sizes))

	mean_infections = np.mean(outbreak_sizes)
	median_infections = np.median(outbreak_sizes)
	upper_quartile = np.percentile(outbreak_sizes,75)
	lower_quartile = np.percentile(outbreak_sizes,25)

	return(mean_infections,median_infections,upper_quartile,lower_quartile)

def size_given_NPI_strategy(priority_of_patches,doses,effectiveness,sim_params,n_sims,time_rescaled=False):
	Sinit = sim_params[3].copy()
	n = Sinit.shape[0]
	R0 = sim_params[1].copy()
	if len(R0)==1:
		R0 = R0*np.ones((n,1))
	initialState = sim_params[4].copy()

	if callable(priority_of_patches):
		priority_order = priority_of_patches(sim_params,n_sims,time_rescaled=time_rescaled)
	else:
		priority_order = priority_of_patches

	allocation = np.zeros((n,1))
	doses_temp = doses
	for patch in priority_order:
		if initialState[patch]==1:
			continue
		elif doses_temp >= Sinit[patch]:
			allocation[patch] = 1
			doses_temp -= Sinit[patch]

	strategy_R0 = R0*(1-effectiveness*allocation)
	strategy_params = (sim_params[0],strategy_R0,
		sim_params[2],sim_params[3],sim_params[4])
	size = infections_function(strategy_params,n_sims,time_rescaled=time_rescaled)
	
	return(allocation,size)

def random_priority(sim_params,infections_function,n_sims,time_rescaled=False):
	n = sim_params[3].shape[0]
	priority = np.arange(n)
	np.random.shuffle(priority)
	return(priority)

def risk_priority(sim_params,infections_function,n_sims,time_rescaled=False):
	if (infections_function==expected_infected) or\
	 (infections_function==infected_summary) or (infections_function==infected_sizes):
		probability = avg_importation(sim_params,n_sims,time_rescaled=time_rescaled)
	else:
		probability = importation_probabilities(sim_params,n_sims)
	return(np.flip(np.argsort(probability[:,0])))

def betweenness_priority(sim_params,infections_function,n_sims,time_rescaled=False):
	travel_rate = sim_params[0]
	with np.errstate(divide='ignore'):
		graph = nx.from_numpy_matrix(1/travel_rate,create_using=nx.DiGraph)
	betweenness = np.array(list(
		(nx.betweenness_centrality(graph,weight='weight')).values()))
	return(np.flip(np.argsort(betweenness)))

# Random-walk percolation centrality
def rwpc_priority(sim_params,infections_function,n_sims,time_rescaled=False):
	n = (sim_params[3].shape)[0]
	initialState = sim_params[4]
	source = np.arange(n)[np.where(initialState==1)[0]]
	centrality = rwpc(sim_params[0],source)
	return(np.flip(np.argsort(centrality)))

def degree_priority(sim_params,infections_function,n_sims,time_rescaled=False):
	travel_rate = sim_params[0]
	with np.errstate(divide='ignore'):
		graph = nx.from_numpy_matrix(travel_rate,create_using=nx.DiGraph)
	weighted_degree = np.array([node_degree[1]
		 for node_degree in list((nx.degree(graph,weight='weight')))])
	return(np.flip(np.argsort(weighted_degree)))

def out_degree_priority(sim_params,infections_function,n_sims,time_rescaled=False):
	travel_rate = sim_params[0]
	with np.errstate(divide='ignore'):
		graph = nx.from_numpy_matrix(travel_rate,create_using=nx.DiGraph)
	weighted_out_degree = np.array([node_degree[1]
		 for node_degree in list((graph.out_degree(weight='weight')))])
	return(np.flip(np.argsort(weighted_out_degree)))

def in_degree_priority(sim_params,infections_function,n_sims,time_rescaled=False):
	travel_rate = sim_params[0]
	with np.errstate(divide='ignore'):
		graph = nx.from_numpy_matrix(travel_rate,create_using=nx.DiGraph)
	weighted_in_degree = np.array([node_degree[1]
		 for node_degree in list((graph.in_degree(weight='weight')))])
	return(np.flip(np.argsort(weighted_in_degree)))

def in_closeness_priority(sim_params,infections_function,n_sims,time_rescaled=False):
	travel_rate = sim_params[0]
	with np.errstate(divide='ignore'):
		graph = nx.from_numpy_matrix(1/travel_rate,create_using=nx.DiGraph)
	closeness = np.array(list(
		(nx.closeness_centrality(graph,distance='weight')).values()))
	return(np.flip(np.argsort(closeness)))

def out_closeness_priority(sim_params,infections_function,n_sims,time_rescaled=False):
	travel_rate = sim_params[0]
	with np.errstate(divide='ignore'):
		graph = nx.from_numpy_matrix(1/travel_rate,create_using=nx.DiGraph)
	closeness = np.array(list(
		(nx.closeness_centrality(graph.reverse(),distance='weight')).values()))
	return(np.flip(np.argsort(closeness)))

def left_eigenvector_priority(sim_params,infections_function,n_sims,time_rescaled=False):
	travel_rate = sim_params[0]
	with np.errstate(divide='ignore'):
		graph = nx.from_numpy_matrix(travel_rate,create_using=nx.DiGraph)
	eigenvector = np.array(list(
		(nx.eigenvector_centrality_numpy(graph,weight='weight')).values()))
	return(np.flip(np.argsort(eigenvector)))

def right_eigenvector_priority(sim_params,infections_function,n_sims,time_rescaled=False):
	travel_rate = sim_params[0]
	with np.errstate(divide='ignore'):
		graph = nx.from_numpy_matrix(travel_rate,create_using=nx.DiGraph)
	eigenvector = np.array(list(
		(nx.eigenvector_centrality_numpy(graph.reverse(),weight='weight')).values()))
	return(np.flip(np.argsort(eigenvector)))

