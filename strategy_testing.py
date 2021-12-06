import sys
import numpy as np
import example_networks as net
import impact_of_strategy as ios
import sizes_and_coords as sc
import simplified_model as sm
import full_SIR as fs

network_name = sys.argv[1]
network_scale = float(sys.argv[2])/1e8
R0 = np.array([float(sys.argv[3])/100])
MU = np.array([float(sys.argv[4])/100])
vacc_doses = eval(sys.argv[5])
vacc_priority_function = eval("ios."+sys.argv[6]+"_priority")
n_sims = int(sys.argv[7])
simulation_function = eval(sys.argv[8])

if network_name=="ggrids": # one-to-one joined grids network
	Sinit = np.ones(18)*200000
	Sinit.shape = (18,1)
	coords = np.array([[0,0],[0,1],[0,2],[1,0],[1,1],[1,2],[2,0],[2,1],[2,2]])
	n,A = net.joined_network('gravity',network_scale,pops=Sinit,coords=coords,
		power=float(sys.argv[9])/100,join_strength=2e5*float(sys.argv[10])/100)
elif network_name=="cnorthernS": # Northwest England commuter network
	Sinit = sc.northern_pops
	n,A = net.network("commuter northern rate",network_scale)
	A = (A + A.transpose())/2
elif network_name=="air_us_urb_rate": # US air passenger traffic network
	Sinit = sc.top19_urb_pops
	n,A = net.network('air us top 19 urb rate',network_scale)
elif network_name=="four_corners": # four corners joined grids network
	Sinit = np.ones(18)*2e5
	Sinit.shape = (18,1)
	Sinit[0] = 4e5
	Sinit[9] = 4e5
	coords = np.array([[0,0],[0,1],[0,2],[1,0],[1,1],[1,2],[2,0],[2,1],[2,2]])
	join_indices = [[0,11],[2,9],[6,17],[8,15]]
	n,A = net.joined_network('gravity',network_scale,pops=Sinit,coords=coords,
		power=float(sys.argv[9])/100,join_indices=join_indices,
		join_strength=2e5*float(sys.argv[10])/400)
elif network_name=="square":
	Sinit = np.ones(4)*1e4
	Sinit.shape = (4,1)
	coords = np.array([[0,0],[0,1],[1,0],[1,1]])
	n,A = net.network('gravity',network_scale,pops=Sinit,coords=coords,
		power=float(sys.argv[9])/100)

sim_params = (A,R0,MU,Sinit)
summary = ios.summary_given_strategy(
	vacc_priority_function,vacc_doses,sim_params,n_sims
	,sizes_function=simulation_function,time_rescaled=True)
print(repr(summary))
