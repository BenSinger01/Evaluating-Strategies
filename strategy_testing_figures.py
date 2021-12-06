import numpy as np
import csv
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib import cm
from scipy import stats
import matplotlib

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Palatino"],
})
plt.rc('text.latex', preamble=
	r'\usepackage[T1]{fontenc}\usepackage{palatino}\usepackage{amsmath}')
matplotlib.verbose.level = 'debug-annoying'

font = {'size'   : 10}

matplotlib.rc('font', **font)

def scale_line_plot(inname,ax,strategy,R0,scale,color,dose_low=0,dose_high=1e16
	,centre='mean',total_scale=1):
	mean_array = np.array([])
	median_array = np.array([])
	upper_array = np.array([])
	lower_array = np.array([])
	doses = np.array([])

	with open(inname) as csvfile:
		readCSV = csv.reader(csvfile, delimiter=';')
		for row in readCSV:
			if row[3] == strategy and int(row[1])==int(R0*100) and\
			 float(row[0])==scale and float(row[2])>=dose_low and\
			  float(row[2])<=dose_high:
				summary = eval(row[4])
				mean_array = np.append(mean_array,summary[0])
				median_array = np.append(median_array,summary[1])
				upper_array = np.append(upper_array,summary[2])
				lower_array = np.append(lower_array,summary[3])
				doses = np.append(doses,float(row[2]))

	if strategy=="rwpc":
		strategy = "RWPC"

	mean_array /= total_scale
	median_array/=total_scale
	upper_array /= total_scale
	lower_array /= total_scale
	doses /= total_scale

	if centre=='mean':
		im=ax.plot(doses,mean_array,label=strategy,color=color)
	elif centre=='median':
		im=ax.plot(doses,median_array,label=strategy,color=color)
		im=ax.fill_between(doses,upper_array,lower_array,color=color,alpha=0.2)

	return(im)

def ggrids_line_plot(inname,ax,strategy,R0,power,scale,color,centre='mean'
	,total_scale=1):
	mean_array = np.array([])
	median_array = np.array([])
	upper_array = np.array([])
	lower_array = np.array([])
	doses = np.array([])

	with open(inname) as csvfile:
		readCSV = csv.reader(csvfile, delimiter=';')
		for row in readCSV:
			if row[4] == strategy and int(row[2])==int(R0*100) and\
			 float(row[1])==float(scale) and int(row[0])==int(power*100):
				summary = eval(row[5])
				mean_array = np.append(mean_array,summary[0])
				median_array = np.append(median_array,summary[1])
				upper_array = np.append(upper_array,summary[2])
				lower_array = np.append(lower_array,summary[3])
				doses = np.append(doses,int(row[3]))

	if strategy=="rwpc":
		strategy = "RWPC"

	mean_array /= total_scale
	median_array/=total_scale
	upper_array /= total_scale
	lower_array /= total_scale
	doses /= total_scale

	if centre=='mean':
		im=ax.plot(doses,mean_array,label=strategy,color=color)
	elif centre=='median':
		im=ax.plot(doses,median_array,label=strategy,color=color)
		im=ax.fill_between(doses,upper_array,lower_array,color=color,alpha=0.2)

	return(im)

def brute_force_plot(inname,ax,R0,power,scale,dose_low=0,dose_high=1e16
	,total_scale=1):
	risk_array = np.array([])
	doses = np.array([])

	with open(inname) as csvfile:
		readCSV = csv.reader(csvfile, delimiter=';')
		for row in readCSV:
			if int(row[2])==int(R0*100) and float(row[1])==float(scale)\
			 and int(row[0])==int(power*100):
				risks = [eval(array) for array in row[slice(4,-1,2)] if array[0]=='n']
				avg_risk = sum([np.sum(risk) for risk in risks])/len(risks)
				risk_array = np.append(risk_array,avg_risk)
				doses = np.append(doses,np.sum(eval(row[5])))

	risk_array /= total_scale

	for i in range(len(doses)):
		if doses[i]>=dose_low and doses[i]<=dose_high:
			im=ax.scatter(doses[i]/total_scale,risk_array[i],color='k',s=8,zorder=3)
			im=ax.text((doses[i]+1.75e4)/total_scale,risk_array[i],(str(i+1)),fontsize=8,
				ha='left',va='center',zorder=3)

	return(im)

def topology_line_plot(inname,ax,strategy,R0,scale,topology,color,centre='mean'
	,total_scale=1):
	mean_array = np.array([])
	median_array = np.array([])
	upper_array = np.array([])
	lower_array = np.array([])
	doses = np.array([])

	with open(inname) as csvfile:
		readCSV = csv.reader(csvfile, delimiter=';')
		for row in readCSV:
			if row[5] == strategy and int(row[3])==int(R0*100)\
			 and row[0]==topology and float(row[2])==float(scale):
				summary = eval(row[6])
				mean_array = np.append(mean_array,summary[0])
				median_array = np.append(median_array,summary[1])
				upper_array = np.append(upper_array,summary[2])
				lower_array = np.append(lower_array,summary[3])
				doses = np.append(doses,int(row[4]))

	if strategy=="rwpc":
		strategy = "RWPC"

	mean_array /= total_scale
	median_array/=total_scale
	upper_array /= total_scale
	lower_array /= total_scale
	doses /= total_scale

	if centre=='mean':
		im=ax.plot(doses,mean_array,label=strategy,color=color)
	elif centre=='median':
		im=ax.plot(doses,median_array,label=strategy,color=color)
		im=ax.fill_between(doses,upper_array,lower_array,color=color,alpha=0.2)

	return(im)

# # ggrids scale figure
# fig, axeses = plt.subplots(3,5,figsize=(8,4.2))
# fig.add_subplot(111,frameon=False)
# centre = 'mean'
# R0s = (1.2,1.5,2,4,8)
# scales = (1e-3,1e-2,1e-1)
# power = 2
# strats = ("risk","betweenness","rwpc")
# c = ['#648FFF','#DC267F','#FFB000','#785EF0','#FE6100']
# for s in range(len(scales)):
# 	axes = axeses[s]
# 	scale = scales[s]
# 	for i in range(len(R0s)):
# 		for j in range(len(strats)):
# 			image = ggrids_line_plot("strategy_testing_summary_ggrids.csv",
# 				axes[i],strats[j],R0s[i],power,scale,c[j],centre=centre,total_scale=1e5)
# 		image = brute_force_plot("brute_force_testing_ggrids_time_rescaled_with_doses.csv",
# 			axes[i],R0s[i],power,scale,dose_low=1e5,dose_high=5e5,total_scale=1e5)

# for i in range(len(R0s)):
# 	axeses[0,i].set_title(r"$R_0 =$"+str(R0s[i]),y=1.02)
# 	for j in range(2):
# 		axeses[j,i].set_xticks([])
# for i in range(len(scales)):
# 	axeses[i,len(R0s)-1].yaxis.set_label_position("right")
# 	axeses[i,len(R0s)-1].set_ylabel(
# 		r'$\lambda = 10^{'+ str(int(np.log10(scales[i]))-8)+r'}$')

# plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
# plt.grid(False)
# plt.ylabel(r"Mean total infections ($\times10^5$)",labelpad=10)
# plt.xlabel(r"Vaccine doses, $n_d$ ($\times10^5$)",labelpad=5)
# fig.tight_layout(pad=0)

# fig.subplots_adjust(top=0.85)
# handles, labels = axeses[0,0].get_legend_handles_labels()
# fig.legend(handles, labels,loc="upper center",ncol=len(strats))

# plt.savefig("ggrids_scale_figure.pdf",bbox_inches='tight')

# topology scale figure
fig, axeses = plt.subplots(3,5,figsize=(8,4.2))
fig.add_subplot(111,frameon=False)
centre = 'mean'
R0s = (1.2,1.5,2,4,8)
scales = (1e-3,1e-2,1e-1)
topology = "four_corners"
strats = ("risk","betweenness","rwpc")
c = ['#648FFF','#DC267F','#FFB000','#785EF0','#FE6100']
for s in range(len(scales)):
	axes = axeses[s]
	scale = scales[s]
	for i in range(len(R0s)):
		for j in range(len(strats)):
			image = topology_line_plot("strategy_testing_summary_ggrids_join_topology_tst.csv",
				axes[i],strats[j],R0s[i],scale,topology,c[j],centre=centre,total_scale=1e5)

for i in range(len(R0s)):
	axeses[0,i].set_title(r"$R_0 =$"+str(R0s[i]),y=1.02)
	for j in range(2):
		axeses[j,i].set_xticks([])
for i in range(len(scales)):
	axeses[i,len(R0s)-1].yaxis.set_label_position("right")
	axeses[i,len(R0s)-1].set_ylabel(
		r'$\lambda = 10^{'+ str(int(np.log10(scales[i]))-8)+r'}$')

plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
plt.grid(False)
plt.ylabel(r"Mean total infections ($\times10^5$)",labelpad=10)
plt.xlabel(r"Vaccine doses, $n_d$ ($\times10^5$)",labelpad=5)
fig.tight_layout(pad=0)

fig.subplots_adjust(top=0.85)
handles, labels = axeses[0,0].get_legend_handles_labels()
fig.legend(handles, labels,loc="upper center",ncol=len(strats))

plt.savefig("ggrids_topology_figure_four_corners_tst.pdf",bbox_inches='tight')

# air us scale figure
fig, axeses = plt.subplots(3,5,figsize=(8,4.2))
fig.add_subplot(111,frameon=False)
centre = 'mean'
R0s = (1.2,1.5,2,4,8)
scales = (1,10,100)
strats = ("risk","betweenness","rwpc")
c = ['#648FFF','#DC267F','#FFB000','#785EF0','#FE6100']
for s in range(len(scales)):
	axes = axeses[s]
	scale = scales[s]
	for i in range(len(R0s)):
		for j in range(len(strats)):
			image = scale_line_plot("strategy_testing_summary_air_us_urb_rate_tst.csv",
				axes[i],strats[j],R0s[i],scale,c[j],centre=centre,total_scale=1e6)

for i in range(len(R0s)):
	axeses[0,i].set_title(r"$R_0 =$"+str(R0s[i]),y=1.02)
	for j in range(2):
		axeses[j,i].set_xticks([])
for i in range(len(scales)):
	axeses[i,len(R0s)-1].yaxis.set_label_position("right")
	axeses[i,len(R0s)-1].set_ylabel(
		r'$s = 10^{'+ str(int(np.log10(scales[i]))-8)+r'}$')

plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
plt.grid(False)
plt.ylabel(r"Mean total infections ($\times10^6$)",labelpad=10)
plt.xlabel(r"Vaccine doses, $n_d$ ($\times10^6$)",labelpad=5)
fig.tight_layout(pad=0)

fig.subplots_adjust(top=0.85)
handles, labels = axeses[0,0].get_legend_handles_labels()
fig.legend(handles, labels,loc="upper center",ncol=len(strats))

plt.savefig("air_us_urb_rate_scale_figure_tst.pdf",bbox_inches='tight')

# cnorthern scale figure
fig, axeses = plt.subplots(3,5,figsize=(8,4.2))
fig.add_subplot(111,frameon=False)
centre = 'mean'
R0s = (1.2,1.5,2,4,8)
scales = (1e4,1e5,1e6)
strats = ("risk","betweenness","rwpc")
c = ['#648FFF','#DC267F','#FFB000','#785EF0','#FE6100']
for s in range(len(scales)):
	axes = axeses[s]
	scale = scales[s]
	for i in range(len(R0s)):
		for j in range(len(strats)):
			image = scale_line_plot("strategy_testing_summary_cnorthernS_tst.csv",
				axes[i],strats[j],R0s[i],scale,c[j],centre=centre,total_scale=1e5)

for i in range(len(R0s)):
	axeses[0,i].set_title(r"$R_0 =$"+str(R0s[i]),y=1.02)
	for j in range(2):
		axeses[j,i].set_xticks([])
for i in range(len(scales)):
	axeses[i,len(R0s)-1].yaxis.set_label_position("right")
	axeses[i,len(R0s)-1].set_ylabel(
		r'$s = 10^{'+ str(int(np.log10(scales[i]))-8)+r'}$')

plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
plt.grid(False)
plt.ylabel(r"Mean total infections ($\times10^5$)",labelpad=10)
plt.xlabel(r"Vaccine doses, $n_d$ ($\times10^5$)",labelpad=5)
fig.tight_layout(pad=0)

fig.subplots_adjust(top=0.85)
handles, labels = axeses[0,0].get_legend_handles_labels()
fig.legend(handles, labels,loc="upper center",ncol=len(strats))

plt.savefig("cnorthernS_scale_figure_tst.pdf",bbox_inches='tight')
