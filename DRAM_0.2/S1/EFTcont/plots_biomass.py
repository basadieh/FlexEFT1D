import matplotlib.pyplot as plt
import numpy as np
def main():
	par = []
	f = open('DET.out')
	par = f.readlines()
	f.close()
	params = [0] * len(par)
	for i in range(0, len(par)):
		par[i] = par[i].split()
		for j in range(2, len(par[i])):
			par[i][j] = float(par[i][j])
		params[i] = par[i][len(par[i])-1]
	fig1 = plt.figure()
	#plt.subplot(3,1,1)
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	plot_params_det = params[300:750]
	'''
	plt.plot(plot_params, 'b-')
	plt.ylabel('Surface Detritus')
	plt.xlabel('Time (Day)')
	plt.title('Daily Surface Detritus')
        month_vec = np.tile(['Nov', 'Dec', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan'], 1)
        plt.xticks(np.arange(0, 450, 30), month_vec)
	plt.axis([0, len(plot_params), 0, max(plot_params)+0.1*(max(plot_params))])
	'''
	par = []
	f = open('PHY 1.out')
	par = f.readlines()
	f.close()
	params = [0] * len(par)
	for i in range(0, len(par)):
		par[i] = par[i].split()
		for j in range(2, len(par[i])):
			par[i][j] = float(par[i][j])
		params[i] = par[i][len(par[i])-1]
	#plt.subplot(3,1,2)
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	plot_params_phy = params[300:750]
	'''
	plt.plot(plot_params, 'b-')
	plt.ylabel('Surface Phytoplankton ($mol m^{-3}$)')
	plt.xlabel('Time (Day)')
	plt.title('Daily Surface Phytoplankton')
        month_vec = np.tile(['Nov', 'Dec', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan'], 1)
        plt.xticks(np.arange(0, 450, 30), month_vec)
	plt.axis([0, len(plot_params), 0, max(plot_params)+0.1*(max(plot_params))])
	'''
	par = []
	f = open('ZOO.out')
	par = f.readlines()
	f.close()
	params = [0] * len(par)
	for i in range(0, len(par)):
		par[i] = par[i].split()
		for j in range(2, len(par[i])):
			par[i][j] = float(par[i][j])
		params[i] = par[i][len(par[i])-1]
	#plt.subplot(3,1,3)
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	plot_params_zoo = params[300:750]
	'''
	plt.plot(plot_params, 'b-')
	plt.ylabel('Surface Zooplankton')
	plt.xlabel('Time (Day)')
	plt.title('Daily Surface Zooplankton')
        month_vec = np.tile(['Nov', 'Dec', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan'], 1)
        plt.xticks(np.arange(0, 450, 30), month_vec)
	plt.axis([0, len(plot_params), 0, max(plot_params)+0.1*(max(plot_params))])
	'''
	x = np.arange(450)
	plt.plot(x, plot_params_det, x, plot_params_phy, x, plot_params_zoo)
	plt.ylabel('Surface Biomass ($mg m^{-3}$)')
	plt.xlabel('Time (days)')
	plt.title('Daily Surface Biomass in Continuous Model') 
        month_vec = np.tile(['Nov', 'Dec', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan'], 1)
	plt.xticks(np.arange(0, 450, 30), month_vec)
	plt.axis([0, len(plot_params_det), 0, max(plot_params_phy)+0.1*max(plot_params_phy)])
	plt.legend(['Detritus', 'Phytoplankton', 'Zooplankton'])
	plt.show()

if __name__ == '__main__':
		main()
