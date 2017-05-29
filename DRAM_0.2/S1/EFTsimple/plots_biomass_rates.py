import matplotlib.pyplot as plt
import numpy as np
def main():
	par = []
	f = open('D_DET.out')
	par = f.readlines()
	f.close()
	params = [0] * len(par)
	for i in range(0, len(par)):
		par[i] = par[i].split()
		for j in range(2, len(par[i])):
			par[i][j] = float(par[i][j])
		params[i] = par[i][len(par[i])-1]
	fig1 = plt.figure()
	plt.subplot(3,1,1)	
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	plot_params = params[300:750]
	plt.plot(plot_params, 'b-')
	plt.ylabel('$D/Ddays(DET)$')
	#plt.xlabel('Time (Day)')
	plt.title('Daily Surface Rate of Detritus Formation')
        month_vec = np.tile(['Nov', 'Dec', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan'], 1)
        plt.xticks(np.arange(0, 450, 30), month_vec)
	plt.axis([0, len(plot_params), 0, max(plot_params)+0.1*(max(plot_params))])
	
	par = []
	f = open('D_ZOO.out')
	par = f.readlines()
	f.close()
	params = [0] * len(par)
	for i in range(0, len(par)):
		par[i] = par[i].split()
		for j in range(2, len(par[i])):
			par[i][j] = float(par[i][j])
		params[i] = par[i][len(par[i])-1]
	plt.subplot(3,1,2)
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	plot_params = params[300:750]
	plt.plot(plot_params, 'b-')
	plt.ylabel('$D/Ddays(ZOO)$')
	#plt.xlabel('Time (Day)')
	plt.title('Daily Surface Rate of Zooplankton Formation')
        month_vec = np.tile(['Nov', 'Dec', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan'], 1)
        plt.xticks(np.arange(0, 450, 30), month_vec)
	plt.axis([0, len(plot_params), 0, max(plot_params)+0.1*(max(plot_params))])
	
	par = []
	f = open('D_P 1.out')
	par = f.readlines()
	f.close()
	params = [0] * len(par)
	for i in range(0, len(par)):
		par[i] = par[i].split()
		for j in range(2, len(par[i])):
			par[i][j] = float(par[i][j])
		params[i] = par[i][len(par[i])-1]
	plt.subplot(3,1,3)
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	plot_params = params[300:750]
	plt.plot(plot_params, 'b-')
	plt.ylabel('$D/Ddays(PHY)$')
	plt.xlabel('Time (Day)')
	plt.title('Daily Surface Rate of Phytoplankton Formation')
        month_vec = np.tile(['Nov', 'Dec', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan'], 1)
        plt.xticks(np.arange(0, 450, 30), month_vec)
	plt.axis([0, len(plot_params), 0, max(plot_params)+0.1*(max(plot_params))])
	
	fig1.subplots_adjust(hspace=0.8)
	plt.show()

if __name__ == '__main__':
	main()
