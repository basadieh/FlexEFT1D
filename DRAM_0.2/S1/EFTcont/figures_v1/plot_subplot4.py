import matplotlib.pyplot as plt
import numpy as np
def main():
	par = []
	f = open('Aks.out')
	par = f.readlines()
	f.close()
	params = [0] * len(par)
	for i in range(0, len(par)):
		par[i] = par[i].split()
		for j in range(2, len(par[i])):
			par[i][j] = float(par[i][j])
		params[i] = par[i][len(par[i])-1]
	fig1 = plt.figure()
	plt.subplot(2,2,1)
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	plot_params_aks = params[300:750]
	plt.plot(plot_params_aks, 'b-')
	plt.ylabel('Mean Aks')
	plt.xlabel('Time (Day)')
	plt.title('Daily Average Aks')
        month_vec = np.tile(['N', 'D', 'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D', 'J'], 1)
        plt.xticks(np.arange(0, 450, 30), month_vec)
	plt.axis([0, len(plot_params_aks), 0, max(plot_params_aks)+0.1*(max(plot_params_aks))])
	
	par = []
	f = open('LNO 1.out')
	par = f.readlines()
	f.close()
	params = [0] * len(par)
	for i in range(0, len(par)):
		par[i] = par[i].split()
		for j in range(2, len(par[i])):
			par[i][j] = float(par[i][j])
		params[i] = par[i][len(par[i])-1]
	plt.subplot(2,2,3)
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	plot_params_lno = params[300:750]
	plt.plot(plot_params_lno, 'b-')
	plt.ylabel('Surface LNO')
	plt.xlabel('Time (Day)')
	plt.title('Daily Surface Nutrient Limitation Index')
        month_vec = np.tile(['N', 'D', 'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D', 'J'], 1)
        plt.xticks(np.arange(0, 450, 30), month_vec)
	plt.axis([0, len(plot_params_lno), 0.3, max(plot_params_lno)+0.1*(max(plot_params_lno))])
	
	par = []
	f = open('NPP_T.out')
	par = f.readlines()
	f.close()
	params = [0] * len(par)
	for i in range(0, len(par)):
		par[i] = par[i].split()
		for j in range(2, len(par[i])):
			par[i][j] = float(par[i][j])
		params[i] = par[i][len(par[i])-1]
	plt.subplot(2,2,4)
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	plot_params_npp = params[300:750]
	plt.plot(plot_params_npp, 'b-')
	plt.ylabel('Surface $NPP_T$')
	plt.xlabel('Time (Day)')
	plt.title('Daily Surface $NPP_T$')
        month_vec = np.tile(['N', 'D', 'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D', 'J'], 1)
        plt.xticks(np.arange(0, 450, 30), month_vec)
	plt.axis([0, len(plot_params_npp), 0, max(plot_params_npp)+0.1*(max(plot_params_npp))])
	
	fig1.subplots_adjust(hspace=0.5, wspace=0.3)
	plt.show()

if __name__ == '__main__':
	main()
