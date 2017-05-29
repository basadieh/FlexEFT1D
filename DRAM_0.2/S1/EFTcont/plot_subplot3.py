import matplotlib.pyplot as plt
import numpy as np
def main():
	par = []
	f = open('muN 1.out')
	par = f.readlines()
	f.close()
	params = [0] * len(par)
	for i in range(0, len(par)):
		par[i] = par[i].split()
		for j in range(2, len(par[i])):
			par[i][j] = float(par[i][j])
		params[i] = par[i][len(par[i])-1]
	fig1 = plt.figure()
	plt.subplot(2,1,1)
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	plot_params_mu = params[300:750]
	plt.plot(plot_params_mu, 'b-')
	plt.ylabel('Surface Net Growth Rate muNet ($day^{-1}$)')
	plt.xlabel('Time (Day)')
	plt.title('Daily Net Growth Rate at Surface in Continuous Model')
        month_vec = np.tile(['N', 'D', 'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D', 'J'], 1)
        plt.xticks(np.arange(0, 450, 30), month_vec)
	plt.axis([0, len(plot_params_mu), 0, max(plot_params_mu)+0.1*(max(plot_params_mu))])
	
	par = []
	f = open('Gra 1.out')
	par = f.readlines()
	f.close()
	params = [0] * len(par)
	for i in range(0, len(par)):
		par[i] = par[i].split()
		for j in range(2, len(par[i])):
			par[i][j] = float(par[i][j])
		params[i] = par[i][len(par[i])-1]
	plt.subplot(2,1,2)
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	plot_params_gra = params[300:750]
	plt.plot(plot_params_gra, 'b-')
	plt.ylabel('Surface Grazing Rate (1/day)')
	plt.xlabel('Time (Day)')
	plt.title('Daily Surface Grazing Rate in Continuous Model')
        month_vec = np.tile(['N', 'D', 'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D', 'J'], 1)
        plt.xticks(np.arange(0, 450, 30), month_vec)
	plt.axis([0, len(plot_params_gra), 0, max(plot_params_gra)+0.1*(max(plot_params_gra))])
	
	fig1.subplots_adjust(hspace=0.4)
	plt.show()

if __name__ == '__main__':
	main()
