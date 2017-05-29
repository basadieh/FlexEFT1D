import matplotlib.pyplot as plt
import numpy as np
def main():
	par = []
	f = open('D_NO3.out')
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
	plot_params_d_no3 = params[300:750]
	plt.plot(plot_params_d_no3, 'b-')
	plt.ylabel('$D/Dt(NO3)$')
	plt.xlabel('Time (Day)')
	plt.title('Daily Surface Rate of Nitrate Formation in Continuous Model')
        month_vec = np.tile(['N', 'D', 'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D', 'J'], 1)
        plt.xticks(np.arange(0, 450, 30), month_vec)
	plt.axis([0, len(plot_params_d_no3), 0, max(plot_params_d_no3)+0.1*(max(plot_params_d_no3))])
	
	par = []
	f = open('QN  1.out')
	par = f.readlines()
	f.close()
	params = [0] * len(par)
	for i in range(0, len(par)):
		par[i] = par[i].split()
		for j in range(2, len(par[i])):
			par[i][j] = float(par[i][j])
		params[i] = par[i][len(par[i])-1]
	plt.subplot(2,2,2)
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	plot_params_qn = params[300:750]
	plt.plot(plot_params_qn, 'b-')
	plt.ylabel('Surface QN $(mol N/mol C)')
	plt.xlabel('Time (Day)')
	plt.title('Daily Surface Cell Quota ($Q_N$) in Continuous Model')
        month_vec = np.tile(['N', 'D', 'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D', 'J'], 1)
        plt.xticks(np.arange(0, 450, 30), month_vec)
	plt.axis([0, len(plot_params_qn), 0.07, max(plot_params_qn)+0.1*(max(plot_params_qn))])
	
	par = []
	f = open('D2N.out')
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
	plot_params_d2n = params[300:750]
	plt.plot(plot_params_d2n, 'b-')
	plt.ylabel('Surface D2N')
	plt.xlabel('Time (Day)')
	plt.title('Daily Surface D2N in Continuous Model')
        month_vec = np.tile(['N', 'D', 'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D', 'J'], 1)
        plt.xticks(np.arange(0, 450, 30), month_vec)
	plt.axis([0, len(plot_params_d2n), 0, max(plot_params_d2n)+0.1*(max(plot_params_d2n))])
	
	par = []
	f = open('CHL_T.out')
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
	plot_params_chl = params[300:750]
	plt.plot(plot_params_chl, 'b-')
	plt.ylabel('Surface CHL ($mol m^{-3}$)')
	plt.xlabel('Time (Day)')
	plt.title('Daily Surface CHL in Continuous Model')
        month_vec = np.tile(['N', 'D', 'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D', 'J'], 1)
        plt.xticks(np.arange(0, 450, 30), month_vec)
	plt.axis([0, len(plot_params_chl), 0, max(plot_params_chl)+0.1*(max(plot_params_chl))])

	fig1.subplots_adjust(hspace=0.4, wspace=0.3)
	plt.show()

if __name__ == '__main__':
	main()
