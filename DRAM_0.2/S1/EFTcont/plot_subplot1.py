import matplotlib.pyplot as plt
import numpy as np
def main():
	#Plotting temp subplot(1,1)
	temp = []
	f = open('Temp.out')
	temp = f.readlines()
	f.close()
	params_temp = [0] * len(temp)
	for i in range(0, len(temp)):
		temp[i] = temp[i].split()
		for j in range(2, len(temp[i])):
			temp[i][j] = float(temp[i][j])
		params_temp[i] = temp[i][len(temp[i])-1]
	fig1 = plt.figure()
	plt.subplot(311)
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	plot_params_temp = params_temp[300:750]
	plt.plot(plot_params_temp, 'b-')
	plt.ylabel('Surface Temperature (Celsius)')
	plt.xlabel('Time (Day)')
	plt.title('Daily Surface Temperature in Continuous Model')
        month_vec = np.tile(['N', 'D', 'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D', 'J'], 1)
        plt.xticks(np.arange(0, 450, 30), month_vec)
	plt.axis([0, len(plot_params_temp), 0, max(plot_params_temp)+0.1*(max(plot_params_temp))])

	#Plotting irradiation subplot(1,2)
	#si = []
	#f = open('SI  1.out')
	#si = f.readlines()
	#f.close()
	#params_si = [0] * len(si)
	#for i in range(0, len(si)):
#		si[i] = si[i].split()
#		for j in range(2, len(si[i])):
#			si[i][j] = float(si[i][j])
#		params_si[i] = si[i][len(si[i])-1]
#	plt.subplot(222)
#	plt.rc('text', usetex=True)
#	plt.rc('font', family='serif')
#	plot_params_si = params_si[300:750]
#	plt.plot(plot_params_si, 'b-')
#	plt.ylabel('Surface Irradiation')
#	plt.xlabel('Time (Day)')
#	plt.title('Daily Surface Irradiance in Continuous Model')
#        month_vec = np.tile(['N', 'D', 'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D', 'J'], 1)
#        plt.xticks(np.arange(0, 450, 30), month_vec)
#	plt.axis([0, len(plot_params_si), 0, max(plot_params_si)+0.1*(max(plot_params_si))])
	
	#Plotting PAR subplot(2,1)
	par = []
	f = open('PAR.out')
	par = f.readlines()
	f.close()
	params_par = [0] * len(par)
	for i in range(0, len(par)):
		par[i] = par[i].split()
		for j in range(2, len(par[i])):
			par[i][j] = float(par[i][j])
		params_par[i] = par[i][len(par[i])-1]
	plt.subplot(312)
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	plot_params_par = params_par[300:750]
	plt.plot(plot_params_par, 'b-')
	plt.ylabel('Surface PAR ($\mu m^{-2} s^{-1}$)')
	plt.xlabel('Time (Day)')
	plt.title('Daily Surface $PAR$ in Continuous Model')
        month_vec = np.tile(['N', 'D', 'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D', 'J'], 1)
        plt.xticks(np.arange(0, 450, 30), month_vec)
	plt.axis([0, len(plot_params_par), 0, max(plot_params_par)+0.1*(max(plot_params_par))])
	
	#Plotting NO3 subplot(2,2)
	no3 = []
	f = open('NO3.out')
	no3 = f.readlines()
	f.close()
	params_no3 = [0] * len(no3)
	for i in range(0, len(no3)):
		no3[i] = no3[i].split()
		for j in range(2, len(no3[i])):
			no3[i][j] = float(no3[i][j])
		params_no3[i] = no3[i][len(no3[i])-1]
	plt.subplot(313)
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	plot_params_no3 = params_no3[300:750]
	plt.plot(plot_params_no3, 'b-')
	plt.ylabel('Surface $NO_3$ ($mol m^{-3}$)')
	plt.xlabel('Time (Day)')
	plt.title('Daily Surface $NO_3$ in Continuous Model')
        month_vec = np.tile(['N', 'D', 'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D', 'J'], 1)
        plt.xticks(np.arange(0, 450, 30), month_vec)
	plt.axis([0, len(plot_params_no3), 0, max(plot_params_no3)+0.1*(max(plot_params_no3))])
	
	fig1.subplots_adjust(hspace=0.5, wspace=0.3)	
	plt.show()

if __name__ == '__main__':
	main()
