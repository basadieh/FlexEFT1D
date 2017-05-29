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
	
	max = -1
	maxIndex = 0
	plot_params = params[300:750]
	for i in range(len(plot_params)):
		if plot_params[i] > max:
			max = plot_params[i]
			maxIndex = i
	params2 = par[maxIndex][2:len(par[maxIndex])]
	depth_vec = par[0][2:len(par[0])]	
	for i in range(len(depth_vec)):
		depth_vec[i] = int(depth_vec[i])	
	fig1 = plt.figure()
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	plt.plot(depth_vec, params2, 'g-')
	plt.ylabel('Temperature')
	plt.xlabel('Depth (meters)')
	plt.title('Temperature vs Depth at Peak')
	
        #plt.xticks(np.arange(len(depth_vec)), depth_vec)
	#plt.axis([0, len(plot_params), 0, max(plot_params)+0.1*(max(plot_params))])
	plt.show()


if __name__ == '__main__':
	main()
