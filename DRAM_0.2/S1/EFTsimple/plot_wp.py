import matplotlib.pyplot as plt
import numpy as np
def main():
	par = []
	f = open('w_p 1.out')
	par = f.readlines()
	f.close()
	#print(len(par))
	means = [0] * len(par)
	params = [0] * len(par)
	for i in range(0, len(par)):
		par[i] = par[i].split()
		for j in range(2, len(par[i])):
			par[i][j] = float(par[i][j])
		params[i] = par[i][len(par[i])-1]
		#means[i] = float(sum(par[i][2:len(par[i])])/(len(par[i])-2))
	#print(means)
	fig1 = plt.figure()
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	plt.plot(params, 'bo')
	plt.ylabel('Surface $W_P$')
	plt.xlabel('Time (Day)')
        month_vec = np.tile(['Jan', 'Apr', 'Jul', 'Oct'], 3)
        plt.xticks(np.arange(0, 1080, 90), month_vec)
	plt.title('Daily Surface $W_P$')
	plt.axis([0, len(params), 0, max(params)])
	plt.show()

if __name__ == '__main__':
	main()


