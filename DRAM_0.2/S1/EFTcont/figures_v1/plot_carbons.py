import matplotlib.pyplot as plt
import numpy as np

def main():
	par = []
	f = open('DET.out')
	par = f.readlines()
	f.close()
	#print(len(par))
	det = [0] * len(par)
	for i in range(0, len(par)):
		par[i] = par[i].split()
		for j in range(2, len(par[i])):
			par[i][j] = float(par[i][j])
		det[i] = par[i][len(par[i])-1]
	par = []
	f = open('ZOO.out')
	par = f.readlines()
	f.close()
	#print(len(par))
	zoo = [0] * len(par)
	for i in range(0, len(par)):
		par[i] = par[i].split()
		for j in range(2, len(par[i])):
			par[i][j] = float(par[i][j])
		zoo[i] = par[i][len(par[i])-1] 
	par = []
        f = open('muN 1.out')
        par = f.readlines()
        f.close()
        #print(len(par))
        mu = [0] * len(par)
        for i in range(0, len(par)):
                par[i] = par[i].split()
                for j in range(2, len(par[i])):
                        par[i][j] = float(par[i][j])
                mu[i] = par[i][len(par[i])-1]
	par = []
	f = open('CHL_T.out')
	par = f.readlines()
	f.close()
	#print(len(par))
	means = [0] * len(par)
	chl = [0] * len(par)
	for i in range(0, len(par)):
		par[i] = par[i].split()
		for j in range(2, len(par[i])):
			par[i][j] = float(par[i][j])
		chl[i] = par[i][len(par[i])-1]
	par = []
	f = open('PHY 1.out')
	par = f.readlines()
	f.close()
	#print(len(par))
	means = [0] * len(par)
	phy = [0] * len(par)
	for i in range(0, len(par)):
		par[i] = par[i].split()
		for j in range(2, len(par[i])):
			par[i][j] = float(par[i][j])
		phy[i] = par[i][len(par[i])-1]
		#means[i] = float(sum(par[i][2:len(par[i])])/(len(par[i])-2))
	#print(means)
	fig1 = plt.figure()
	plt.subplot(2,1,1)
	#x = np.linspace(0,1080,num=len(phy))
	x = np.arange(1080)
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	plt.plot(x,phy[0:1080],x,zoo[0:1080],x,det[0:1080],x,chl[0:1080])
	plt.ylabel('Mean Surface Values ($mg m^{-3}$)')
	#plt.xlabel('Month')
	plt.legend(["Phytoplankton", "Zooplankton", "Detritus", "Chlorophyll"], prop={'size':8})
        month_vec = np.tile(['Jan', 'Apr', 'Jul', 'Oct'], 3)
        plt.xticks(np.arange(0, 1080, 90), month_vec)
	plt.title('Daily Averaged Mass Values at Surface')
	plt.axis([0, 1080, 0, max(phy)+0.1*max(phy)])
	#plt.show()
	
	chlp = [0] * len(phy)
	chlc = [0]*len(phy)
	for i in xrange(len(phy)):
		chlc[i] = chl[i]*1.0/(phy[i] + zoo[i] + det[i])
		chlp[i] = (chl[i]*1.0)/phy[i]
	par = []
        f = open('THE 1.out')
        par = f.readlines()
        f.close()
        #print(len(par))
        the = [0] * len(par)
        for i in range(0, len(par)):
                par[i] = par[i].split()
                for j in range(2, len(par[i])):
                        par[i][j] = float(par[i][j])
                the[i] = par[i][len(par[i])-1]
	#fig1 = plt.figure()
	plt.subplot(2,1,2)
        #x = np.linspace(0,1080,num=len(phy))
        x = np.arange(1080)
	plt.rc('text', usetex=True)
        plt.rc('font', family='serif')
        plt.plot(x,the[0:1080],x,chlc[0:1080],x,mu[0:1080])
        plt.ylabel('Chl/C (and Growth)')
        #plt.xlabel('Month')
        plt.legend(["Model Theta", "Chl : (Phy, Det, Zoo)", "Growth (muN)"], prop={'size':10})
        month_vec = np.tile(['Jan', 'Apr', 'Jul', 'Oct'], 3)
        plt.xticks(np.arange(0, 1080, 90), month_vec)
        plt.title('Daily Averaged Mass Values at Surface')
        plt.axis([0, 1080, 0, max(chlc)+0.1*max(chlc)])
	fig1.subplots_adjust(hspace=0.5)

        plt.show() 

	fig2 = plt.figure()
	x = np.arange(450)	
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	plot_params_chlp = chlp[300:750]
	plt.plot(x, plot_params_chlp, x, the[300:750])
	plt.ylabel('Surface Thetas (g/mol)')
	plt.xlabel('Time (Day)')
	plt.title('Daily Surface Theta Comparison')
	plt.legend(['Chl/Phy', 'Model Theta'])
        month_vec = np.tile(['Nov', 'Dec', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan'], 1)
        plt.xticks(np.arange(0, 450, 30), month_vec)
	plt.axis([0, len(plot_params_chlp), 0, max(plot_params_chlp)+0.1*(max(plot_params_chlp))])
	plt.show()
		
if __name__ == '__main__':
	main()


