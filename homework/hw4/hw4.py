"""template code for hw4
"""
import numpy as np
import matplotlib.pyplot as plt
from ns import netstats as ns #assumes netstats module has been compiled as ns.so with f2py

def convergence(n0,l,nt,m,numthreads,display=False):
	"""test convergence of qmax_ave"""
	
	qnetm,qmaxm,qvarm=ns.stats(n0,l,nt,m)
	qmax_ave=np.zeros(m)
	for j in range(1,m+1):
	    for i in range(j):
	        qmax_ave[j-1]=qmax_ave[j-1] + max(qnetm[:,i])
	    qmax_ave[j-1]=qmax_ave[j-1]/j
	
	if display==True:
	    #l=np.polyfit(np.arange(1,m+1),np.log(qmax_ave),2)
	    plt.figure()
	    x=np.log(np.arange(1,m+1))
	    #plt.plot(x,np.exp(x*l[0]+l[1]),color="yellow")
	    plt.plot(x,np.log(qmax_ave),color="red")
            #Legend of the plot
            plt.legend(['Log-log plot'])
            plt.xlabel('M')
            plt.ylabel('Qmax_ave')
            plt.title('Evgeniia Gleizer. Created by convergence.')
            plt.show()    

	
	
def speedup(n0,l,ntarray,marray):
	"""test speedup of stats_par"""
	


def degree_average(n0,l,nt,m,display=False):
	"""compute ensemble-averaged degree distribution"""	




#if __name__ == '__main__'
#Add code to here to call functions above
