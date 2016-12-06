#Evgeniia Gleizer, CID: 00948999
import numpy as np
import matplotlib.pyplot as plt
from ns import netstats as ns #assumes netstats module has been compiled as ns.so with f2py

def convergence(n0,l,nt,m,numthreads,display=False):
	"""test convergence of qmax_ave"""
	#Calling stats, assigning the output to the variables
	qnetm,qmaxm,qvarm=ns.stats(n0,l,nt,m)
	#Assigning qmax_ave array
	qmax_ave=np.zeros(m)
	#Loop to calculate qmax_ave
	for j in range(1,m+1):
	    for i in range(j):
	        qmax_ave[j-1]=qmax_ave[j-1] + max(qnetm[:,i])
	    qmax_ave[j-1]=qmax_ave[j-1]/j
	#Assign x for our plot    
	x=np.arange(1,m+1)
	#Polifit to do the estimate 
	k,l=np.polyfit(np.log(x),np.log(qmax_ave),1)
	#If loop for the display flag
	if display==True:
	    plt.figure()
	    #Loglog plot of qmax_ave 
	    plt.loglog(x,qmax_ave,color="red")
	    #Loglog plot of the estimate
	    plt.loglog(x,np.exp(l)*x**k)
            #Legend of the plot
            plt.legend(['Log-log plot'])
            plt.xlabel('M')
            plt.ylabel('Qmax_ave')
            plt.title('Evgeniia Gleizer. Created by convergence.')
            plt.show()    
	
	
def speedup(n0,l,ntarray,marray):
	"""test speedup of stats_par"""
	#Finding length of the given arrays
	m=len(marray)
	nt=len(ntarray)
	#walltime=np.zeros(m)
	#walltimeomp=np.zeros(m)
	k=np.zeros((m,nt))
	#Loop to find times for all m and nt
	for i in range(m):
	    for j in range(nt):
	        q1=ntarray[j]
	        q2=marray[i]
	        walltime=ns.test_stats_omp(n0,l,q1,q2,1)
	        walltimeomp=ns.test_stats_omp(n0,l,q1,q2,2)
	        k[i,j]=(1.0*walltime/walltimeomp)
	        
	plt.figure()
	#Loop for our plot
	for g in range(nt):
	    plt.plot(marray,k[:,g])
	plt.show()
	#Legend of the plot
        plt.xlabel('M')
        plt.ylabel('time ratio')
        plt.title('Evgeniia Gleizer. Created by speedup.')

#The plotts indicate that the higher m is, the bigger gain we get from parralelisation
#It also shows the growth of the gain is not infinite and not directly
#proportional to m (which we saw in lectures)
def degree_average(n0,l,nt,m,display=False):
	"""compute ensemble-averaged degree distribution"""
	#Calling stats and assigning its output to varuables
	qnetm,qmaxm,qvarm=ns.stats(n0,l,nt,m)
	P=np.zeros(n0+nt)
	k=np.amax(qnetm)
	bins=np.zeros(k)
	#For loop to find the probabilities	
	for l in range(m):
            bins = np.bincount(qnetm[:,l])
            z=bins.size
            P[0:z]=P[0:z]+(1.0*bins)/len(qnetm[:,l])
            bins=np.zeros(k)
        #Averaging our P
        P=P/m
        x=np.zeros(len(P)+1)
        t=0
        #Loop to get rid of zeros
        for i in range(len(P)):
            if P[i]>0:
                x[t]=i
                t=t+1
        print len(P)
        x = x[x>0]      
        P = P[P>0]
        print len(x)
        #Polifit for logarithm of probability 
        k,l = np.polyfit(np.log(x),np.log(P),1)
        print k
        #Display flag check
        if display==True:
            
 	    plt.figure()
	    #Loglog plot of probability
	    plt.loglog(x,P,'o',color="red")
	    #Loglog plot of the estimate
	    plt.loglog(x,np.exp(l)*x**k)
            #Legend of the plot
            plt.legend(['Log-log plot'])
            plt.xlabel('X')
            plt.ylabel('P')
            plt.title('Evgeniia Gleizer. Created by degree_average.')
            plt.show()              
        return k

if __name__ == '__main__':
#Calling all the functions
    convergence(5,2,400,1000,1,display=True)
    degree_average(5,2,400,1000,display=True)
    u=np.array([2,3,4])
    o=np.array([3,4,5])
    speedup(5,2,np.array([10,100,200]),np.arange(300))