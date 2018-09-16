# -*- coding: utf-8 -*-
"""Project, part 1"""
import numpy as np
import matplotlib.pyplot as plt
#import rw #Assumes rwmodule has been compiled with f2py to produce rw.so
from rw import rwmodule as rw


def analyze_rnet(Ntime,Nm,X0,N0,L,Nt,display):
	"""Input variables:
	Ntime: number of time steps
    	Nm: number of walks
    	X0: initial node, (node with maximum degree if X0=0)
    	N0,L,Nt: recursive network parameters
    	"""
    	#Call rwnet
    	X=rw.rwnet(Ntime,Nm,X0,N0,L,Nt,1)
    	X0=X[1,1]
    	#Create an empty F function and K array for maximum values
    	F=np.zeros(((Ntime+1),(N0+Nt)))
    	K=np.zeros(Ntime+1)

        #Loop over the time and then steps to fill in thr F(t,n)
    	for i in range(Ntime+1):
    	    for j in range(Nm):
    	        #Increasing the appropriate elements of F(t,n)
    	        a=X[i,j]-1
    	        F[i,a]=F[i,a]+1
            #Making F a probability matrix
    	    F[i,:]=F[i,:]/(1.0*Nm)
    	
    	#Loop over all times to get the greatest number of walkers at each step    
        for k in range(Ntime+1):
            K[k] = np.argmax(F[k,:])
            
        #Assiging x axis for the plot
    	x=np.arange(1,(Ntime+2))

        #If check for the display flag
    	if (display==True):
    	    #Plotting the figure
    	    plt.figure()
    	    plt.plot(x, K)
    	    #Legend of the plot
            plt.legend(['Maximum amount of walkers'])
            plt.xlabel('x range')
            plt.ylabel('F(t,n)')
            plt.title('Evgeniia Gleizer. Created by analyze_rnet.')
    	    plt.show()
        return F



def convergence_rnet(Ntime,Nm,X0,N0,L,Nt,display):
    	"""Input variables:
	Ntime: number of time steps
    	m: number of walks
    	X0: initial node, (node with maximum degree if X0=0)
    	N0,L,Nt: recursive network parameters
    	"""
    	#I am presenting here a graph of the average fraction of the walkers 
    	#Only of the nodes which fraction of the walkers is non-zero
    	#against m. It is easy to see that with Ntime=105 and m increasing,
    	#it tendes to converge to around 0
    	
    	#M=np.array([10,40,50,60,100])
    	l=len(Nm)
    	x=np.zeros(l)
    	y=np.zeros(l)
    	for j in range(l):
    	    F=analyze_rnet(105,Nm[j],0,5,2,200, False)
    	    #rw.rwnet(Ntime,Nm,X0,N0,L,Nt,1)
    	    #print sum(F[:,(X0-1)])
    	    for i in range (Ntime+1):
    	        if (F[i,(X0-1)]>0):
    	            x[j]=x[j]+1
    	    y[j]=sum(F[:,(X0-1)])/x[j]  
    	    
    	if (display==True):
    	    plt.figure()
    	    plt.ylim(0,0.2)
    	    plt.plot(Nm, y)
    	    #Legend of the plot
            plt.legend(['Average amount of walkers against Nm'])
            plt.xlabel('Nm')
            plt.ylabel('Av. walkers')
            plt.title('Evgeniia Gleizer. Created by convergence_rnet.')
    	    plt.show()  	
#This function helps us see that the maxes of F tend to converge to a
#relatively small number  	
    	
if __name__== '__main__':
#add code here to call functions and generate figures
     #Ntime,Nm,X0,N0,L,Nt,display)
    Ntime = 1e5
    Nm=1e3
    N0= 5
    L=2
    Nt=200
    analyze_rnet(Ntime,Nm,0,N0,L,Nt,True)
    convergence_rnet(Ntime,[10,20,30,40,100,1000],0,N0,L,Nt,True)