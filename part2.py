"""Project, part 2"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from n1 import network as net
from p1 import flunet as fn
from time import clock,time

def initialize(N0,L,Nt,pflag):
    """Generate network, and initial conditions
    If pflag is true, construct and return transport matrix
    """
    #Call the generate function 
    qmax,qnet,enet=net.generate(N0,L,Nt)
    #assigning the length of qnet
    l=len(qnet)
    q=0
    #Loop for finding maximum valur of qnet
    for i in range(l):
        if (qnet[i]>q):
            q=qnet[i]
            n=i
    #Setting the initial conditions       
    y0=np.zeros(3*l)
    y0[0:l]=1
    y0[n-1]=0.1
    y0[l+(n-1)]=0.05
    y0[2*l+(n-1)]=0.05
    
    #Checking the display flag condition
    if pflag==True:
        #Calling the adjacency matrix fir anet
        anet=net.adjacency_matrix(N0+Nt,enet)
        #Setting the size of the transport matrix
        P=np.zeros((len(anet),len(anet)))
        #Loop over it's length
        for j in range(len(anet)):
            summ=0
            #Loop for the denominatoe
            for m in range(len(anet)):
                #Computing the denominator
                summ=summ+(qnet[m]*anet[m,j])
            #Loop over the width of P
            for i in range(len(anet)):
                #Assigning P
                P[i,j]=float(qnet[i]*anet[i,j])/summ
        return y0,n, P
    else:
        return y0,n


    
def solveFluNet(P,T,Ntime,a,b0,b1,g,k,w,y0,option):
    """Simulate network flu model at Ntime equispaced times
    between 0 and T (inclusive). Add additional input variables
    as needed clearly commenting what they are and how they are  
    used
    """
    #Setting the size of dy
    s=len(P[0,:])
    dy=range(s)

    #add input variables to RHS functions if needed
    def RHSnet(y,t,a,b0,b1,g,k,w):   
        """RHS used by odeint to solve Flu model"""
        #Setting the initial values
        S=y[0:s]
        E=y[s:2*s]
        C=y[2*s:3*s] 
        
        #Creating the ODE, using vectors
        b = b0 + b1*(1.0+np.cos(2.0*np.pi*t))
        dS = k*(1-S) - b*C*S + w*(np.dot(P,S)-S)
        dE = b*C*S - (k+a)*E + w*(np.dot(P,E)-E)
        dC = a*E - (g+k)*C + w*(np.dot(P,C)-C)
        #Assigning the output to one array
        dy[0:s]= dS
        dy[s:2*s]=dE
        dy[2*s:3*s]=dC
        return dy

    def RHSnetF(y,t,a,b0,b1,g,k,w):
        """RHS used by odeint to solve Flu model"
        Calculations carried out by fn.rhs
        """
        #Calling the rhs from Fortran file
        dy=fn.rhs(y,t,a,b0,b1,g,k,w,P)
        return dy

    def RHSnetFomp(y,t,a,b0,b1,g,k,w):
        """RHS used by odeint to solve Flu model
        Calculations carried out by fn.rhs_omp
        """
        #Calling parallised rhs from Fortran file
        dy=fn.rhs_omp(y,t,a,b0,b1,g,k,w,P,2)
        return dy
        
    #Assigning timespan
    t=np.linspace(0,T,Ntime)

    #Creating all options for compiling
    if (option==1):

        Y=odeint(RHSnet,y0,t,args=(a,b0,b1,g,k,w))
    if (option==2):
        Y=odeint(RHSnetF,y0,t,args=(a,b0,b1,g,k,w))

    if (option==3):
        Y=odeint(RHSnetFomp,y0,t,args=(a,b0,b1,g,k,w))
        
    #Assigning solution to three arrays
    Ssol=np.zeros((Ntime,s))
    Esol=np.zeros((Ntime,s))
    Csol=np.zeros((Ntime,s))
    
    Ssol=Y[:,0:(s)]
    Esol=Y[:,(s):(2*s)]
    Csol=Y[:,(2*s):(3*s)]
    return t, Ssol, Esol, Csol


def analyze(N0,L,Nt,T,Ntime,a,b0,b1,g,k,threshold,warray,display=False):
    """analyze influence of omega on: 
    1. maximum contagious fraction across all nodes, Cmax
    2. time to maximum, Tmax
    3. maximum fraction of nodes with C > threshold, Nmax    
    Input: N0,L,Nt: recursive network parameters 
           T,Ntime,a,b0,b1,g,k: input for solveFluNet
           threshold: use to compute  Nmax
           warray: contains values of omega
    """
    #Assigning the arrays for the outputs
    s=len(warray)
    Cmax=np.zeros(s)
    Tmax=np.zeros(s)
    Csol2=np.zeros(Ntime)
    Nmaxa=np.zeros(Ntime)
    Nmax=np.zeros(s)
    #Loop to go over all w's
    for i in range(s):
        #Calling initialise and solveFluNet
        y0,n,P=initialize(N0,L,Nt,True)
        t1, Ssol1, Esol1, Csol1=solveFluNet(P,T,Ntime,a,b0,b1,g,k,warray[i],y0,1)
        #Assigning Cmax
        Cmax[i]=np.amax(Csol1)
        #Assigning Tmax
        Tmax[i]=t1[np.argmax(np.amax(Csol1,1))]
        #Loop to assign Nmax
        for j in range(Ntime):
            Csol2=Csol1[j,:]
            Nmaxa[j]=float(len(Csol2[Csol2>threshold]))/(N0+Nt)
        Nmax[i]=np.amax(Nmaxa)
        
    #Check for display flag      
    if (display==True):
        #Plotting the figure
    	plt.figure()
    	plt.plot(warray, Cmax)
    	#Legend of the plot
        plt.legend(['Cmax against warray'])
        plt.xlabel('w')
        plt.ylabel('Cmax')
        plt.title('Evgeniia Gleizer. Created by analyze.N0=10')
    	plt.show() 
    	   #New plot, for Tmax
    	plt.figure()
    	plt.plot(warray, Tmax)
    	#Legend of the plot
        plt.legend(['Tmax against warray'])
        plt.xlabel('w')
        plt.ylabel('Tmax')
        plt.title('Evgeniia Gleizer. Created by analyze.N0=10')
    	plt.show() 
    	  #New plot, for Nmax
    	plt.figure()
    	plt.plot(warray, Nmax)
    	#Legend of the plot
        plt.legend(['Nmax against warray'])
        plt.xlabel('w')
        plt.ylabel('Nmax')
        plt.title('Evgeniia Gleizer. Created by analyze.N0=10')
    	plt.show() 
    return Cmax,Tmax,Nmax
    

def visualize(enet,C,threshold):
    """Optional, not for assessment: Create figure showing nodes with C > threshold.
    Use crude network visualization algorithm from homework 3
    to display nodes. Contagious nodes should be red circles, all
    others should be black"""
    
    return None


def performance(P,T,Ntime,a,b0,b1,g,k,w,y0):
    
    """function to analyze performance of python, fortran, and fortran+omp approaches
        Add input variables as needed, add comment clearly describing the functionality
        of this function including figures that are generated and trends shown in figures
        that you have submitted
    """
    #Assigning all the arrays we will use
    s1=np.zeros(len(Ntime))
    s2=np.zeros(len(Ntime))
    s3=np.zeros(len(Ntime))
    e1=np.zeros(len(Ntime))
    e2=np.zeros(len(Ntime))
    e3=np.zeros(len(Ntime))
    t1=np.zeros(len(Ntime))
    t2=np.zeros(len(Ntime))
    t3=np.zeros(len(Ntime))
    sp1=np.zeros(len(Ntime))
    sp2=np.zeros(len(Ntime))
    
    #Loop over time
    for i in range(len(Ntime)):
        #Timing all the three versions of solveFluNet
        s1[i]=time()
        solveFluNet(P,T,Ntime[i],a,b0,b1,g,k,w,y0,1)
        e1[i]=time()

        s2[i]=time()
        solveFluNet(P,T,Ntime[i],a,b0,b1,g,k,w,y0,2)
        e2[i]=time()
      
        s3[i]=time()
        solveFluNet(P,T,Ntime[i],a,b0,b1,g,k,w,y0,3)
        e3[i]=time()
        #Finding the time
        t1[i]=e1[i]-s1[i]
        t2[i]=e2[i]-s2[i]
        t3[i]=e3[i]-s3[i]
        #Finding the speed ups
        sp1=t1/t2
        sp2=t2/t3
    #Plotting the speed ups
    plt.figure()
    line1, = plt.plot(Ntime,sp1, label='Python to Fortran code')
    line2, = plt.plot(Ntime,sp2, label='Fortran code with parallelisation to without one')
    #Legend of the plot
    plt.legend(handles=[line1,line2])
    plt.xlabel('Ntime')
    plt.ylabel('Relations of different RHS')
    plt.ylim(0,5)
    plt.title('Evgeniia Gleizer. Created by performance')
    plt.show()
    
    #Plotting the speeds
    plt.figure()
    line3, = plt.plot(Ntime,t1, label='Python code')
    line4, = plt.plot(Ntime,t2, label='Fortran code with parallelisation')
    line5, = plt.plot(Ntime,t3, label='Fortran code without parallelisation')
    #Legend of the plot
    plt.legend(handles=[line3,line4,line5])
    plt.xlabel('Ntime')
    plt.ylabel('Speeds')
    plt.title('Evgeniia Gleizer. Created by performance')
    plt.show()
    
    
    #It can be seen from the graphs that the first speed up is around two times which 
    #corresponds to existence of two cores. We also can see the changes after 
    #adding parallelisation and conclude that withing time the speed up tends to
    #be approximately the same and does not change that much.
    #The speed though goes down with time which is naturally what we would expect in the beginning
    #when going from small to bigger values.
if __name__ == '__main__':            
   a,b0,b1,g,k,w = 45.6,750.0,0.5,73.0,1.0,0.1
   warray=[0,1e-2,1e-1,0.2,0.5,1.0]
   N0 = 5 
   L=2
   Nt=500
   T=2, 
   threshold=0.1
   y0,n,P=initialize(5,3,3,True)
   Ntime=[5,10,20,30,40,50,100,200]
   T=2
   performance(P,T,Ntime,a,b0,b1,g,k,w,y0)