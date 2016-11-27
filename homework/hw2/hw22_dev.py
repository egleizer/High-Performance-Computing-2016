"""Homework 2, part 2"""

"""Evgeniia Gleizer 00948999"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

def solveFlu(T,NT,a,b0,b1,g,k,y0):
    #giving the parameters for the differential equation
    def RHS(y,t):
        #initial values
        s=y[0]
        e=y[1]
        c=y[2]
        r=y[3]
        
        #equations
        b=b0+b1*(1+np.cos(2*np.pi*t))
        ds = k*(1-s)-b*c*s
        de = b*c*s-(k+a)*e
        dc = a*e-(g+k)*c
        dr = g*c-k*r
    
        dy = ds,de,dc,dr
        return dy
        
    #timepsan
    t=np.linspace(0,T,NT)
        
    #integrating ODE specified by RHS
    Y = odeint(RHS,y0,t)
    return Y[:,0],Y[:,1],Y[:,2],Y[:,3],t


def displaySolutions(t,y):
    #plotting the first figure
    plt.figure()
    #obtaining only such t's that t>1
    t1=t[t>1]
    #as C increases, we take elements starting from (len t - len t1)
    d=len(t) - len(t1)
    plt.plot(t1,Y[2,d:])
    plt.legend(['C to t'])
    plt.xlabel('t')
    plt.ylabel('C')
    plt.title('Evgeniia Gleizer. Created by displaySolution.')
    plt.show()

    #plotting the second figure
    plt.figure()
    #as we already now the gap from the left, we ammend C and E straight away
    plt.plot(Y[2,d:],Y[1,d:])
    plt.legend(['E to C'])
    plt.xlabel('C')
    plt.ylabel('E')
    plt.title('Evgeniia Gleizer. Created by displaySolution.')
    plt.show()
    
    
def linearFit(b0=5):
    #input all the given data
    a,b1,k,g=1.0,0.2,0.1,0.2
    S,E,C,R,t=solveFlu(10,100,a,b0,b1,g,k,[0.997,0.001,0.001,0.001])
    #ammend C to only have C such that it's element are less then 0.1
    C1=C[C<0.1]
    #find the length of C1
    k=len(C1)
    #amment t according to the length of C1
    t1=t[0:k]
    #use polyfit to find the coefficients of the approximation polynomial
    m=np.polyfit(t1,np.log(C1),1)

    
    plt.figure() 
    #plotting the log function againt time
    plt.plot(t,np.log(C))
    #plotting the approximation polynomial
    plt.plot(t,t*m[0]+m[1])
    plt.legend(['C1 to t1'])
    plt.xlabel('C1')
    plt.ylabel('t1')
    plt.title('Evgeniia Gleizer. Created by linearFit. Comparison of solution and expected exponential growth')
    plt.show()
    print C1, t1,m
    return


def fluStats(t,y,i1,i2):
    #checking the initial conditions
    if i1<i2:
        #calculating appropriate yn
        yn=y[i1-1:i2-1,:]
        #obtaining varience
        v=np.array([np.var(yn[:,0]),np.var(yn[:,1]),np.var(yn[:,2]),np.var(yn[:,3])])
        #obtaining mean
        k=np.array([np.mean(yn[:,0]),np.mean(yn[:,1]),np.mean(yn[:,2]),np.mean(yn[:,3])])
    else:
        print 'error'
    return v,k
 
    
    
if __name__ == '__main__':
    S,E,C,R,t=solveFlu(5, 1000, 45.6,750.0,1000.0,73.0,1.0,[0.1,0.05,0.05,0.8])
    Y=np.array([S,E,C,R])
    displaySolutions(t,Y)
    Nt=len(t)
    fluStats(t,Y,int(Nt/2),Nt)
    