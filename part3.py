import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from p3 import sync as sy#assumes that fortran module sync has been compiled with f2py as p3.so



def oscillator(Nt,T,N,c,mu,s):
    """Simulate fully-coupled network of oscillators
    Compute simulation for Nt time steps between 0 and T (inclusive)
    input: N: number of oscillators
           c: coupling coefficient
           mu,sigma: distribution parameters for omega_i
    """
    #Assigning random normal w
    sy.w = np.random.normal(mu, s, N)
    omega=np.zeros(N)
    omega=sy.w
    #Assigning random uniform theta
    theta = np.random.uniform(0,2*np.pi,N)
    t0=0
    #Assigning dt and t 
    dt=(1.0*T)/Nt
    t=np.linspace(t0,T,dt)
    
    #Calling the rk4 function
    y=np.zeros((len(theta),Nt+1))
    order=np.zeros(Nt)
    y,order=sy.rk4(t0,theta,dt,Nt)

    theta0=theta[0]
    return t,omega,theta0,theta,order

if __name__ == '__main__':
    n,c,mu,s = 101,10.0,1.0,0.1
    Nt,T = 500,100
 
    t,omega,theta0,theta,order= oscillator(Nt,T,n,c,mu,s)
    print order
