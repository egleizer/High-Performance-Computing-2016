"""An illustrative example for integrating an ODE
Uses scipy's odeint function to solve d^2y/dt^2 = -w^2 y
After solving the ODE, the numerical and exact solutions
are compared and plotted
"""

#import modules
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

#function provides RHS to odeint
def RHS(y,t,w=1.0):

    y1 = y[0]
    y2 = y[1]
    
    dy1 = y2
    dy2 = -(w**2)*y1
    
    dy = dy1,dy2
    return dy
    

#timepsan
t=np.linspace(0,2*np.pi,1000)

#initial condition
y0 = 1.0,0.0

#set parameter, w
w = 2.0

#integrate ODE specified by RHS
Y = odeint(RHS,y0,t,args=(w,))

#exact solution and error
Yexact = y0[0]*np.cos(w*t)
error = np.abs(Yexact-Y[:,0].T) #the .T takes the transpose of Y
print "max error=",error.max()

#make figure showing solutions
plt.figure()
plt.plot(t,Y[:,0])
plt.plot(t,Yexact,'r--')
plt.legend(['computed','exact'])
plt.xlabel('t')
plt.ylabel('y')
plt.title('solution to d^2y/dt^2 = -w^2 y')
plt.xlim(0,t[-1])
plt.show()