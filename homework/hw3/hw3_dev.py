"""Template code for M3C 2016 Homework 3
"""
import numpy as np
import matplotlib.pyplot as plt
from n1 import network as net
def visualize(e):    
    """Visualize a network represented by the edge list
    contained in e
    """
    #Set an all-zeros array
    p=np.zeros((9,1))
    for i in e:
        p[i-1,0]=i
    #Find random values for x and y        
    x=np.random.random(len(p))
    y=np.random.random(len(p))
    plt.figure()
    #Plot the nodes
    for l in range(len(p)+1):
        if p[l-1]!=0:
            plt.plot(x[l-1],y[l-1],'o')
    #Plot the edjes for appropriate connections   
    for z in range(len(e)+1):
        a,b=e[z-1,0],e[z-1,1]
        c,d=(x[a-1],x[b-1]),(y[a-1],y[b-1])
        plt.plot(c,d,'-')
        plt.show()
    return
   
def degree(N0,L,Nt,display=False):
    """Compute and analyze degree distribution based on degree list, q,
    corresponding to a recursive network (with model parameters specified
    as input.
    """
    #Call our generate function
    net.generate(N0,L,Nt)
    #find the size of the array we will need
    o=max(net.qnet)
    k,p,d,u=np.zeros(o),np.zeros(o),np.zeros(o),np.zeros(o)
    #Loop to increase k for every connection appearing in qnet
    for l in net.qnet:
        k[l-1]=k[l-1]+1
    al=sum(k) 
    #Assign the probability 
    for i in range(o):
        u[i-1]=k[i-1]/al  
    t=0
    #Loop to assign all the non-zero numbers to d  and delete all the zero 
    #elements from p
    for j in range(o):
        if u[j]!=0:
            p[t]=u[j]
            d[t]=j
            t=t+1
    d=d[0:t]
    p=p[0:t]
    #If statement to plt or not plot the graph
    if display==True:
        #Try to find appropriate fit
        p2=p[0:int(len(p)/2)]
        d2=d[0:int(len(p)/2)]
        m=np.polyfit(d2,np.log(p2),1)
        width = 1/1.5
        #Plot the figure
        plt.figure()
        plt.bar(d, p, width,color="blue")
        plt.plot(d2,np.exp(d2*m[0]+m[1]),color="red")
        #Legend of the plot
        plt.legend(['D to P'])
        plt.xlabel('D')
        plt.ylabel('P')
        plt.title('Evgeniia Gleizer. Created by degree.')
        plt.show()     
    print d
    print p
    return
