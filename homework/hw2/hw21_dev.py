"""Homework 2, part 1"""
"""Evgeniia Gleizer 00948999"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import misc
if __name__ == '__main__':
    print""
def crop(F,T=0,B=0,L=0,R=0):
    F = misc.imread("image.jpg")
    #checking that T is valid
    if T>0:
        #calculating the gap from above
        t=int(F.shape[0]*T)-1
    #calculating the gap from below
    b=F.shape[0]-int(F.shape[0]*B)-1
    #checking that L is valid
    if L>0:
        #calculating the gap from left
        l=int(F.shape[1]*L)-1
    #calculating the gap from right
    r=F.shape[1]-int(F.shape[1]*R)-1
    F=F[t:b,l:r]
    #plot the new image
    plt.figure()
    plt.imshow(F)
    plt.title('Evgeniia Gleizer. Created by crop.')
    plt.show()
    return F

def smooth(F):
    F = misc.imread("image.jpg")
    #establishing the length and width of the original image
    b=F.shape[0]-1
    r=F.shape[1]-1
    #creating an array full of zeros to further fill it in with the nex data
    K=np.zeros(shape=(b+1,r+1,3))
    for n in range(3):
        #first we fill in all the corners
        K[0,0,n]=np.sum(F[0:2,0:2,n])/4.0
        K[b,0,n]=np.sum(F[b-1:b+1,0:2,n])/4.0
        K[b,r,n]=np.sum(F[b-1:b+1,r-1:r+1,n])/4.0
        K[0,r,n]=np.sum(F[0:2,r-1:r+1,n])/4.0
        #we create an array like F but with another layor of zeros on the border 
        #we do this to make calculations easier as otherwise we have to specify 
        #different regions that we are summing
        L=np.pad(F,(1,1), mode='constant')
        for x in range(b):
          for y in range(r):
           #now we are calculating the average for the outline border of our F
           if x==0 or y==0:
                K[x,y,n]=np.sum(L[x-2:x+1,y-2:y+1,n])/6.0
           #calculating averagefor all the other elements
           else:
                K[x,y,n]=np.sum(F[x-1:x+2,y-1:y+2,n])/9.0
    #normalising the matrix
    S=(K-K.min())/(K.max()-K.min())
    #plotting the nex figure
    plt.figure()
    plt.imshow(S)
    plt.title('Evgeniia Gleizer. Created by smooth.')
    plt.show() 
    return 
      
def compress(F2d,C):
  F = misc.imread("image.jpg")
  F2d=F[:,:,0]
  #computing svd algorithm
  M,S,N = np.linalg.svd(F2d, full_matrices=True)
  #diagonalising the middle matrix
  P=np.diag(S)

  #we want the sum of the sizes of the three matrixes to be less then C*m*n. 
  #So we nee to solve m*x+x*x+x*n=< C*m*n where m and are as follows
  m=F.shape[0]
  n=F.shape[1]
  #x is the solution of the quadratic above
  d=(m+n)**2-4*(1-C*m*n)
  x=int(-m-n+d**(0.5))/2 -1
  #calculation M,N,P of the new size
  P=P[0:x,0:x]
  M=M[:,0:x]
  N=N[0:x,:]
  return M,P,N


def reconstruct(a1,a2,a3,flag_display=False): 
    #multiplying the matrices back together to get F
    if flag_display==True:
        F=np.dot(a1,np.dot(a2,a3))
        return F


if __name__ == '__main__':
    #Read and display image
    F = misc.imread("image.jpg")
    plt.figure()
    plt.imshow(F)
    plt.show()


# The code below can be used to test if
# the code is functioning in a "reasonable" manner.

#    T,B,L,R = 0.1,0.2,0.3,0.4
#    G = crop(F,T,B,L,R)
#    plt.figure()
#    plt.imshow(G)
    
#    G = smooth(F)
#    plt.figure()
#    plt.imshow(G)
    
#    a1,a2,a3 = compress(F.mean(2),0.1)
#    G = reconstruct(a1,a2,a3)
#    plt.figure()
#    plt.imshow(G,'gray')
