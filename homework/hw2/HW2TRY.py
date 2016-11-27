import numpy as np
import matplotlib.pyplot as plt
from scipy import misc
from scipy import linalg
if __name__ == '__main__':
    print""
def compress(F2d,C):
  """ Construct a compression of 2-d image matrix F2d which
  requires C times the memory needed for F2d"""
  F = misc.imread("image.jpg")
  F2d=F[:,:,0]
  M,S,N = np.linalg.svd(F2d, full_matrices=True)
  P=np.diag(S)

  #we want the sum of the sizes of the three matrixes to be less then C*m*n. 
  #So m*x+x*x+x*n=< C*m*n
  m=F.shape[0]
  n=F.shape[1]
  d=(m+n)**2-4*(1-C*m*n)
  x=int(-m-n+d**(0.5))/2 -1
  P=P[0:x,0:x]
  M=M[:,0:x]
  N=N[0:x,:]
  return M,P,N
