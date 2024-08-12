'''
use vandermond matrix to find a polynomal that interpolates a set of points
'''
import numpy as np
from math import floor
from math import ceil
from Polynomial import Polynomial as Poly

def vandermonde(x : np.ndarray) -> np.ndarray:
    '''
    create the vandermonde matrix for a set of points
    '''
    n = len(x)
    A = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            A[i,j] = x[i]**j
    return A

def interpolate(x: np.ndarray, y:np.ndarray , name = 'p', ) -> Poly:
    '''
    interpolate a set of points
    '''
    if len(x) != len(y):
        raise ValueError('x and y must have the same length')
    for i in range(len(x)):
        for j in range(len(x)):
            if i != j:
                if x[i] == x[j] and y[i] != y[j]:
                    raise ValueError('That\'s not a function!')
    
    A = vandermonde(x)
    coeff =  np.linalg.inv(A) @ y
    return Poly(coeff, name)

def factorial(n : int) -> int:
    '''
    compute the factorial of a number
    '''
    if n == 0:
        return 1
    else:
        return n*factorial(n-1)

def  vand_exp(x : np.ndarray) -> np.ndarray:
    '''
    create the vandermonde matrix weighted by factorials for a set of points
    '''
    n = len(x)
    A = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            A[i,j] = x[i]**(j)/factorial(j)
    return A

x = np.array([n for n in range(1,21)])
y = np.array([1,1,1,2,2,2,3,3,4,4,5,5,6,7,8,9,10,12,15,20])
t = interpolate(x,y,'t' )
int_coeffs = [round(c,8) for c in t.coeffs]
T = Poly(int_coeffs, 'T')
#print(t.poly_str)

# prime counting function
def prime_count(n : int) -> int:
    '''
    count the number of primes less than or equal to n
    '''
    if n < 2:
        return 0
    primes = [2]
    for i in range(3,n+1):
        is_prime = True
        for p in primes:
            if i % p == 0:
                is_prime = False
                break
        if is_prime:
            primes.append(i)
    return len(primes)

print([n - prime_count(n) + floor(n/9)*(ceil(np.cos(n/np.e)+np.sin(n/420))-1)+floor(n/17)*((n-16)**2+2*0 **((n)%4)) for n in range(1,21)])