class Z:
    def __init__(self, mod : int, num : int = 0, ):
        ''' A class for representing finite cyclic groups.'''
        assert type(mod) == int and mod > 0, 'modulus is positive integer'

        self.repr = num % mod
        self.mod = mod
        self.inv = ... # sentinel since `None` will mean not inverible
    
    def __add__(self, other) -> object:
        assert self.mod == other.mod, 'adding in cyclic group '\
                                      'needs to have same modulus'
        result = self.repr + other.repr % self.mod
        return Z(self.mod, result)
    
    def __eq__(self, other) -> bool:
        try: # if other is also a Z(n) object
            if self.repr == other.repr and self.mod == other.mod:
                return True
            else:
                return False
        except: # allow comparison to integer
            if type(other) == int: 
                return (self.repr - other % self.mod) == 0
            else: # all others return False
                return False
    
    def __mul__(self, other) -> object:
        try:
            if self.mod == other.mod:
                return Z(self.mod, self.repr * other.repr % self.mod)
            else:
                print('need to have same modulus to multiply')
                return None
        except:
            assert type(other) == int, 'multiplication is supported with another'\
            'element of cyclic group or with an integer'
            return Z(self.mod, self.repr * other % self.mod)

    def __pow__(self, power : int) -> object:
        return Z(self.mod, pow(self.repr, power, self.mod))

    def __repr__(self):
        return str(self.repr)
    
    def __sub__(self, other):
        return Z(self.mod, self.repr - other.repr % self.mod)

    def inverse(self) -> int | None:
        ''' Returns the inverse if it exists or None if it doesn't'''
        if self.inv == ... : # not already computed
            try:
                self.inv = pow(self.repr, -1, self.mod)
            except:
                self.inv = None
        return self.inv
    
def miller_rabin(n, k=20):
    """
    Miller-Rabin primality test
    Returns True if n is probably prime, False if n is composite
    k is the number of rounds of testing
    """
    import random

    if n < 2:
        return False
    if n == 2 or n == 3:
        return True
    if n % 2 == 0:
        return False
    
    # Write n-1 as d * 2^r
    d = n - 1
    r = 0
    while d % 2 == 0:
        d //= 2
        r += 1
    
    # Perform k rounds of testing
    for _ in range(k):
        a = random.randint(2, n - 2)
        x = pow(a, d, n)
        
        if x == 1 or x == n - 1:
            continue
        
        for _ in range(r - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                break
        else:
            return False
    
    return True



def is_prime(n):
    """
    Deterministic primality test for small numbers, probabilistic for large numbers
    """
    if n < 2:
        return False
    small_primes = {
                2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 
                61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 
                131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193,
                197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 
                263, 269, 271, 277, 281, 283, 293
                }
    if n in small_primes:
        return True
    for p in small_primes:
        if n % p == 0: 
            return False
    
    # Use deterministic test for small numbers
    if 300 < n < 5000:
        for i in range(3, int(n**0.5) + 1, 2):
            if n % i == 0:
                return False
        return True
    
    # Use Miller-Rabin for larger numbers
    return miller_rabin(n)

def integer_nth_root(n, k):
    """
    Compute the integer k-th root of n using binary search
    Returns the largest integer r such that r^k <= n
    """
    if n == 0:
        return 0
    if n == 1:
        return 1
    if k == 1:
        return n
    
    # Binary search for the k-th root
    low = 1
    high = int(n**(1.0/k)) + 2  # Add buffer for floating point errors
    
    while low <= high:
        mid = (low + high) // 2
        mid_k = mid ** k
        
        if mid_k == n:
            return mid
        elif mid_k < n:
            low = mid + 1
        else:
            high = mid - 1
    
    return high

def is_power(n, k):
    """
    Check if n is a perfect k-th power
    Returns the k-th root if it exists, False otherwise
    """
    if n <= 0:
        return None
    if n == 1:
        return 1
    
    root = integer_nth_root(n, k)
    if root ** k == n:
        return root
    return False

def int_log(N : int, base : int) -> int:
    ''' Returns p such that base ^ p <= N < base ^ (p+1)'''
    def level(N : int, p : int = 2) -> int:
        '''
        Returns the largest 'level' L (w.r.t. p) such that 
        p ** (2 ** L) <= N
        '''
        L = 0
        if N < p:
            return 0
        while N // (p ** (1 << L)) != 0 :
            L += 1
        return L-1
    
    if N < base:
        return 0
    # for the highest power 'exp' s.t. base^exp <= N, find
    # the largest power of 2 that is less than or equal to exp
    total = 1 << level(N, base)
    # now divide N to recursively find the binary expansion of exp
    N = N // (base ** total)
    while N != 0:
        total += 1 << level(N, base)
        N = N // (base ** (1 << level(N, base)))
    return total - 1

def prime_power(n : int) -> list | bool:
    """
    Determine if n = p^k for some prime p and integer k > 0
    Returns (p, k) if n is a prime power, None otherwise
    
    Algorithm:
    1. Handle special cases (n <= 1)
    2. For each possible exponent k from 2 to log₂(n):
       - Check if n is a perfect k-th power
       - If so, test if the k-th root is prime
    3. If no prime power found, check if n itself is prime (k=1 case)
    """
    if n <= 1:
        return False
    
    if n == 2:
        return [2, 1]
    
    # Check for prime powers p^k where k >= 2
    max_exponent = int_log(n, 2)
    
    for k in range(1, max_exponent + 1):
        root = is_power(n, k)
        if root:
            if is_prime(root):
                return [root, k]
    return False