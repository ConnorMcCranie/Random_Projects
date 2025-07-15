class Poly:
    def __init__(self, coeff, zero : int | object = 0, one : int | object =1):
        ''' Class for representing polynomials via lists. The coefficients
        can be any class which has __add__(), sub, __mul__() , _eq__( , 0) and 
        __pow__( , -1) dunder methods. If you want coefficients in regular
        integers or reals or complex, it's much more efficient to use 
        numpy, sympy, etc. instead
        
        self.zero is the zero element of the ring'''
        # trim off excess zeroes
        n = len(coeff)
        if n == 0: # treat empty list as zero polynomial
            coeff.append(zero)
            n = 1 
        is_zero = False
        for c in reversed(coeff):
            if c == zero:
                if n > 1:
                    n -= 1
                    continue
                is_zero = True
                break
            else:
                break          
        self.coeff = coeff[:n]
        self.deg = n-1
        self.is_zero = is_zero
        self.zero = zero
        self.one = one

    def __add__(self, other):
        assert self.zero == other.zero, 'need to have the same coefficient field'
        small, big = sorted([self, other], key=lambda x: x.deg)
        sum = [self.coeff[n] + other.coeff[n] for n in range(small.deg + 1)]
        sum += big.coeff[small.deg + 1:]
        return Poly(sum, zero=self.zero, one=self.one)

    def __call__(self, num):
        result = self.zero
        for a in reversed(self.coeff):
            result = result * num
            result = result + a
        return result

    def __eq__(self, other):
        try:
            return self.coeff == other.coeff and self.zero == other.zero
        except:
            return False

    def __floordiv__(self, other):
        assert self.zero == other.zero, 'need to have the same coefficient field'
        numerator = Poly(self.coeff, self.zero, self.one) # copy
        quotient = Poly([self.zero], self.zero, self.one)
        def monomial(n, const):
            coeff = [self.zero for _ in range(n+1)]
            coeff[n] = const
            return Poly(coeff, self.zero, self.one)
        while numerator.deg >= other.deg and not numerator.is_zero:
            coeff = numerator.coeff[-1] * other.coeff[-1] ** -1
            term = monomial(numerator.deg - other.deg, coeff)

            quotient = quotient + term
            numerator = numerator - term * other
        return quotient

    def __getitem__(self, item):
        return self.coeff[item]

    def __iter__(self):
        return self.coeff.__iter__()

    def __mod__(self, other):
        return self - other * (self // other)

    def __mul__(self, other):
        try: # scalar multiplication
            result = [self.coeff[n] * other for n in range(self.deg + 1)]
            return Poly(result, self.zero, self.one)
        except:
            deg = self.deg + other.deg
            p, q = self.coeff, other.coeff
            prod = [self.zero for _ in range(deg + 1)]
            for n in range(deg + 1):
                for k in range(n + 1):
                    if k < len(p) and n - k < len(q):
                        prod[n] = prod[n] + p[k] * q[n - k]
            return Poly(prod, self.zero, self.one)     

    def __pow__(self, p : int):
        copy = self
        if p >= 0: # binary exponentiation by squaring
            result = Poly([self.one], self.zero, self.one) # identity
            while p > 0:
                if p & 1:
                    result = result * copy
                copy = copy * copy
                p >>= 1
            return result
        else:
            raise ValueError('power needs to be a non-negative integer')      

    def __repr__(self):
        return str(self.coeff)
        
    def __sub__(self, other):
        assert self.zero == other.zero, 'need to have the same coefficient field'
        small, big = sorted([self, other], key=lambda x: x.deg)
        pad = [small.coeff[n] for n in range(small.deg + 1)]
        pad += [self.zero for _ in range(big.deg - small.deg)]
        if self.deg <= other.deg: 
            diff = [pad[n] - other.coeff[n] for n in range(len(pad))]
        else:
            diff = [self.coeff[n] - pad[n] for n in range(len(pad))]
        return Poly(diff, self.zero, self.one)
    
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
    
class PolyMod(Poly):
    ''' Polynomials with coefficients in finite cyclic group'''
    def __init__(self, data : list[int] | int, mod : int):
        if type(data) == list:
            mod_coeff = [Z(mod, n) for n in data]
        elif type(data) == int:
            assert data >= 0, 'integer input should be non-negative'
            def base_decomp(base : int, num : int) -> list[int]:
                ''' write an integer in a given base'''
                decomp = []
                while num > 0:
                    decomp.append(num % base)
                    num //= base
                return decomp
            mod_coeff = [Z(mod, n) for n in base_decomp(mod, data)]
        else:
            raise ValueError('Input is either list or integer')
        zero, one = Z(mod, 0), Z(mod, 1)
        super().__init__(mod_coeff, zero, one)
        reduced_coeff = [c.repr for c in mod_coeff]

        self.id = Poly(reduced_coeff)(mod)
        self.mod = mod
    
    def __add__(self, other):
        assert self.mod == other.mod, 'need same base field'
        sum = super().__add__(other)
        coeff = [c.repr for c in sum.coeff]
        return PolyMod(coeff, self.mod)
    
    def __floordiv__(self, other):
        assert self.mod == other.mod, 'need same base field'
        div = super().__floordiv__(other)
        coeff = [c.repr for c in div.coeff]
        return PolyMod(coeff, self.mod)
    
    def __mod__(self, other):
        assert self.mod == other.mod, 'need same base field'
        remainder = super().__mod__(other)
        coeff = [c.repr for c in remainder.coeff]
        return PolyMod(coeff, self.mod)
    
    def __mul__(self, other):
        assert self.mod == other.mod, 'need same base field'
        prod = super().__mul__(other)
        coeff = [c.repr for c in prod.coeff]
        return PolyMod(coeff, self.mod)
    
    def __pow__(self, p : int):
        power = super().__pow__(p)
        coeff = [c.repr for c in power.coeff]
        return PolyMod(coeff, self.mod)
    
    def __sub__(self, other):
        assert self.mod == other.mod, 'need same base field'
        sub = super().__sub__(other)
        coeff = [c.repr for c in sub.coeff]
        return PolyMod(coeff, self.mod)       

    def factor(self, degree : int):
        '''uses brute force approach to find the least monic polynomial 
        of a given degree dividing `self`, or `None` if no such factor
        exists '''
        assert type(degree)==int and degree >=0, (
            'degree of factor should be non-negative integer')
        if degree == 0: return PolyMod([1], self.mod)
        elif degree > self.deg: return None
        else:
            p, d = self.mod, degree
            for n in range(p**d + 1, 2 * p**d):
                candidate = PolyMod(n, p)
                remainder = self % candidate
                if remainder.is_zero:
                    return candidate
            return None

    def factors(self, degree : int) -> list:
        result = []
        term = self
        if degree == 0: return [PolyMod([1], self.mod)]
        for iter in range(self.deg // degree):
            factor = term.factor(degree)
            if factor != None:
                result.append(factor)
                term = term // factor
            else:
                break
        return result
    
    def gcd(self, other) -> list:
        ''' Find the bezout coefficents x * self + y * other = gcd
        returns [x, y, gcd]'''
        # Initialize: [x, y, value]
        zero, one = [PolyMod([i], self.mod) for i in range(2)]
        x0, y0, r0 = one, zero, self
        x1, y1, r1 = zero, one, other
        
        while not r1.is_zero: # do (extended) euclidean algorithm
            q = r0 // r1
            x0, x1 = x1, x0 - q * x1
            y0, y1 = y1, y0 - q * y1
            r0, r1 = r1, r0 - q * r1        
        return [x0, y0, r0]
    
    def is_irred(self) -> bool:
        # requires prime modulus
        n, p = self.deg, self.mod
        if n == 0: return False
        elif n == 1: return True
        else:
            divs = [k for k in range(1, int(n**.5 + 1)) if n % k == 0]
            def frob(k): # x^(p^k) - x
                poly = [0 for _ in range(p**k + 1)]
                poly[p**k], poly[1] = 1, -1
                return PolyMod(poly, p)
            for div in divs:
                poly = frob(div)
                d = self.gcd(poly)[2]
                if d.deg > 0: return False
            return (frob(n) % self).is_zero
