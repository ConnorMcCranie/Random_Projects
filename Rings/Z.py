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
    
