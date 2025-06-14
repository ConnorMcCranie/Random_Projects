class Perm:
    ''' permutations represented by lists'''
    def __init__(self, func : list):
        assert set(func) == {n for n in range(len(func))}, (
            'That list is not a permutation')
        self.func = func

    def __getitem__(self, item : int) -> int:
        return self.func[item]

    def __len__(self) -> int:
        return len(self.func)

    def __mul__(self, other):
        '''multiplication of permutations (right first, then left)'''
        assert len(self) == len(other), (
        'permutations need to be the same length to multiply')
        prod = [self[other[n]] for n in range(len(self)) ]
        return Perm(prod)

    def __pow__(self, pow : int):
        copy = Perm(self.func)
        if pow >= 0: # binary exponentiation by squaring
            result = Perm([n for n in range(len(self))]) # identity
            while pow > 0:
                if pow & 1:
                    result = result * copy
                copy = copy * copy
                pow >>= 1
            return result
        elif pow == -1: # inverse
            inv = [0 for _ in range(len(self))]
            for i, n in enumerate(self.func):
                inv[n] = i
            return Perm(inv)
        elif pow < -1:
            return (self ** -1) ** -pow
        else:
            raise ValueError('power needs to be an integer')

    def __repr__(self) -> str:
        return str(self.func)