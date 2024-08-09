
def is_power_of_two(n : int) -> bool:
    if n == 0:
        return False
    return (n & (n - 1)) == 0

def is_fermat_power(n : int) -> bool:
    if n < 4_000_000_000:
        return (n == 2) | (n == 4) | (n == 16) | (n == 256) | (n == 65536)
    elif not is_power_of_two(n):
        return False
    else:
        return is_power_of_two(n.bit_length() - 1)

def nim_sum(n : int,m : int) -> int:
    return n ^ m

def nim_prod(x : int,y : int) -> int:
    # first handle trivial cases
    if x == 0 | y == 0:
        return 0
    elif x == 1:
        return y
    elif y == 1:
        return x
    
    # next, use the rule for multiplying Fermat powers F=2 ** (2 ** n)
    # if x < F, then nim_prod(x,F)= x*F   (ordinary product)
    # and nim_prod(F,F) = 3*F/2
    m = min(x,y)
    M = max(x,y)
    
    if is_fermat_power(M):
        # nim product of Fermat power with smaller is ordinary product
        if m < M:
            return m*M
        else:
            # nim square of fermat power x is 3x/2
            return 3*M >> 1
    elif is_power_of_two(M):
        # if exponent is not power of 2, factor out 2's until it is
        # M = (factored) * (2 ** exponent)
        # we know at least one 2 needs to be pulled out        
        exponent = 1
        factored = M >> 1
        while not is_fermat_power(factored):
            factored >>= 1
            exponent += 1
        # now use formula for fermat power and associativity
        # m* M = (m * factored) * (2 ** exponent) = (m * (2 ** exponent)) * factored
        # we have to re-order the parentheses carefully to avoid infinite loop
        intermediate = nim_prod(m, factored)
        if intermediate == M:
            intermediate = nim_prod(m, 1 << exponent)
            return nim_prod(intermediate, factored)
        return nim_prod(intermediate, 1 << exponent)
    else:
        # otherwise, write it as the sum of powers of 2 and distribute
        sum = 0
        for index in range(M.bit_length()):
            if M >> index & 1 == 1:
                sum ^= nim_prod(m, 1 << index)
        return sum 
