
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

def least_fermat(N : int) -> int:
    '''
    Find the least Fermat power <= N 
    2 ** (2 ** n) <= N
    '''
    if N < 2:
        raise ValueError('The lowest Fermat power is 2')
    bit_log = N.bit_length() - 1
    bit_log_log = bit_log.bit_length()-1
    return 1 << (1 << bit_log_log)

def nim_product(a : int, b : int) -> int:
    # first handle trivial cases
    if a == 0 or b == 0:
        return 0
    elif a == 1:
        return b
    elif b == 1:
        return a
    elif a == 2 and b == 2:
        return 3
    else:
        # do euclidean division by greatest possible fermat power 
        # a = q_a * F_a + r_a and b = q_b * F_b + r_b
        F_a = least_fermat(a); F_b = least_fermat(b)
        q_a = a // F_a ; q_b = b // F_b
        r_a = a % F_a ; r_b = b % F_b

        # if one the Fermat powers is greater than the other, then
        # nim multiplication by it is the same as ordinary multiplication
        if F_a < F_b:
            return nim_product(a,q_b)*F_b ^ nim_product(a,r_b)
        elif F_a > F_b:
            return nim_product(q_a,b)*F_a ^ nim_product(r_a,b)
        else:
            # otherwise we have to distribute and use F_n ** 2 = 3 * F_n / 2
            p_1 = nim_product(q_a,q_b)
            p_2 = nim_product(r_a,r_b)
            p_3 = nim_product(q_a ^ r_a, q_b ^ r_b)
            p_4 = nim_product(p_1, F_a >> 1)
            p_5 = p_3 ^ p_2
            return p_5 * F_a ^ p_2 ^ p_4

def nim_power(x : int,n : int) -> int:
    if n == 0:
        return 1
    elif n == 1:
        return x
    elif x == 0:
        return 0
    elif x == 1:
        return 1
    else:
        i = 1
        prod = x
        while i < n:
            prod = nim_product(prod,x)
            i += 1
        return prod

def fermat_divisors(n : int, include_one : bool = False ) -> list:
    '''
    Find the divisors of a Mersenne number 2 ** (2 ** n) - 1
    By default does not include 1
    '''
    # for now, this only works for n < 6
    # could potentiall go up to n = 11 using known factors on wikipedia 
    # no one knows the factors of 2 ** (2 ** 11) + 1
    if n >= 6:
        raise ValueError('This function only works for n < 6')
    else:
        divisors = []
        for i in range(0 + int(not include_one),2 ** n):
            product = 1
            for j in range(n):
                if i >> j & 1:
                    product *= 2 ** (2 ** j) + 1
            divisors.append(product)
        return divisors
       

def nim_order(n : int) -> int:
    if n == 0:
        return 0
    elif n == 1:
        return 1
    elif n == 2:
        return 3
    elif n == 3:
        return 3
    elif n < 1 << (1 << 5):
        # make more efficient by only checking possible orders
        # use Lagrange's theorem
        # find the smallest field containing n i.e. smallest F_k > n
        exp = (n.bit_length() - 1).bit_length()
        # find the order of n must divide F_k - 1 which factors by difference of squares
        divisors = fermat_divisors(exp)
        for factor in divisors[:-1]:
            if nim_power(n,factor) == 1:
                return factor
        else:
            return divisors[-1] 
    else:
        for factor in fermat_divisors(5):
            if nim_power(n,factor) == 1:
                return factor
            # brute force: will probably loop forever
            i = 1 << (1 << 5)
            while True:
                if nim_power(n,i) == 1:
                    return i
                i += 1

nim_order(10)