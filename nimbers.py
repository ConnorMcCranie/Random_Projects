
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
            p_4 = nim_product(q_a, F_a >> 1)
            p_5 = p_3 ^ p_2
            return p_5 * F_a ^ p_2 ^ p_4

print(nim_product(2,4)) #