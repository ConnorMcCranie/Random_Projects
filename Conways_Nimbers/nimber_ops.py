def exp(level : int) -> int:
    return 1 << level

def mask(level : int) -> int:
    return exp(exp(level)) - 1

def high_part(x : int, level : int) -> int:
    return x >> exp(level)

def low_part(x : int, level : int) -> int:
    return x & mask(x)

def combine(x : int, y : int, level : int) -> int:
    return x << exp(level) ^ y

def level(x : int) -> int:
    i = 0
    while high_part(x, i) != 0:
        i += 1
    return i

def half_mult(x : int, level : int) -> int:
    level -= 1
    if level < 0: return x
    high = high_part(x, level)
    low = low_part(x, level)
    return combine(half_mult(high ^ low, level),
                   half_mult(half_mult(high, level), level), level)
    
def nim_square(x : int, level : int) -> int:
    level -= 1
    if level < 0: return x
    high = high_part(x, level)
    low = low_part(x, level)
    high_square = nim_square(high, level)
    return combine(high_square, 
                   half_mult(high_square, level) ^ nim_square(low, level), 
                   level)
    
def nim_sqrt(x : int, level : int) -> int:
    level -= 1
    if level < 0: return x
    high = high_part(x, level)
    low = low_part(x, level)
    return combine(nim_sqrt(high, level),
            nim_sqrt(half_mult(high,level) ^ low, level),
            level)
    
def nim_times(x : int, y : int, level : int) -> int:
    level -= 1
    if level < 0: return x & y
    x_high = high_part(x, level)
    x_low = low_part(x, level)
    y_high = high_part(y, level)
    y_low = low_part(y, level)
    low_mult = nim_times(x_low, y_low, level)
    return combine(nim_times(x_high ^ x_low, y_high ^ y_low, level) ^ low_mult,
            half_mult(nim_times(x_high,y_high,level),level) ^ low_mult,
            level)
    
def nim_inverse(x : int, level : int) -> int:
    level -= 1
    if level < 0: return x
    high = high_part(x, level)
    low = low_part(x, level)
    ad = nim_times(high ^ low, low, level)
    bc = half_mult(nim_square(high, level), level)
    inv_det = nim_inverse(ad ^ bc, level)
    return combine(nim_times(high, inv_det, level),
            nim_times(high ^ low, inv_det, level),
            level)