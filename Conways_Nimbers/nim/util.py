from transfinite import Ordinal


def base(n: int, b: int, length: int = 0) -> list[int]:
    """returns the base b expansion of n as a  list
    n = sum_j^N a_j*b^j returns [a_0, a_1, ..., a_N]
    optional length parameter to include leading zeros if desired,
    otherwise minimum possible length where highest digit is non-zero"""
    assert n >= 0 and b > 1, (
        "positional base expansion defined for"
        "non-negative integers and positive bases"
    )
    if length == 0:
        if n == 0:
            return [0]
        coeffs = []
        while n > 0:
            coeffs.append(n % b)
            n //= b
        return coeffs
    else:
        coeffs = base(n, b)
        while len(coeffs) < length:
            coeffs.append(0)
        return coeffs


def base_eval(coeffs: list[int], b: int) -> int:
    total = 0
    for coeff in coeffs[::-1]:
        total = coeff + b * total
    return total


def ord_decomp(ordinal: Ordinal | int) -> list:
    """returns [inf_1, inf_2, ..., inf_n, finite term]"""
    if isinstance(ordinal, int):
        assert ordinal >= 0
        return [ordinal]
    high = Ordinal(ordinal.exponent, ordinal.coefficient)
    terms = [high]
    remainder = ordinal.addend
    while isinstance(remainder, Ordinal):
        terms.append(Ordinal(remainder.exponent, remainder.coefficient, 0))
        remainder = remainder.addend
    assert isinstance(remainder, int)
    terms.append(remainder)
    return terms


def ord_recomp(list: list) -> Ordinal | int:
    result = 0
    terms = sorted(list, reverse=True)
    for term in list:
        result = result + term
    return result


w = Ordinal()

small_primes = [
    2,
    3,
    5,
    7,
    11,
    13,
    17,
    19,
    23,
    29,
    31,
    37,
    41,
    43,
    47,
    53,
    59,
    61,
    67,
    71,
    73,
    79,
    83,
    89,
    97,
    101,
    103,
    107,
    109,
    113,
    127,
    131,
    137,
    139,
    149,
    151,
    157,
    163,
    167,
    173,
    179,
    181,
    191,
    193,
    197,
    199,
    211,
    223,
    227,
    229,
    233,
    239,
    241,
    251,
    257,
    263,
    269,
    271,
    277,
    281,
    283,
    293,
    307,
    311,
    313,
    317,
    331,
    337,
    347,
    349,
    353,
    359,
    367,
    373,
    379,
    383,
    389,
    397,
    401,
    409,
    419,
    421,
    431,
    433,
    439,
    443,
    449,
    457,
    461,
    463,
    467,
    479,
    487,
    491,
    499,
    503,
    509,
    521,
    523,
    541,
    547,
    557,
    563,
    569,
    571,
    577,
]

alpha_p = {
    2: 3,
    3: 2,
    5: 4,
    7: w + 1,
    11: w**w + 1,
    13: w + 4,
    17: 16,
    19: w**3 + 4,
    23: w**w**3 + 1,
    29: w**w**2 + 4,
    31: w**w + 1,
    37: w**3 + 4,
    41: w**w + 1,
    43: w**w**2 + 1,
    47: w**w**7 + 1,
    53: w**w**4 + 1,
    59: w**w**8 + 1,
    61: w**w + w,
    67: w**w**3 + w,
    71: w**w**2 + w**w,
    73: w**3 + 1,
    79: w**w**4 + 1,
    83: w**w**11 + 1,
    89: w**w**3 + 1,
    97: w + 256,
    101: w ** (w * 5) + 1,
    103: w**w**5 + w,
    107: w**w**14 + 1,
    109: w**3 + 4,
    113: w**w**2 + 4,
    127: w**w**2 + 1,
    131: w**w**4 + w**w,
    137: w**w**5 + 1,
    139: w**w**7 + w,
    149: w**w**10 + 1,
    151: w**w + w,
    157: w**w**4 + 1,
    163: w**27 + 4,
    167: w**w**21 + 1,
    173: w**w**12 + 4,
    179: w**w**22 + 1,
    181: w**w + w**3,
    191: w**w**6 + w**w,
    193: w + 65536,
    197: w ** (w**2 * 7) + 4,
    199: w**w**3 + w**3,
    211: w**w**2 + w**w,
    223: w**w**10 + 1,
    227: w**w**28 + 1,
    229: w**w**6 + 1,
    233: w**w**8 + 1,
    239: w**w**5 + w**w**2,
    241: w + 16,
    251: w ** (w * 5) + 1,
    257: 256,
    263: w**w**30 + 1,
    269: w**w**17 + 1,
    271: w**w + w**9,
    277: w**w**7 + 1,
    281: w**w**2 + w**w,
    283: w**w**13 + 1,
    293: w**w**19 + 4,
    307: w**w**5 + w,
    311: w**w**9 + 1,
    313: w**w**4 + 1,
    317: w**w**20 + 1,
    331: w**w + w,
    337: w**w**2 + 1,
    347: w**w**38 + 1,
    349: w**w**8 + 1,
    353: w**w**3 + 16,
    359: w**w**39 + 1,
    367: w**w**16 + 1,
    373: w**w**9 + w,
    379: w**w**2 + w**9,
    383: w**w**41 + 1,
    389: w**w**23 + 1,
    397: w**w**3 + 1,
    401: w ** (w * 5) + 16,
    409: w**w**5 + w,
    419: w**w**6 + w**w**3,
    421: w**w**2 + w**w,
    431: w**w**12 + 1,
    433: w**3 + 16,
    439: w**w**19 + 1,
    443: w**w**5 + w**w**4,
    449: w**w**2 + 65536,
    457: w**w**6 + 1,
    461: w**w**7 + 1,
    463: w**w**3 + w**w**2,
    467: w**w**49 + 1,
    479: w**w**50 + 1,
    487: w**81 + 1,
    491: w ** (w**2 * 7) + w**w,
    499: w**w**21 + 1,
    503: w**w**52 + 1,
    509: w**w**29 + 4,
    521: w**w**4 + w**w,
    523: w**w**8 + w**3,
    541: w**w + w**9,
    547: w**w**4 + w**w**2,
    557: w**w**32 + 1,
    563: w**w**58 + 1,
    569: w**w**18 + 1,
    571: w**w**6 + 1,
    577: w**3 + 256,
}
