from transfinite import Ordinal
import pandas as pd
from util import base, base_eval, ord_decomp, ord_recomp, small_primes, alpha_p

W = Ordinal

# utility functions: collect into util.py later


class Nim:
    """nimbers"""

    def __init__(self, n: int | Ordinal) -> None:
        """ordinal considered an a field element in On_2
        val = ordinal
        field = smallest x > n such that x is a field
        base = largest y < n such that y is a field (or base = 0 for n < 2)
        write n = high * base + low, where low,high < base
        """

        def exp2(level: int) -> int:
            return 1 << level

        if isinstance(n, int):
            self.val = abs(n)
            self.isfinite = True
            if self.val < 2:
                self.field = 2  # smallest
                self.base = 0
                self.high = 0
                self.low = self.val
                self.level = 0
            else:
                level = 0
                while n >> (1 << level) > 0:
                    level += 1
                self.level = level
                self.field = exp2(exp2(level))
                self.base = exp2(exp2(level - 1))
                self.high = self.val // self.base
                self.low = self.val - self.high * self.base
        else:
            self.val = n
            self.isfinite = False
            # these attributes can be found, but not as useful for infinite nums
            self.field = None
            self.base = None
            self.high = None
            self.low = None
            # assert isinstance(
            #     n, Ordinal), 'An infinite nimber must be an ordinal'
            # self.val = n
            # self.isfinite = False
            # # find the smallest field containing n
            # if (k := n.exponent) < Ordinal():  # n = omega^k, k fintie => cubic extension
            #     power = 0
            #     while 3 ** power <= k:
            #         power += 1
            #     self.field = Ordinal(3**power)
            #     exp = 3**(power-1)
            #     self.base = Ordinal(exp)
            #     high = Ordinal(n.exponent - exp, n.coefficient,
            #                 0) if exp < n.exponent else n.coefficient
            #     remainder = n.addend
            #     while isinstance(remainder, Ordinal) and remainder.exponent > exp:
            #         high += Ordinal(remainder.exponent - exp,
            #                         remainder.coefficient, 0)
            #         remainder = remainder.addend
            #         if isinstance(remainder, Ordinal) and remainder.exponent == exp:
            #             high += remainder.coefficient
            #             remainder = remainder.addend
            #     self.high = high
            #     self.low = remainder

    def __add__(self, other):
        if self.isfinite and other.isfinite:
            return Nim(self.val ^ other.val)  # finite nim sum is bitwise XOR
        if self.val == other.val:
            return Nim(0)
        ord1, ord2 = self.val, other.val
        terms1, terms2 = ord_decomp(ord1), ord_decomp(ord2)
        sum = {}
        for term in terms1[:-1]:
            sum[term.exponent] = term.coefficient
        for term in terms2[:-1]:  # nim sum the coeffs if any terms with same exp
            try:
                sum[term.exponent] = sum[term.exponent] ^ term.coefficient
            except:
                sum[term.exponent] = term.coefficient

        keys = sorted(sum.keys())  # ordinal addition not commutative
        result = terms1[-1] ^ terms2[-1]  # nim sum of finite part
        for key in keys:
            if sum[key] > 0:  # add bigger on the left
                result = Ordinal(key, sum[key]) + result
        return Nim(result)

    def __eq__(self, other):
        return self.val == other.val

    def __hash__(self):
        return self.val.__hash__()

    def __mul__(self, other):
        assert isinstance(other, Nim)
        x, y = self.val, other.val
        if x == 0 or y == 0:
            return Nim(0)
        if x == 1:
            return other
        if y == 1:
            return self
        if self.isfinite and other.isfinite:
            if self.val == other.val:
                return self.sq()

            def nim_product(a: int, b: int) -> int:
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
                    F_a, q_a, r_a = Nim(a).base, Nim(a).high, Nim(a).low
                    F_b, q_b, r_b = Nim(b).base, Nim(b).high, Nim(b).low

                    # if one the Fermat powers is greater than the other, then
                    # nim multiplication by it is the same as ordinary multiplication
                    if F_a < F_b:
                        return nim_product(a, q_b) * F_b ^ nim_product(a, r_b)
                    elif F_a > F_b:
                        return nim_product(q_a, b) * F_a ^ nim_product(r_a, b)
                    else:
                        # otherwise we have to distribute and use F_n ** 2 = 3 * F_n / 2
                        p_1 = nim_product(q_a, q_b)
                        p_2 = nim_product(r_a, r_b)
                        p_3 = nim_product(q_a ^ r_a, q_b ^ r_b)
                        p_4 = nim_product(p_1, F_a >> 1)
                        p_5 = p_3 ^ p_2
                        return p_5 * F_a ^ p_2 ^ p_4

            return Nim(nim_product(x, y))
        elif self.isfinite and not other.isfinite:
            terms = ord_decomp(other.val)  # distribute to each coefficient
            terms[-1] = (self * Nim(terms[-1])).val
            for term in terms[:-1]:
                term.coefficient = (self * Nim(term.coefficient)).val
            return Nim(ord_recomp(terms))
        elif not self.isfinite and other.isfinite:
            return other * self
        else:  # both infinite
            terms1 = ord_decomp(self.val)
            if self.val == other.val and (len(terms1) > 2 or terms1[-1] > 0):
                return self.sq()
            terms2 = ord_decomp(other.val)
            inf1, fin1 = terms1[:-1], terms1[-1]
            inf2, fin2 = terms2[:-1], terms2[-1]
            if fin1 == 0 and fin2 == 0:
                to_sum = {Nim(0)}
                # result = Nim(0)

            else:
                # start by "FOIL-ing" to handle the terms where one is finite
                to_sum = (
                    {Nim(fin1) * other} ^ {self * Nim(fin2)} ^ {Nim(fin1) * Nim(fin2)}
                )
                # result = Nim(fin1) * other + self * Nim(fin2) \
                # + Nim(fin1) * Nim(fin2)
            for x in inf1:  # expand and distribute the purely infinite terms
                for y in inf2:
                    # calculate (w^N * a) x (w^M * b)
                    X, Y = sorted([x, y], reverse=True)  # X >= Y
                    N, M = X.exponent, Y.exponent  # N >= M
                    a, b = X.coefficient, Y.coefficient
                    coeff = Nim(a) * Nim(b)

                    if N < Ordinal() and M < Ordinal():  # handle finite case
                        N_tern = base(N, 3)  # write exponents in ternary
                        M_tern = base(M, 3)
                        K = max([len(N_tern), len(M_tern)])

                        # invert the powers of 3 since w^(3^k) = 2^(3^{-k-1})
                        def phi(n):
                            return base_eval(base(n, 3, K)[::-1], 3)  # reversed 3s

                        exp_p = phi(N) + phi(M)  # the exponent 2^(3^{-K} * exp_2)

                        next_power = 3**K
                        q, r = exp_p // next_power, exp_p % next_power
                        coeff = coeff * Nim(2) ** q
                        W = Nim(Ordinal(phi(r))) if r > 0 else Nim(1)
                        to_sum ^= {coeff * W}
                        # result = result + coeff * W
                    elif N < Ordinal() and not M < Ordinal():
                        W = Nim(Ordinal(N + M))
                        to_sum ^= {coeff * W}
                        # result = result + coeff * W
                    else:
                        assert max(N, M) < Ordinal(
                            105
                        ), "Nimber multiplication is only implemented for ordinals < w**w**105"
                        N_decomp, M_decomp = ord_decomp(N), ord_decomp(M)
                        # multiply all omega terms using exponent rules
                        N_fin_exp, M_fin_exp = N_decomp[-1], M_decomp[-1]
                        N_inf_exp, M_inf_exp = N_decomp[:-1], M_decomp[:-1]
                        # handle finite exponent first
                        omega_N = Ordinal(N_fin_exp) if N_fin_exp != 0 else 1
                        omega_M = Ordinal(M_fin_exp) if M_fin_exp != 0 else 1
                        # start multiplying omegas
                        prod: Ordinal = (Nim(omega_N) * Nim(omega_M)).val
                        # group the exponents by like terms
                        N_dict = {exp.exponent: exp for exp in N_inf_exp}
                        M_dict = {exp.exponent: exp for exp in M_inf_exp}
                        N_exp_exp = set(N_dict)
                        M_exp_exp = set(M_dict)
                        like_terms = N_exp_exp & M_exp_exp
                        unlike_terms = (N_exp_exp | M_exp_exp) - like_terms
                        for t in sorted(unlike_terms):  # normal ordinal product
                            mult = N_dict[t] if t in N_dict else M_dict[t]
                            prod = Ordinal(mult) * prod  # pull through
                        for n in like_terms:
                            i = N_dict[n].coefficient
                            j = M_dict[n].coefficient
                            p = small_primes[n + 1]
                            alpha = Nim(alpha_p[p])
                            i_p = base(i, p)  # write coeff in base p
                            j_p = base(j, p)
                            K = max([len(i_p), len(j_p)])

                            # invert the powers of p
                            def phi(n):
                                return base_eval(base(n, p, K)[::-1], p)

                            exp_p = phi(i) + phi(j)

                            next_power = p**K
                            q, r = exp_p // next_power, exp_p % next_power
                            coeff = coeff * alpha**q
                            term = Ordinal(Ordinal(n) * phi(r)) if r > 0 else 1
                            # need recursive because might be new like terms
                            prod = (Nim(prod) * Nim(term)).val
                        to_sum ^= {coeff * Nim(prod)}
                        # result = result + coeff * Nim(prod)
            return Nim(0).sum(*to_sum)

    def __pow__(self, p):
        """
        Compute x**n using exponentiation by squaring.

        """
        if p >= 0:  # binary exponentiation by squaring
            result = Nim(1)
            nimber = self
            while p > 0:
                if p & 1:
                    result = result * nimber
                nimber = nimber.sq()
                p >>= 1
            return result
        elif p == -1:
            if self.field == 2:
                return self
            a, b, F, f = (
                Nim(self.high),
                Nim(self.low),
                Nim(self.base),
                Nim(self.base >> 1),
            )
            det = (a + b) * b + a * a * f
            return det ** (-1) * (a * F + (a + b))
        else:
            inv = self ** (-1)
            return inv ** (-p)

    def __repr__(self) -> str:
        if self.isfinite:
            return str(self.val)
        else:
            return self.val.__repr__()

    def _repr_latex_(self):
        """
        Special method for Jupyter to render LaTeX.
        """
        if self.isfinite:
            # No special LaTeX for integers, just return the string
            return f"${self.val}$"
        else:
            # Delegate to the Ordinal's LaTeX representation
            return self.val._repr_latex_()

    def deg(self) -> int:
        """returns the degree of the minimal polynomial"""
        if self.isfinite:
            return self.field.bit_length() - 1
        else:
            d = 1
            x = self * self
            while x != self:
                x = x * x
                d += 1
            return d

    def det(self):
        """det(N) = determinant of the multiplication by N matrix
        This is the same as the field norm over the next smallest field.
        det(x) = x^(x.base + 1) = x^(2^2^n + 1)  if 2^2^n <= x < 2^2^(n+1)
        """
        if self.isfinite:
            if self.field == 2:
                return self
            a, b, F = Nim(self.high), Nim(self.low), Nim(self.base >> 1)
            return (a + b) * b + a.sq() * F  # much faster than self**(self.base+1)

    def det_star(self):
        """number of times x -> det(x) is repeated until landing in F_2"""
        iter = 0
        if self.level == 0:
            return iter
        d = self.det()
        while d != Nim(1):
            iter += 1
            d = d.det()
        return iter + 1

    def inv(self):
        return self ** (-1)

    def is_gen(self) -> bool:
        """True iff order(self) == self.field - 1"""
        if self.isfinite:
            """is_gen ==> det* = self.level
            det* = self.level ==> is_gen if level < 6"""
            if self.val == 0:
                return False
            if self.level != self.det_star():
                # is_gen ==> det* = self.level
                return False
            elif self.level < 6:
                # self.base = 2^2^(n-1) + 1 is prime for n < 6
                return True
            elif self.level == 6:  # F_5 = 641 × 6,700,417
                level_small = (self**641).level
                level_big = (self**6700417).level
                return min(level_small, level_big) == 6  # both are divisors of order
            elif self.level == 7:  # F_7 = 274,177 × 67,280,421,310,721
                level_small = (self**274177).level
                level_big = (self**67280421310721).level
                return min(level_small, level_big) == 7 and self.det().is_gen()
            # could keep going, but the highest factored is 'only' F_11 anyways
            else:
                return self.order() == (1 << (1 << self.level)) - 1

    def order(self):
        if self.isfinite:
            if self.level <= 7:
                if self.is_gen():  # avoid infinite loop
                    return (1 << (1 << self.level)) - 1
            # def fermat_divisors(n: int, include_one: bool = False) -> list:
            #     '''
            #     Find the divisors of a Mersenne number 2 ** (2 ** n) - 1
            #     By default does not include 1
            #     '''
            #     # for now, this only works for n < 6
            #     # could potentiall go up to n = 11 using known factors on wikipedia
            #     # no one knows the factors of 2 ** (2 ** 11) + 1
            #     if n >= 6:
            #         raise ValueError('This function only works for n < 6')
            #     else:
            #         divisors = []
            #         for i in range(0 + int(not include_one), 2 ** n):
            #             product = 1
            #             for j in range(n):
            #                 if i >> j & 1:
            #                     product *= 2 ** (2 ** j) + 1
            #             divisors.append(product)
            #         return divisors

            n = self.val
            L = self.level
            if L == 0:
                return n
            elif L == 1:
                return 3  # ord(2) = ord(3) = 3
            elif L <= 5:  # order divides 3 * 5 * 17 * 257 * 65_537
                return (self.base + 1) * self.det().order()
            else:
                factor = 1
                nimber = self
                while nimber.level > 5:
                    nimber = nimber.det()
                    factor *= nimber.base + 1

                prime_divs = nimber.order()
                nimber_to_primes = self**prime_divs  # order divides F_5*...*F_{L-1}
                factor *= prime_divs
                max_non_gen = (1 << (1 << L) - 1) // 3
                print(
                    f"Warning: Order of {self} may take forever to calculate:\
                        too difficult to factor {factor}"
                )
                for i in range(641, max_non_gen // prime_divs + 1, 2):
                    if (factor // prime_divs) % i:
                        continue
                    if (nimber_to_primes**i).val == 1:
                        return i * prime_divs
                return 1 << (1 << L) - 1
                # # make more efficient by only checking possible orders
                # # use Lagrange's theorem
                # # find the smallest field containing n i.e. smallest F_k > n
                # exp = (n.bit_length() - 1).bit_length()
                # # find the order of n must divide F_k - 1 which factors by difference of squares
                # divisors = fermat_divisors(exp)
                # for factor in divisors[:-1]:
                #     if (Nim(n)**factor).val == 1:
                #         return factor
                # else:
                #     return divisors[-1]
            # else:
            #     for factor in fermat_divisors(5):
            #         if (Nim(n)**factor).val == 1:
            #             return factor
            #         # brute force: will probably loop forever
            #         i = 1 << (1 << 5) + 1
            #         while True:
            #             if (Nim(i)**factor).val == 1:
            #                 return i
            #             i += 2
        else:
            ...  # inifite case is hard..

    def sum(self, *args):
        if len(args) == 0:
            return self
        elif len(args) == 1:
            return self + args[0]
        else:
            *inf, fin = ord_decomp(self.val)
            terms = {term.exponent: term.coefficient for term in inf}
            terms[0] = fin
            for arg in args:
                *inf_, fin_ = ord_decomp(arg.val)
                terms[0] ^= fin_
                for term_ in inf_:
                    try:
                        terms[term_.exponent] = (
                            terms[term_.exponent] ^ term_.coefficient
                        )
                    except:
                        terms[term_.exponent] = term_.coefficient
            result = terms[0]
            for key in sorted(terms.keys())[1:]:
                if terms[key] > 0:
                    result = Ordinal(key, terms[key]) + result
            return Nim(result)

    def sq(self):
        """returns the Nimber's square using Freshman's Dream"""
        if self.isfinite:
            a, b, base = self.high, self.low, self.base
            # x = a *2^2^n + b
            # x^2 = (a^2)*(2^2^n + 2^(2^n-1)) + b^2
            if self.base == 0:  # either 0 or 1
                return self
            else:
                term = base + (base >> 1)
                return Nim(a).sq() * Nim(term) + Nim(b).sq()
        else:
            terms = ord_decomp(self.val)
            sum = Nim(terms[-1]).sq()
            for term in terms[:-1]:
                sum = sum + Nim(term) * Nim(term)
            return sum

    def sqrt(self):
        if self.isfinite:
            if self.field == 2:
                return self
            term = self.sq() + self
            return term.sqrt() + self
        else:
            ...  # not sure how to implement...
