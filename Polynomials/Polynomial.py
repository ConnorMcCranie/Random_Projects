'''
create a class which represents polynomials
'''
import numpy as np





class Polynomial:
    '''
    create a class which represents polynomials of a single variable
    ''' 
    def __init__(self, coeffs = [], name : str = 'p'):

        def reduced_coeffs(coeffs : list) -> list:
            '''
            remove trailing zeros from a list of coefficients
            '''
            reduced_coeffs = coeffs
            while len(reduced_coeffs)>1:
                if reduced_coeffs[-1] == 0.0:
                    reduced_coeffs.pop()
                else:
                    return reduced_coeffs
            return [0]
        
        self.coeffs = reduced_coeffs(coeffs)
    
        def reduced_degree(coeffs : list) -> int:
            '''
            return the degree of a polynomial given its reduced coefficients
            '''
            if len(coeffs) == 1 and coeffs[0] == 0.0:
                return -1*np.inf
            else:
                return len(reduced_coeffs(coeffs)) - 1

        self.deg = reduced_degree(self.coeffs)
        fancy_name = f'{name}(x)'
        self.name = name
        self.full_name = fancy_name

        def write_poly(coeffs : list) -> str:
            '''
            write a polynomial as a string
            '''
            # check if the polynomial is the zero polynomial
            n = len(self.coeffs)
            if n == 1 and self.coeffs[0] == 0.0:
                return f'{name} = 0'
            
            # if not, find the first non-zero coefficient
            non_zero : int = 0
            while self.coeffs[non_zero] == 0.0:
                non_zero += 1

            def omit_ones(a : int) -> str:
                '''
                omit the coefficient if it is 1
                '''
                if a == 1.0:
                    return ''
                elif a == -1.0:
                    return '-'
                else:
                    return str(a)

            # write x^0 as a constant term and x^1 as x
            if non_zero == 0:
                poly_str = f'{self.full_name} = {coeffs[non_zero]}'
            elif non_zero == 1:
                poly_str = f'{self.full_name} = {omit_ones(coeffs[non_zero])}x'
            else:
                poly_str = f'{self.full_name} = {omit_ones(coeffs[non_zero])}x^{non_zero}'
            
            # write the rest of the polynomial
            for i in range(non_zero + 1, n):
                if coeffs[i] != 0.0:
                    if i == 1:
                        poly_str += f' + {omit_ones(coeffs[i])}x'
                    else:
                        poly_str += f' + {omit_ones(coeffs[i])}x^{i}'
            return poly_str
        
        self.poly_str = write_poly(self.coeffs)
        
    
   
    
def poly_add(p : Polynomial, q : Polynomial ) -> Polynomial:
    '''
    add two polynomials
    '''
    n = max(p.deg, q.deg)
    p_coeffs = p.coeffs + [0]*(n - p.deg)
    q_coeffs = q.coeffs + [0]*(n - q.deg)
    sum_coeffs = [p_coeffs[i] + q_coeffs[i] for i in range(n+1)]
    return Polynomial(sum_coeffs, f'({p.name}+{q.name})' )

def poly_mult(p : Polynomial, q : Polynomial ) -> Polynomial:
    '''
    multiply two polynomials
    '''
    n = p.deg
    m = q.deg
    prod_coeffs = [0]*(n+m+1)
    for i in range(n+1):
        for j in range(m+1):
            prod_coeffs[i+j] += p.coeffs[i]*q.coeffs[j]
    return Polynomial(prod_coeffs, f'({p.name}{q.name})')

# debugging
'''
p = Polynomial([1, 1.0, 1,0.0], 'p')
q = Polynomial([0.0,-1, 1, 0,0,0,0,0],'q')


    
r = poly_add(p, q)
print(f'The sum of {p.poly_str} and {q.poly_str} is {r.poly_str} and their product is {poly_mult(p, q).poly_str}')
'''