''' First, define some funtions for integer to bit conversion '''

def numtobits(n : int) -> list[int]:
    # given n, return a list of which bits are 1 in the binary
    index = [i for i in range(n.bit_length()) if (n >> i) & 1]
    return index
def bitstonum(bits : list[int]) -> int:
    # given array of which bits are 1 in binary, return the integer
    return sum([1 << bit for bit in bits])

def subfaces(face : int, nonempty : bool = True) -> set[int]:
        # given face, return all subfaces
        pts = numtobits(face)
        start = 1 if nonempty else 0
        subs = [numtobits(sub) for sub in range(start, 1 << len(pts))]
        subface = [[pts[i]   for i in sub] for sub in subs]
        return {bitstonum(bits) for bits in subface}  

class Simp:
    ''' Encodes an abstract simplicial complex using bitwise operations 
    on integers. For example, the line {{x0}, {x1}, {x0, x1}} would be 
    {1,10,11} in binary. Only need to feed in the maximal simplices, all
    others can be infered. Right now only works for simplices with less
    than 65 vertices. '''

    
    def __init__(self, max_faces : int | list[int]) -> None:
        # determine which points occur in the complex
        total = 0
        for face in max_faces:
             total |= face
        self.points = {1 << bit for bit in numtobits(total)}
        max_faces = [max_faces] if type(max_faces) == int else max_faces


        self.dim : int = max([face.bit_count() for face in max_faces]) - 1
        self.maxs = max_faces
        self.faces : dict[set[int]] | None = None # compute only if necessary
        
          
    def getFaces(self, ordered : bool = False) -> dict[set | list]: 
        # compute all faces from maximal ones
        if self.faces: 
            if not ordered:
                return self.faces
            else: 
                return {n: sorted(list(self.faces[n])) for n in self.faces}
        ''' a bit faster if these funcs are in local scope'''
        def numtobits(n : int) -> list[int]:
            # given n, return a list of which bits are 1 in the binary
            index = [i for i in range(n.bit_length()) if (n >> i) & 1]
            return index
        def bitstonum(bits : list[int]) -> int:
            # given array of which bits are 1 in binary, return the integer
            return sum([1 << bit for bit in bits])
        def subfaces(face : int, nonempty : bool = True) -> set[int]:
                # given face, return all subfaces
                pts = numtobits(face)
                start = 1 if nonempty else 0
                subs = [numtobits(sub) for sub in range(start, 1 << len(pts))]
                subface = [[pts[i]   for i in sub] for sub in subs]
                return {bitstonum(bits) for bits in subface}  
        faces : list[set] = [subfaces(max_face) for max_face in self.maxs]
        self.faces : dict = {-1: set(), 0: self.points} # {} is -1 simplex :)
        for n in range(1, self.dim + 1):
             self.faces[n] = {k for k in set.union(*faces) 
                              if k.bit_count()==n+1}
        if not ordered:
            return self.faces # {n: {simplices of dimension n}}
        # else, {n: [simplices of dimension n in increasing order]}
        return {n: sorted(list(self.faces[n])) for n in self.faces}
            
    def bd(self, simplex : int, orient : bool = False) -> set[int]:
         ''' Returns the maximal subfaces of the simplex. For now I don't
         care about doing the orientation'''

         bits = [i for i in range(simplex.bit_length()) if (simplex >> i & 1)]
         if not orient:
            return {simplex - (1 << bit) for bit in bits}
         
         else: # if oriented, the orientation is encoded by negative signs
              return {(-1) ** i * (simplex - (1 << bit)) 
                      for i, bit in enumerate(bits)}

    def complex(self, reduced : bool = True,
                orient : bool = False) -> dict[np.ndarray]:
         ''' Returns a a dict {n: d_n}, where d_n is the nth boundary map,
         expressed as a matrix d_n : F_2{n-dim faces} -> F_2{(n-1)-dim faces}

         If reduced is True, then the empty set is considered a (-1)-dim face,
         in which case d_0 = [1, 1, ..., 1]. Otherwise, d_0 = [0, 0, ..., 0]
         '''
         num_pts = len(self.points)
         faces = self.getFaces(ordered=True)
         bdry = {0: np.full((1,num_pts), reduced, dtype=bool)}
         if not orient: # +1 = -1 %2  , so we can just do True and False
            for n in range(1, self.dim + 1):
                 simplex, face = faces[n], faces[n-1]
                 dom, rang = len(simplex), len(face)
                 bdry[n] = np.array([[face[j] in self.bd(simplex[i]) 
                                for i in range(dom)] 
                                for j in range(rang)], dtype=bool)
         else: # need to keep track of negative signs
              for n in range(1, self.dim + 1):
                   pass
         return bdry
                    
            