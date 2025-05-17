class Top:
    ''' A class for representing a finite topological space. 
        ------------
        Attributes |
        ----------------------------------------------------------------
        points: a set {x | x in X } of the points of the top. space
                The elements should all be hashable

        opens: a dictionary of the basic open sets. It is of the form
        { x : { y | x <= y} } where x <= y means x in closure({y})

        closures:   dictionary of the closures of each point
                    { x : cl({x}) }
                    Starts as empty because computing all is O(n^2)
                    Run getClosures() method if all are needed. Else,
                    just call X.cl(x) or X[x] for closure of 1 point
        ----------------------------------------------------------------
    
    '''
    def __init__(self, data, discrete=False, boolean=False):
        # can pass a networkx graph as the data
        # import networkx as nx

        if type(data) == dict: # dict of open sets {x: {y | y <= x} }
            self.opens = data
            self.points = set(data)
        
        elif type(data) == int:
            assert data >= 0, 'Topology of integer can\'t be negative'
            assert (discrete and boolean)==False, ('boolean lattice is '
                                                    'not discrete')
            self.points = {n for n in range(data)}            
            if not discrete: 
                if not boolean: # order linearly
                    self.opens = {n: {k for  k in range(n+1)} 
                                for n in range(data)}
                else: # make power set on n elements <= under inclusion
                    # use binary expansion to encode subsets
                    pts = range(1 << data)
                    self.points = {pt for pt in pts}
                    self.opens = {pt: {x for x in pts if x & pt == x}
                                  for pt in pts}

            else: # discrete topology
                self.opens = {n: {n} for n in range(data)}

        ''' optional functionality with networkx directed graphs, but I
        currently am not using this, so I don't want to bother importing '''

        # elif type(data) == nx.DiGraph:
        #     opens = {}
        #     for node in data.nodes:
        #         opens[node] = nx.descendants(data,node) | {node}
        #     self.opens = opens
        #     self.points = set(data.nodes)
        ''' The following attributes won't be filled in until computed 
        because they can be expensive to compute for large # of points'''
        self.closures = {}      # can specify closures and ordering,
                                # call self.getClosures()
        self.ordering = None    # but left empty until needed
                                # call self.order()
        self.maxs = None    # max/minimal elements w.r.t. partial ordering
                            # call self.getMaxs()
        self.mins = None    # (must be T_0 space to make sense)
                            # call self.getMins()
        self.T0 = None      # has been checked to be T0
                            # call self.isT0()
        self.beats = None   # points with a unique succesor or predecesor
                            # call self.getBeats()


    def __repr__(self):
        return str(self.opens)
    
    def __len__(self):
        return len(self.points)
    
    def __contains__(self, x):
        return x in self.points
    
    def __iter__(self):
        return iter(self.points)

    def __call__(self, x):
        # Top(x) returns the same as self.opens[x] = {y | y <= x}
        return self.opens[x]
        
    def __getitem__(self, x):
        # Top[x] returns the same as Top.cl(x) or Top.closures[x]
        # This convention mirrors the notation (a,b) for open interval 
        # and [a,b] for closed interval
        try:
            return self.closures[x]
        except:
            self.closures[x] = self.cl(x)
            return self.closures[x]
        
    def __mul__(self, other):
        # product topology
        return self.product(other)
    
    def __add__(self, other):
        # disjoint union topology
        if not self.points & other.points: # already disjoint
            return Top(self.opens | other.opens)
        else: # force them to be disjoint
            point0 = Top({0: {0}})
            point1 = Top({1: {1}})
            return Top((point0 * self).opens | (point1 * other).opens)
        
    def __eq__(self, other):
        return self.opens == other.opens
    
    def __le__(self, other):
        # is a subspace
        if self.points & other.points != self.points:
            return False # isn't a subset
        for point in self.points:
            if other(point) & self.points != self(point):
                return False # not the same open sets
        return True
    
    def __lt__(self, other):
        # proper subspace
        return (not self == other) and self <= other
    
    def __ge__(self, other):
        # contains other as a subspace
        return other <= self
    
    def __gt__(self, other):
        # proper super space
        return (not self == other) and self >= other
    
    def __truediv__(self, subspace):
        # quotient by a subspace
        assert self >= subspace, 'Quotient defined for a subspace'
        quotient = set(self.points - subspace.points)

        # the quotient projection
        def pi(x, pt): 
            if x in self.points - subspace.points:
                return x
            elif x in subspace.points:
                return pt
            else:
                raise ValueError('projection map domain error')
        
        # adjoin a new point not in X \ A, make sure key isn't taken
        pt = '*'; idx = 0
        while pt in quotient:
            idx += 1
            pt = '*' + str(idx)
        # now make the open sets for each x in X \ A
        opens = {x: {pi(y, pt) for y in self(x)} for x in quotient}
        # add in the open set around the new extra point 
        opens[pt] = set()
        for a in subspace.points:
            opens[pt] |= {pi(x, pt) for x in self(a)}
        return Top(opens)
    
    def __pow__(self, other):
        return other.hom(self)

    def __invert__(self):
        return self.op()

    def cl(self, x):
        # returns the closure of the singleton {x}
        assert x in self.points,\
        f'{x} is not a point in the topological space'
        close = {y for y in self if x in self(y)}
        # close = { y | x <= y } = { y | x in self.opens[y] }
        return close
        
    def getClosures(self):
        # populates Top.closures for all elements
        for x in self:
            try: # don't bother if it's already computed
                self.closures[x]
            except:
                self.closures[x] = self.cl(x)

    def isleq(self, x, y):
        # partial order defined by x<=y iff x in self.opens[y]
        # or equivalently iff y in self.closures[x]
        return x in self(y)
    
    def isT0(self):
        # True iff self.opens[x] = self.opens[y] implies x = y
        if self.T0 == None:
            checked = set()
            for x in self:
                for y in self.points - (checked | {x}):
                    if self(x) == self(y):
                        self.T0 = False
                        return False
                checked |= {x} # don't check U_x=U_y and U_y=U_x separately
            self.T0 = True
        return self.T0
    
    def isT1(self):
        # for finite spaces, this is the same as being discrete
        for x in self:
            if self(x) != {x}:
                return False
        return True
    
    def subspace(self, subset):
        # subspace topology
        assert subset & self.points == subset, 'needs to be a subset'
        return Top({x: self(x) & subset for x in subset})
    
    def product(self, other):
        # cartesian product
        def prod(set1, set2):
            return [(x, y) for x in set1 for y in set2]
        ord1 = self.order()
        ord2 = other.order()
        points = prod(ord1, ord2)
        opens = {point : set(prod(self(point[0]), other(point[1]))) \
                 for point in points}
        return Top(opens)
    
    def op(self):
        # reverse the ordering
        self.getClosures()
        op_opens = self.closures
        X_op = Top(op_opens)
        X_op.closures = self.opens
        X_op.maxs = self.mins
        X_op.mins = self.maxs
        if self.ordering:
            op_order = list(self.ordering).copy()
            op_order.reverse()
            X_op.ordering = tuple(op_order)
        return (X_op)
    
    def order(self, update=True):
        # returns a topological ordering of self.points, meaning that if
        # x <= y, then x will appear first in the list (but converse may 
        # be false if x and y aren't comparable)

        nodes = {x: {'start':0, 'end':0, 'parent':None, 'visited':False}
                 for x in self.points} # nodes for depth first search
        time = 0 
        def dfsVisit(pt): # one iteration of depth first search
            nonlocal time
            time += 1
            node = nodes[pt]
            node['start'] = time
            node['visited'] = True

            for child in sorted(list(self[pt]), key=hash, reverse=True): # pt <= child
                childnode = nodes[child]
                if childnode['visited'] == False:
                    childnode['parent'] = node
                    dfsVisit(child)
            time += 1
            node['end'] = time
            node['visited'] = 'done'

        # now iterate over all points 
        ordered_pts = sorted(list(self.points), key=hash, reverse=True)
        for pt in ordered_pts:
            if nodes[pt]['visited'] == False:
                dfsVisit(pt)

        # after depth first search, order the points by finish time
        top_order = sorted(nodes, reverse=True,
                           key=lambda x: nodes[x]['end'] )
        if update:
            self.ordering = tuple(top_order)
        return tuple(top_order)
    
    def to(self, other, inj=False):
        ''' returns a list of continuous fnunctions self -> other. Each
        function is represented as a dictionary {x: f(x)}.
        X.to(Y) returns [{x: f(x)} for f: X --> Y continuous]'''
        P_ordered = self.order()          # topologically sorted points
        results = []

        # helper function to iteratively modify an existing function
        def backtrack(index, current_map):  
            if index == len(P_ordered): # complete function
                results.append(current_map.copy())
                return 
            x = P_ordered[index]
            for q in other.points: # potential values for f(x)
                valid = True # make sure (x <= y) ==> ( f(x) <= f(y) )
                for y in current_map:
                    fy = current_map[y]
                    if y in self.opens[x] and fy not in other.opens[q]:
                        valid = False
                        break
                if valid:
                    current_map[x] = q # assign f(x) = q
                    backtrack(index + 1, current_map) # progress 1 step
                    # reached end of valid mapping
                    del current_map[x] # remove last assignment, iterate

        backtrack(0, {})
        return results
    
    def hom(self, other, domain_order=False):
        # X.hom(Y) returns the topological space Hom(X, Y) with the com-
        # pact-open topology: f <= g iff f(x) <= g(x) for all x
        funcs = self.to(other)
        opens = {}
        order = self.ordering if self.ordering else self.order()
        def funcTuple(func):
            # use top ordering to write func as (f(x_0), ..., f(x_n))
            # so that it is hashable            
            return tuple([func[x] for x in order])

        for f1 in funcs:
            key = funcTuple(f1)
            opens[key] = {funcTuple(f2) for f2 in funcs
                          if all(f2[x] in other.opens[f1[x]] 
                                        for x in self.points)    
                            }
        if domain_order:
            return order, Top(opens)
        return Top(opens)
    
    def relabel(self):
        # returns a homeomorphic top space with integer labels
        order = self.ordering if self.ordering else self.order()
        opens = {n:{m for m in range(len(self)) 
                    if order[m] in self(order[n])}
                for n in range(len(self))}
        return Top(opens)
    
    def succ(self, x, ordered=False):
        ''' Given point x in X, returns the immediate successors 
            {z in X | x < z and (x<y<=z ==> y=z) }  '''
        assert x in self.points, f'{x} is not a point in the topological space'
        succ = self[x] - {x} # remove x from the set of successors
        for y in succ.copy(): # if y > x
            if y not in succ: # y may have been removed proviously
                continue
            for z in (self[y] - {y}) & succ: # and z > y
                succ -= {z} # z can't be a successor of x
        if ordered:
            succ = self.sort(succ)       
        return succ
    
    def susp(self):
        ''' Returns the (discrete) suspension '''
        north, south  = 'N', 'S' # poles to attach lines to 
        if north in self:
            N_num = 1
            while north + str(N_num) in self: 
                N_num += 1
            north += str(N_num)
        if south in self:
            S_num = 1 
            while south + str(S_num) in self:
                S_num +=1        
            south += str(S_num)
        # unique str for poles found   
         
        opens = (self.opens).copy() # original space embeds into suspension
        opens[north] = self.points | {north}
        opens[south] = self.points | {south} # both poles maximal
        return Top(opens)
        
    def pred(self, x, ordered=False):
        ''' Given point x in X, returns the immediate predecessors 
            {z in X | z < x and (y<z<=x ==> y=z) }  '''
        assert x in self.points, f'{x} is not a point in the topological space'
        pred = self(x) - {x} # remove x from the set of successors
        for y in pred.copy(): # if y > x
            if y not in pred: # y may have been removed proviously
                continue
            for z in (self(y) - {y}) & pred: # and z > y
                pred -= {z} # z can't be a successor of x
        if ordered:
            pred = self.sort(pred)       
        return pred


    def getMaxs(self, checkT0=False):
        '''Returns the set of maximal elements wrt the partial ordering.
        Needs to be T_0 space otherwise can have infinite chains. '''
        if checkT0: # O(n^2) to check, so better to avoid if possible            
            assert self.isT0(), ('maximal elements aren\'t defined for pre-'
            'orders which aren\'t partial orders') 
        if self.maxs == None: # not already computed
            self.maxs = {x for x in self if self[x] == {x}}
        return self.maxs
    
    def getMins(self, checkT0=False):
        '''Returns the set of maximal elements wrt the partial ordering.
        Needs to be T_0 space otherwise can have infinite chains. '''
        if checkT0: # O(n^2) to check, so better to avoid if possible            
            assert self.isT0(), ('maximal elements aren\'t defined for pre-'
            'orders which aren\'t partial orders') 
        if self.mins == None: # not already computed
            self.mins = {x for x in self if self(x) == {x}}
        return self.mins
    
    def isBeat(self, pt) -> bool:
        if len(self.succ(pt)) == 1: # upbeat
            return True
        return  len(self.pred(pt)) == 1 # true if downbeat else false
    
    def getBeats(self) -> set:
        if self.beats: # already computed
            return self.beats
        beats = set()
        for pt in self.points:
            if len(self.succ(pt)) == 1: beats |= {pt} # unique succesor
            elif len(self.pred(pt)) == 1: beats |= {pt} # unique predecesor
        self.beats = beats
        return self.beats

    def maxChains(self, checkT0=False):
        ''' Returns a list of the maximal chains in a T0 top. sp. '''
        if checkT0:
            assert self.isT0(), ('maximal chains aren\'t defined for pre-'
            'orders which aren\'t partial orders')
        
        order = self.order() # for consistency of ordering
        def index(pt): # pass as key to sorted() to sort by top ordering
            for i, x in enumerate(order):
                if x == pt:
                    return i
        maxs = sorted(self.getMaxs(), key=index)
        mins = sorted(self.getMins(), key=index) 
        chains = []
        

        def extendChain(chain):
            # given a non-empty chain (as list), extend it further is possible
            curr_max = chain[-1]
            bigger = sorted(self.succ(curr_max), key=index)
            if curr_max in maxs: # already maximal, add to list
                chains.append(tuple(chain))
                return
            for pt in bigger:
                new_chain = chain.copy()
                new_chain.append(pt)
                extendChain(new_chain)
        
        # extend all chains starting at minimal elements
        for min in mins:
            extendChain([min])
        return chains # will be ordered lexographically
        
    def grading(self, ordered=False):
        ''' Returns a dictionary {n : {points n levels up from bottom}}
        grading[0]={minimals}, grading[1]={x | y < x ==> y minimal}, etc
        
        If (X, <=) comes from a simplicial complex wrt inclusion, then
        X.grading()[n] is the set of n-faces. If not, then the grading 
        isn't necessarily unique, so '''
        grade = {}
        level = 0
        grade[0] = self.getMins()
        remaining = self.points - grade[0]
        while remaining != set():
            level += 1
            grade[level] = set()
            for pt in grade[level - 1]:
                grade[level] |= self.succ(pt)
            remaining -= grade[level]
        if ordered:
            order = self.order()
            def index(pt): # pass as key to sorted() to sort by top ordering
                for i, x in enumerate(order):
                    if x == pt:
                        return i
            for n in grade: # return ordered tuples instead of sets
                in_order = sorted(list(grade[n]), key=index)
                grade[n] = tuple(in_order)
        return grade

    def index(self, pt): # pass as key to sorted() to sort by top ordering
        order = self.ordering if self.ordering else self.order()
        for i, x in enumerate(order):
            if x == pt:
                return i
    
    def sort(self, subset, check_subset=True) -> tuple:
        ''' given a subset A <= X, returns A in sorted order accoring to 
        the topological ordering of X. '''
        if check_subset:
            assert subset & self.points == subset, 'needs to be a subset'
        pts = list(subset)
        return tuple(sorted(pts, key=self.index))

    def T0closure(self):
        t0 = self.T0 if self.T0 != None else self.isT0()
        if t0: return self # already T0

        # self / ~    to do!

    def core(self, relabel=False): 
        '''Returns a Top space in the same homotopy class with the min
        number of points. Removing a beat point is a retraction'''
        order = self.ordering if self.ordering else self.order()
        core = self
        while True:
            order = core.ordering
            # keep removing beat points until none remain
            for pt in reversed(order): # keep smallest to avoid relabeling
                if core.isBeat(pt):
                    core = core.subspace(core.points - {pt})
                    core.ordering = list(order).remove(pt)
                    continue
            break
        return core.relabel() if relabel else core
