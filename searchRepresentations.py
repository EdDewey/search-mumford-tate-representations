from __future__ import division
from sage.all import *

'''
This file contains code to do two things:
    1) For any n, list all irreducible representations of semisimple
    Lie algebras, with dimension n, that admit a symmetric invariant form
    (this is the function symmetric_representations).
    Representations are listed by giving a highest weight.
    2) For any highest weight representation and any list of Hodge numbers, say
    whether this could be a Mumford-Tate representation with those Hodge
    numbers.  I am not certain whether what it checks is a sufficient
    condition, but it is necessary.  This is the function search_for_circle

    For example, here is how we verified that the Mumford-Tate group of the
    primitive cohomology of an abelian fivefold has semisimple type:
    1) run symmetric_representations(72)
    2) run search_for_circle on each representation we found to exclude all but
    D36

    Tested in sagemath version 
'''
import itertools
import copy

############### list representations
@parallel
def irreps_one_system(cartan_type, n): 
    # output all irreducible characters of the weyl character ring L with
    # degree at most n.
    space = RootSystem(cartan_type).ambient_space()
    cartan_type = space.cartan_type() # avoid formatting problems
    funW = [f for f in space.fundamental_weights() 
            if space.weyl_dimension(f) <= n] # List relelvant fundamental weights
    output = {}
    last_weights = [0]
    while(True):
        weights_to_check = ((e+f, space.weyl_dimension(e+f))  
                             for e in last_weights 
                             for f in funW) 
        valid_weights = [p for p in weights_to_check if p[1] <= n]
        if not valid_weights:
            break
        last_weights = [p[0] for p in valid_weights]
        # Add the new reps we found to the output
        for p in valid_weights:
            dim = p[1]
            if dim in output:
                output[dim] = output[dim] + [p[0]]
            else:
                output[dim] = [p[0]]
    # Clean up
    if cartan_type[0] == 'A':
        for d in output:
            output[d] = [w.coerce_to_sl() for w in output[d]] 
    output = {d:output[d] for d in output if output[d]} # remove empyt list
    if not output.keys(): 
        return {}
    wcr = WeylCharacterRing(cartan_type)
    for d in output:
        maybe0 = list(set(output[d])) # remove redundancies
        maybe1 = [m for m in maybe0 if wcr(m).is_irreducible()] # irreps only
        output[d]=maybe1
    return output
   
def generate_list(n):
    # Create a list of relevant root systems
    typeA_list = ['A'+str(j) for j in range(1,n)]
    big_list = ['F4', 'G2', 'E6', 'E7', 'E8'] + typeA_list
    families = ['B', 'C', 'D']
    start_num = {'A':1, 'B':2, 'C':3, 'D':4}
    for foo in families:
        k = start_num[foo]
        while(True):
            L = RootSystem([foo,k]).ambient_space()
            mindim = min(L.weyl_dimension(f) for f in L.fundamental_weights())
            if mindim > n:
                break
            else:
                big_list = big_list + [foo + str(k)]
                k = k+1
    return big_list

def irreps_of_simple_groups(n, systems = [] ): 
    # List all irreducible representations of all simple lie algebras of
    # dimension at most n.  The output is in the form of a dict {dimension:
    # [highest weight vectors of that dimension]}.  If systems is non-empty,
    # only return representations of lie algebras with root system in systems. 
    output = {}
    if not systems:
        systems = generate_list(n)
    parallel_input = zip(systems, itertools.cycle([n]))
    parallel_output = irreps_one_system(parallel_input)
    for foo in parallel_output:
        irreps_of_L = foo[1]
        for d in irreps_of_L: # d is a dimension
            if d in output:
                output[d] = output[d] + irreps_of_L[d]
            else:
                output[d] = irreps_of_L[d]
    # remove empty entries
    concise_output = {d:output[d] for d in output if output[d]}
    return concise_output


# Now we need to handle semsimple Lie algebras.

def all_factorizations(n):
    # Output all factorizations of the integer n, ignoring singletons
    # (since semisimple Lie groups have no characters)
    if is_prime(n):
        return [(n,)]
    divs = divisors(n)
    output = []
    for d in divs[1:][:-1]: # 1 is always listed as the first divisor
        m = int(n/d)
        d_output = [tuple(sorted((d,) + fact_m)) for fact_m in all_factorizations(m)]
        output = output + d_output
    return list(set(output))+[(n,)]

def tensor_of_weights(ws): 
    # Given a bunch of highest weight representations of simple Lie algebras,
    # output a tensor representation of the product of those Lie algebras
    # I haven't found documentation on coercing to a representation of the
    # product group like we do here, but I tested it and it seems to be doing
    # the right thing.
    cartan_type = [w.parent().cartan_type() for w in ws]
    space = RootSystem(cartan_type).ambient_space()
    vectorize = [vector(w) for w in ws]
    bigTuple = tuple(itertools.chain.from_iterable(vectorize))
    return space(bigTuple)
    
def irreps_factored(factorization, reps):
    # list irresps with a tensor factorization of shape factorization. reps
    # should be the output of irreps_of_simple_groups.
    tensor_factors = [reps[d] for d in factorization]
    possibles = itertools.product(*tensor_factors)
    return [tensor_of_weights(p) for p in possibles]

def irreps(n, reps = []): 
    # list all irreps semisimple lie algebras of dimension equal to n, whose
    # tensor factors are in reps.  reps should be a dictionary of the form { m:
    # highest weight vector of dimension m}.  It is assumed that all
    # representations in reps are of simple lie algebras.
    if not reps:
        reps = irreps_of_simple_groups(n)
    is_useable = lambda f: all([d in reps.keys() for d in f])
    factorizations = (f for f in all_factorizations(n) if is_useable(f))
    output = (w for f in factorizations for w in irreps_factored(f, reps))
    return sorted(list(set(output)))

@parallel
def is_self_dual(w, frobenius_schur_indicator = None):
    # Return True if the representation is self-dual.  Can also require that it
    # admits a symmetric/antisymmetric form by specifying the probenius-schur
    # indicator as 1/-1
    L = w.parent()
    cartan_type = L.cartan_type()
    wcr = WeylCharacterRing(cartan_type)
    true_indicator = wcr(w).frobenius_schur_indicator()
    if frobenius_schur_indicator:
        return true_indicator == frobenius_schur_indicator
    else:
        return true_indicator != 0

def symmetric_representations(n,  sign = 1, irreps_list = [],):
    if not irreps_list:
        irreps_list = irreps(n)
    parallel_input = zip(irreps_list, itertools.cycle([sign]))
    is_self_dual_output = is_self_dual(parallel_input)
    good_reps = [output[0][0][0] for output in is_self_dual_output if output[1]]
    good_reps.sort(key = lambda w: w.parent().rank()) # Just for convenience
    return good_reps

############## Linear algebra

def move_basis_to_end(list_of_vectors):
    # Reorder a list of vectors so that the terminal segment is a basis for the
    # subspace that they generate.  Also output the rank of the space they
    # generate.
    initial_segment = copy.copy(list_if_vectors)
    terminal_segment = [] 
    V = VectorSpace(QQ, len(vectors[-1]))
    target_rank = V.subspace(functionals).dimension()
    current_rank = 0
    while current_rank < target_rank:
        new_vector = initial_segment.pop()
        test_segment = terminal_segment + [new_vector]
        test_rank = V.subspace(test_segment).dimension()
        if test_rank > current_rank:
            terminal_segment = test_segment
            current_rank = test_rank
        else:
            terminal_segment = [new_vector] + terminal_segment
    return (initial_segment + terminal_segment, target_rank)

def coordinate_vectors(vectors):
    # Replace a list of vectors with coordinate vectors with repsect to a basis
    # consisting of elements of the list
    (reordered, rank) = move_basis_to_end(functionals)
    ambient_rank = len(reordered[-1])
    V = VectorSpace(QQ, ambient_rank)
    W = V.subspace_with_basis(reordered[-rank:])
    return [W.coordinate_vector(f) for f in reordered]

@parallel
def check_cocharacter(chi, vectors, outputs_with_mults):
    allowed_outputs = outputs_with_mults.keys()
    current_multiplicities = {key: 0 for key in allowed_outputs}
    for v in vectors:
        t = vector(chi).dot_product(v)
        if t not in allowed_outputs:
            return False
        else:
            current_multiplicities[t] = current_multiplicities[t]+1
    for p in allowed_outputs:
        if current_multiplicities[p] != outputs_with_mults[p]:
            return False
    return True

def find_cocharacter(vectors, outputs_with_mults, just_check=True):
    # Find all vectors chi such that the number of times that the dot product
    # (chi, v) equals p for v in vectors is given by
    # outputs_with_multiplicities[v].  Assumes that vectors includes the
    # standard basis.  just_check controls whether it actually remembers the
    # cocharacters or merely checks whether one exists
    rank = len(vectors[-1])
    allowed_outputs = outputs_with_mults.keys()
    to_check = itertools.product(*rank*[allowed_outputs])
    input = zip(to_check, 
                itertools.cycle([vectors]), 
                itertools.cycle([outputs_with_mults])  )
    output = check_cocharacter(input)
    if just_check:
        return any( foo[1] for foo in output)
    else:
        return (foo[0][0][0] for foo in output if foo[1])

######## Hodge stuff

def get_weights(w):
    # Get the weights of the highest weight representation with highest weight
    # w.  Outputs as a list of vectors.
    cartan_type = w.parent().cartan_type()
    wcr = WeylCharacterRing(cartan_type)
    bar = wcr(w).weight_multiplicities()
    return [vector(u) for u in bar for _ in range(0,bar[u])]

def search_for_circle(w, hodge_numbers, just_check=True):
    # searchs the highest weight representation of w for a cocharacter that
    # could be S^1 inside the Deligne torus.  hodge_numbers should be input as
    # a dict {p, h^{p,q}}.  (This code pays no attention to the weight of the
    # Hodge structure).  just_check controls whether it actually remembers the
    # cocharacters or merely checks whether one exists
    hodge_degrees = hodge_numbers.keys()
    twice_average = max(hodge_degrees)+min(hodge_degrees)
    shifted_hodge_numbers = {2*p-twice_average: hodge_numbers[p] for p in
            hodge_degrees}
    return find_cocharacter(get_weights(w), shifted_hodge_numbers, just_check)
