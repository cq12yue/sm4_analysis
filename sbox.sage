from sage.crypto.boolean_function import BooleanFunction
from sage.modules.free_module import VectorSpace

V = VectorSpace(GF(2), 8)

R.<x> = PolynomialRing(GF(2))
f = x^8 + x^7 + x^6 + x^5 + x^4 + x^2 + 1
F.<a> = GF(2^8, modulus=f)

def byte_to_poly(byte):
    coeffs = [(byte >> i) & 1 for i in range(8)]
    return F(coeffs)

def poly_to_byte(poly):
    coeffs = poly.polynomial().coefficients(sparse=False)
    return sum(int(c) << i for i, c in enumerate(coeffs))

M = MatrixSpace(GF(2), 8, 8)

A1 = M ([
    [1,0,1,0,0,1,1,1],
    [0,1,0,0,1,1,1,1],
    [1,0,0,1,1,1,1,0],
    [0,0,1,1,1,1,0,1],
    [0,1,1,1,1,0,1,0],
    [1,1,1,1,0,1,0,0],
    [1,1,1,0,1,0,0,1],
    [1,1,0,1,0,0,1,1]
   ])
 
A2 = M([
   [1,1,0,0,1,0,1,1],
   [1,0,0,1,0,1,1,1],
   [0,0,1,0,1,1,1,1],
   [0,1,0,1,1,1,1,0],
   [1,0,1,1,1,1,0,0],
   [0,1,1,1,1,0,0,1],
   [1,1,1,1,0,0,1,0],
   [1,1,1,0,0,1,0,1]
  ])

C1 = V([1, 1, 0, 0, 1, 0, 1, 1])  
C2 = V([1, 1, 0, 1, 0, 0, 1, 1])

# sbox constants
S_M     = 8                      #input bit length
S_N     = 8                      #output bit length
S_SIZE = 2 ^ S_M          #the number of elements

###########################################################
def sm4_sbox(byte):
    v = V([(byte >> i) & 1 for i in range(S_M-1, -1, -1)])

    v1 = A1 * v + C1
    r_byte = sum(int(v1[i]) << i for i in range(S_M-1, -1, -1))
    elem = byte_to_poly(r_byte)

    if elem != 0:
       inv_elem = elem^-1
       inv = poly_to_byte(inv_elem)
    else:
       inv = 0

    v2 = V([(inv >> i)  & 1 for i in range(S_M)])  
    r = A2 * v2 + C2
    return sum(int(r[i]) << (S_M-1-i) for i in range(S_M))

def sm4_sbox_create():
    table = [sm4_sbox(i) for i in range(S_SIZE)]
    return table

def sm4_sbox_print(S):
    print(f"const uint8_t sm4_sbox[{S_SIZE}] = {{")

    for i in range(0, S_SIZE, 16):
        row = ", ".join(f"0x{s:02X}" for s in S[i:i+16])
        print("    " + row + ",")

    print("};", end="\n\n")

###########################################################
def sbox_balance(S):
    for j in range(S_N):
       cnt = sum((S[i]>>j)&1 for i in range(S_SIZE))
       print(f"Output bit {j}: {cnt} ones")

###########################################################
def max_nonlinearity(n):
    if n % 2 == 0:
        return 2^(n-1) - 2^(n//2 - 1)
    else:
        return 2^(n-1) - 2^((n-1)//2)

def sbox_boolfun_property(S):
    min_nl = infinity

    for j in range(S_N):
        bf = BooleanFunction([(S[i]>>j)&1 for i in range(S_SIZE)])
        deg = bf.algebraic_degree()
        nl = bf.nonlinearity()
        min_nl = min(min_nl, nl)
        walsh_max = max(abs(w) for w in bf.walsh_hadamard_transform())
        print(f"Bit {j}: degree={deg}, nonlinearity={nl}, max|Walsh|={walsh_max}")

    print(f"the minimum nonlinearity is {min_nl}, theory max nonlinearity is {max_nonlinearity(S_N)}")

###########################################################
def sbox_fixed_points(S):
    fps = []
    for x in range(S_SIZE):
        if S[x] == x:
            fps.append(hex(x))
    return fps

###########################################################
def flip_bit(x, i):
    return x ^^ (1 << i)

def sbox_sac(S):
    sac_matrix = matrix(QQ, S_M, S_N, 0)
    for x in range(S_SIZE):
        for i in range(S_M):          
            xp = flip_bit(x, i)
            dx = S[xp] ^^ S[x]      
        
            for j in range(S_N):     
                if (dx >> j) & 1:
                    sac_matrix[i, j] += 1

    # Normalize to probability
    sac_matrix = sac_matrix / S_SIZE
    return sac_matrix

def sbox_check_pck(S, k):
    bool_funcs = []
    bf_satisfy_pcks = []

    for i in range(S_N):
        bf = BooleanFunction([(S[x] >> i) & 1 for x in range(S_SIZE)])
        bool_funcs.append(bf)
        bf_satisfy_pcks.append(True)

    for i, f in enumerate(bool_funcs):
        w = f.walsh_hadamard_transform()

        for a in range(1,  S_SIZE):
            if bin(a).count('1') <= k:
                if w[a] != 0:
                    # D = f.derivative(a)             
                    # if not D.is_balanced():
                    bf_satisfy_pcks[i] = False
                    break               
        r = bf_satisfy_pcks[i]
        print(f"bf[{i}] satify PC({k}):{r}")

###########################################################
def sbox_dict_print(Dict):
    col_h = " ".join(f"0x{b:02x}" for b in range(S_SIZE))
    print(" ", end="");  print(col_h)
   
    for a in range(S_SIZE):
        row = f"0x{a:02x} "
        for b in range(S_SIZE):
            if (a,b) in Dict:
                row += str(Dict[(a,b)])
            else:
                row += "-"
            row += " "
        print(row)

###########################################################

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--out_basic", action="store_true")
parser.add_argument("--out_sbox", action="store_true")
args = parser.parse_args()

sbox_table = sm4_sbox_create()

if args.out_sbox:
    sm4_sbox_print(sbox_table)

if args.out_basic:
    sbox_balance(sbox_table)
    print()
   
    sbox_boolfun_property(sbox_table)
    print()
    
    fps = sbox_fixed_points(sbox_table)
    print(f"has {len(fps)} fixed points: {fps}", end="\n\n")

    print(sbox_sac(sbox_table))
    print()

    for k in range(1, 4):
        sbox_check_pck(sbox_table, k);  print()
