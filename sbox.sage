from sage.crypto.boolean_function import BooleanFunction

# define GF(2^8)
R.<x> = PolynomialRing(GF(2))
f = x^8 + x^7 + x^6 + x^5 + x^4 + x^2 + 1
F.<a> = GF(2^8, modulus=f)

BF = GF(2) 

def byte_to_poly(byte):
    coeffs = [(byte >> i) & 1 for i in range(8)]
    return F(coeffs)

def poly_to_byte(poly):
    coeffs = poly.polynomial().coefficients(sparse=False)
    return sum(int(c) << i for i, c in enumerate(coeffs))

A1=matrix(BF,[
   [1,0,1,0,0,1,1,1],
   [0,1,0,0,1,1,1,1],
   [1,0,0,1,1,1,1,0],
   [0,0,1,1,1,1,0,1],
   [0,1,1,1,1,0,1,0],
   [1,1,1,1,0,1,0,0],
   [1,1,1,0,1,0,0,1],
   [1,1,0,1,0,0,1,1]
 ])
 
A2=matrix(BF,[
  [1,1,0,0,1,0,1,1],
  [1,0,0,1,0,1,1,1],
  [0,0,1,0,1,1,1,1],
  [0,1,0,1,1,1,1,0],
  [1,0,1,1,1,1,0,0],
  [0,1,1,1,1,0,0,1],
  [1,1,1,1,0,0,1,0],
  [1,1,1,0,0,1,0,1]
 ])

C1 = vector(BF, [1, 1, 0, 0, 1, 0, 1, 1])  
C2 = vector(BF, [1, 1, 0, 1, 0, 0, 1, 1])

# sbox constants
S_M     = 8                      #input bit length
S_N     = 8                      #output bit length
S_SIZE = 2 ^ S_M          #the number of elements

###########################################################
def sm4_sbox(byte):
    v = vector(BF, [(byte >> i) & 1 for i in range(S_M-1, -1, -1)])

    v1 = A1 * v + C1
    r_byte = sum(int(v1[i]) << i for i in range(S_M-1, -1, -1))
    elem = byte_to_poly(r_byte)

    if elem != 0:
       inv_elem = elem^-1
       inv = poly_to_byte(inv_elem)
    else:
       inv = 0

    v2 = vector(BF,[(inv >> i)  & 1 for i in range(S_M)])  
    r = A2 * v2 + C2
    return sum(int(r[i]) << (S_M-1-i) for i in range(S_M))

sm4_sbox_table = [sm4_sbox(i) for i in range(S_SIZE)]

print(f"const uint8_t sm4_sbox[{S_SIZE}] = {{")
for i in range(0, S_SIZE, 16):
    row = ", ".join(f"0x{s:02X}" for s in sm4_sbox_table[i:i+16])
    print("    " + row + ",")
print("};", end="\n\n")

###########################################################
def sbox_balance(S):
    for j in range(S_N):
       cnt = sum((S[i]>>j)&1 for i in range(S_SIZE))
       print(f"Output bit {j}: {cnt} ones")

sbox_balance(sm4_sbox_table);  print()

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

sbox_boolfun_property(sm4_sbox_table);  print("")

###########################################################
def sbox_differential_uniformity(S):
    ddt = {}
    for dx in range(1, S_SIZE):
        for x in range(S_SIZE):
            dy = S[x] ^^ S[x^^dx]
            ddt[(dx, dy)] = ddt.get((dx, dy), 0) + 1

    delta = max(ddt.values())
    print("Differential Uniformity =", delta)

sbox_differential_uniformity(sm4_sbox_table);  print()

###########################################################
def sbox_fixed_points(S):
    fps = []
    for x in range(S_SIZE):
        if S[x] == x:
            fps.append(hex(x))
    return fps

fps = sbox_fixed_points(sm4_sbox_table)
print(f"has {len(fps)} fixed points: {fps}", end="\n\n")

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

print(sbox_sac(sm4_sbox_table));  print()

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

for k in range(1, 4):
    sbox_check_pck(sm4_sbox_table, k);  print()
