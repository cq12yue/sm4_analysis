load("sbox.sage")

from sage.crypto.sbox import SBox

# Note: this method is rather slow
def sbox_lat_make(S):
    lat = {}
    D = 2^(S_M-1)   

    for a in range(S_SIZE):
       v_a = V([(a >> i) & 0x1 for i in range(S_M-1, -1, -1)])
       for b in range(S_SIZE):
            v_b = V([(b >> i) & 0x1 for i in range(S_M-1, -1, -1)])
            for x in range(S_SIZE):
                v_x = V([(x >> i) & 0x1 for i in range(S_M-1, -1, -1)])
                y = S[x]
                v_y = V([(y >> i) & 0x1 for i in range(S_N-1, -1, -1)])
                if v_a.dot_product(v_x) == v_b.dot_product(v_y):
                    lat[(a, b)] = lat.get((a,b), 0) + 1
            lat[(a,b)] -= D
    
    return lat

# LAT = sbox_lat_make(sbox_table)
# sbox_dict_print(LAT)

###########################################################
def sbox_best_lat_item(LAT):
    max_val = 0
    best_a, best_b = None, None

    for a in range(1, 256):
        for b in range(1, 256):
            val = abs(LAT[a,b])
            if val > max_val:
               max_val = val
               best_a, best_b = a, b

    print(f"The best lat item is [0x{best_a:02x},0x{best_b:02x}], |bias| is {max_val/256}")

S = SBox(sbox_table)
LAT = S.linear_approximation_table()
print(LAT); print() 

sbox_best_lat_item(LAT)
