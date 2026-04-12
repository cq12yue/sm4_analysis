load("sbox.sage")

def sbox_ddt_make(S):
     ddt = {}
     for dx in range(1, S_SIZE):
         for x in range(S_SIZE):
             dy = S[x] ^^ S[x^^dx]
             ddt[(dx, dy)] = ddt.get((dx, dy), 0) + 1
     return ddt

def sbox_differential_uniformity(DDT):
    delta = max(DDT.values())
    print("Differential Uniformity =", delta)

###########################################################
DDT = sbox_ddt_make(sbox_table)
sbox_dict_print(DDT); print()

sbox_differential_uniformity(DDT)
