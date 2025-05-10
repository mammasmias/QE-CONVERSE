# diagonalize the g-tensor from QE-CONVERSE outputs, according to the
# prescription of Weyl, Bolton "Electron Paramagnetic Resonance", Chap. 4.4
import numpy as np
import math
import sys

def read_tensor_from_files(file1, file2, file3):
    def find_line_with_keyword(file, keyword):
        with open(file, 'r') as f:
            for line in f:
                if keyword in line:
                    return line.strip()
        raise ValueError(f"Keyword '{keyword}' not found in {file}")

    l1 = find_line_with_keyword(file1, 'delta_g total')
    l2 = find_line_with_keyword(file2, 'delta_g total')
    l3 = find_line_with_keyword(file3, 'delta_g total')

    tens = []
    tens.append(list(map(float, l1.split()[-3:])))
    tens.append(list(map(float, l2.split()[-3:])))
    tens.append(list(map(float, l3.split()[-3:])))
    
    return np.array(tens)

def diag_tensor(g):
    return np.linalg.eigh(np.dot(g, g.T))

def print_tensor(g):
    np.savetxt(sys.stdout, g, fmt="%12.2f")

if len(sys.argv) != 4:
    sys.stderr.write("usage: %s file1.out file2.out file3.out\n" % (sys.argv[0]))
    sys.exit(1)

file1, file2, file3 = sys.argv[1], sys.argv[2], sys.argv[3]

try:
    g_tensor = read_tensor_from_files(file1, file2, file3)
except ValueError as e:
    sys.stderr.write(str(e) + "\n")
    sys.exit(1)

print("Delta g-tensor from files:")
print_tensor(g_tensor)
print()

g_e = 2.00231930436182
print("using g_e =", g_e)
print()

np.set_printoptions(precision=6, suppress=True)

print("Diagonal + principal components (cartesian):")
g_tensor = g_e * np.eye(3) + g_tensor / 1e6
g2, pc = diag_tensor(g_tensor)
print("%12.6f =>" % (math.sqrt(g2[0])), pc[:,0])
print("%12.6f =>" % (math.sqrt(g2[1])), pc[:,1])
print("%12.6f =>" % (math.sqrt(g2[2])), pc[:,2])
print()

