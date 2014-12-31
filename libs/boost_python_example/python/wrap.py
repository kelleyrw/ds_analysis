import sys
import numpy as np;
sys.path.append("/Users/rwk7t/Development/ds_analysis/bin/release")

import example
z1 = np.zeros((5,6), dtype=float)
print(z1)
example.fill1(z1)
print(z1)
example.fill2(z1)
print(z1)
