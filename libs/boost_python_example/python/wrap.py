import sys
sys.path.append("/Users/rkelley/Development/rovi/ds_analysis/bin/release")

import numpy as np;
import example

z1 = np.zeros((5,6), dtype=float)
print(z1)
example.fill1(z1)
print(z1)
example.fill2(z1)
print(z1)
