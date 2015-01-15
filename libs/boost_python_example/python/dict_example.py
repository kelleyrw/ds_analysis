#!/usr/bin/env python


import numpy as np
a1 = np.array([1,2,3])
a2 = np.array([4,5,6,7])
a3 = np.array([[9,10,11,13],[14,15,16,17]])
d = {'one': a1, 'two' : a2, 'three': a3}

import dict_example as de
de.pass_dict(d)
de.pass_dict_np(d)
