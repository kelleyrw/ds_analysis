import dict_example as de
print(dir(de))

# my_dict = {'pythondict': 1, 'foo': 'bar', 'Bach': 'Goldberg Variation' }
# de.pass_dict(my_dict)

import numpy as np
a1 = np.array([1,2,3])
a2 = np.array([4,5,6,7])
a3 = np.array([[9,10,11,13],[14,15,16,17]])
d = {'one': a1, 'two' : a2, 'three': a3}
de.pass_dict_np(d)
