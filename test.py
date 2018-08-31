import numpy as np
import global_functions

a = np.array([1,2,3,4,5])
b =np.flip(a,-1)
print(a[1:])
print(a[:-1])
print(global_functions.calc_dif(b))