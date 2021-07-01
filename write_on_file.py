import numpy as np
a = np.ones((2,2))
b = a/100000000000000
np.savetxt('result.txt', b,'%.2e')