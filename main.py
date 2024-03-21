from build.polynomial import Polynomial as poly
import random
a, b = [], []
for i in range(100000):
    a.append(random.uniform(-1000,1000))
    b.append(random.uniform(-1000,1000))
a = poly(a)
b = poly(b)
import time
t1 = time.time()
print((a*b)[1])
print(time.time()-t1)
print(a[1]*b[1])