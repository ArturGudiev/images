from __future__ import print_function
from functools import reduce
from numpy import sin, pi, cos
import numpy as np
import sys

N = 2
p1 = [-N, N]
p2 = [N, -N]
step = 1
n1 = reduce((lambda x, y: int(abs(x - y)/step) + 1), p1)
n2 = reduce((lambda x, y: int(abs(x - y)/step) + 1), p2)
N = n1 * n2
nums = np.zeros((n1, n2))
i = 0
num = 0
f = open('workfile', 'w')
sys.stdout = f
print(N)
for y in np.arange(p1[1], p2[1] - step, -step):
    for x in np.arange(p1[0], p2[0] + step, step):
        print(x, y)
        nums[i / n1, i % n1] = i
        i += 1

func = lambda x, y: -1 if x < 0 else (0 if x == 0 else 1)
# func = lambda x, y: sin(pi*x)*sin(pi*y)

for y in np.arange(p1[1], p2[1]-step, -step):
    l = []
    for x in np.arange(p1[0], p2[0]+step, step):
        l.append(func(x, y))
    print(*l, sep=' ')
    # print(l)

nums = nums.astype(int)

print((n1-1)*(n2-1)*2)
for i in range(0, n1-1):
    for j in range(0, n2-1):
        p1 = nums[i, j]
        p2 = nums[i, j + 1]
        p3 = nums[i+1, j]
        p4 = nums[i+1, j+1]
        print(p1, p2, p3)
        print(p3, p4, p2)