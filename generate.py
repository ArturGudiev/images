from __future__ import print_function
from functools import reduce
from numpy import sin, pi
import numpy as np
import sys

N = 1
p1 = [-N, N]
p2 = [N, -N]
step = 1

n1 = reduce((lambda x, y: step*abs(x - y) + 1), p1)
n2 = reduce((lambda x, y: step*abs(x - y) + 1), p2)
N = n1 * n2
nums = np.zeros((n1, n2))
i = 0
num = 0
f = open('workfile', 'w')
sys.stdout = f
print(N)
for y in range(p1[1], p2[1] - 1, -step):
    for x in range(p1[0], p2[0] + 1, step):
        print(x, y)
        nums[i / n1, i % n1] = i
        i += 1

for y in range(p1[1], p2[1]-1, -1):
    l = []
    for x in range(p1[0], p2[0]+1):
        l.append(sin(pi*x)*sin(pi*y))
        #l.append(-1 if x < 0 else (0 if x == 0 else 1))
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