from __future__ import division
import numpy as np
import sys
from numpy.linalg import inv, solve
from scipy import integrate
from numpy import square as sq
from numpy import linalg as LA
from numpy import linalg
from math import sqrt, sin, cos, pi

def triangleValue(f1, f2, f3):
    return f1 + f2 + f3


def dJ(u, v, p1, p2, p3):
    x1, y1 = p1
    x2, y2 = p2
    x3, y3 = p3
    dxdu = ((1 - v) * x2 + v * x3 - x1)
    dxdv = (u * x3 - u * x2)
    dydu = ((1 - v) * y2 + v * y3 - y1)
    dydv = (u * y3 - u * y2)
    return np.abs(dxdu * dydv - dxdv * dydu)


def tridblquad(integrand, p1, p2, p3):
    '''
    Perform double quadtrature integration on triangular domain.
    Input: function to integrate, points of triangle as tuples.
    Output: integral and estimated absolute error as a tuple.
    '''
    x1, y1 = p1;
    x2, y2 = p2;
    x3, y3 = p3
    # transformation to the unit square
    g = lambda u, v, c1, c2, c3: (1 - u) * c1 + u * ((1 - v) * c2 + v * c3)

    # transformation for the integrand,
    # including the Jacobian scaling factor
    def h(u, v):
        x = g(u, v, x1, x2, x3)
        y = g(u, v, y1, y2, y3)
        I = integrand(x, y)
        I *= dJ(u, v, p1, p2, p3)
        return I

    # perfrom the double integration using quadrature in the transformed space
    integral, error = integrate.dblquad(h, 0, 1, lambda x: 0, lambda x: 1, epsrel=1e-6, epsabs=0)
    return integral, error


def readData(filename):
    with open(filename) as f:
        [num] = [int(x) for x in next(f).split()]
        numbers = []
        for line in f:
            numbers.extend([float(x) for x in line.split()])
        # array = [[int(x) for x in line.split()] for line in f]
        # print numbers
        # exit(0)
        edges = []
        for i in range(0, num):
            edges.append((numbers.pop(0), numbers.pop(0)))
        F = []
        for i in range(0, num):
            F.append(numbers.pop(0))
        triangles_num = int(numbers.pop(0))
        triangles = []
        for i in range(0, triangles_num):
            a = int(numbers.pop(0))
            b = int(numbers.pop(0))
            c = int(numbers.pop(0))
            triangles.append((a, b, c))
            # print edges
    return (edges, F, triangles)


def printEigens(A):
    w, v = LA.eig(A)
    print w, ' = e.v. \n'


def getEigens(A):
    w, v = LA.eig(A)
    return w


# edges, F, triangles = readData('/home/artur/Programming/Python/images/c.in')
edges, F, triangles = readData('/home/artur/Programming/Python/images/workfile')

N = len(edges)
alpha = 1
betta = 1
A = np.zeros((N, N))
b = np.zeros(N)
np.set_printoptions(suppress=True)
np.set_printoptions(precision=3)
for triangle in triangles:
    (p1, p2, p3) = triangle
    (x1, y1), (x2, y2), (x3, y3) = edges[p1], edges[p2], edges[p3]
    f = triangleValue(F[p1], F[p2], F[p3])

    M = np.array([[x1, y1, 1], [x2, y2, 1], [x3, y3, 1]])
    rM = inv(M)
    r = lambda i, j: rM[i-1, j-1]
    # r = inv(M)

    Ix, _ = tridblquad(lambda x, y: x, (x1, y1), (x2, y2), (x3, y3))
    Iy, _ = tridblquad(lambda x, y: y, (x1, y1), (x2, y2), (x3, y3))
    Ixy, _ = tridblquad(lambda x, y: x*y, (x1, y1), (x2, y2), (x3, y3))
    Ix2, _ = tridblquad(lambda x, y: x**2, (x1, y1), (x2, y2), (x3, y3))
    Iy2, _ = tridblquad(lambda x, y: y**2, (x1, y1), (x2, y2), (x3, y3))
    # print 'Ixy', Ix, Iy
    delta, _ = tridblquad(lambda x, y: 1, (x1, y1), (x2, y2), (x3, y3))
    # delta = 1
    # f = 1
    r11, r12, r13, r21, r22, r23, r31, r32, r33 = r(1,1),r(1,2),r(1,3),r(2,1),r(2,2),r(2,3),r(3,1),r(3,2),r(3,3)
    A[p1, p1] += r11**2 * (alpha*delta + betta*Ix2) + r21**2 * (alpha*delta + betta*Iy2) + r31**2 * betta *delta \
                 + r11*r21*Ixy + r11*r31*Ix + r21*r31*Iy

    A[p2, p2] += r12**2 * (alpha*delta + betta*Ix2) + r22**2 * (alpha*delta + betta*Iy2) + r32**2 * betta *delta \
                 + r12*r22*Ixy + r12*r32*Ix + r22*r32*Iy

    A[p3, p3] += r13**2 * (alpha*delta + betta*Ix2) + r23**2 * (alpha*delta + betta*Iy2) + r33**2 * betta *delta \
                 + r13*r23*Ixy + r13*r33*Ix + r23*r33*Iy

    A[p1, p2] += r11*r12*(alpha*delta + betta*Ix2) + r21*r22*(alpha*delta + betta*Iy2) + r31*r32*betta*delta \
                 + Ixy*(r11*r22 + r12*r21) + Ix*(r11*r32 + r12*r31) + Iy*(r21*r32+r22*r31)

    A[p1, p3] += r11*r13*(alpha*delta + betta*Ix2) + r21*r23*(alpha*delta + betta*Iy2) + r31*r33*betta*delta \
                 + Ixy*(r11*r23 + r13*r21) + Ix*(r11*r33 + r13*r31) + Iy*(r21*r33+r23*r31)

    A[p2, p3] += r12*r13*(alpha*delta + betta*Ix2) + r22*r23*(alpha*delta + betta*Iy2) + r32*r33*betta*delta \
                 + Ixy*(r12*r23 + r13*r22) + Ix*(r12*r33 + r13*r32) + Iy*(r22*r33+r23*r32)

    A[p2, p1] = A[p1, p2]
    A[p3, p1] = A[p1, p3]
    A[p3, p2] = A[p2, p3]

    b[p1] += r31 * 2*f * betta*delta + r11*f*Ix + r21*f*Iy
    b[p2] += r32 * 2*f * betta*delta + r12*f*Ix + r22*f*Iy
    b[p3] += r33 * 2*f * betta*delta + r13*f*Ix + r23*f*Iy


b = np.array(b)
for i in range(0, len(edges)):
    A[i, i] *= 2

w = getEigens(A)


excluded = -1
for i in range(0, len(w)):
    if abs(w[i]) < 0.0000001:
        print 'Exclude ', i
        excluded = i
        A = np.delete(A, i, 0)
        A = np.delete(A, i, 1)
        b = np.delete(b, i, 0)

# u = np.linalg.solve(A, b)
# print u

L = linalg.cholesky(A)
Lt = L.transpose()

y = np.linalg.solve(L, b)
u = np.linalg.solve(Lt, y)
# print u

if excluded > -1:
    u = np.insert(u, excluded, 0)

# print u
d = int(sqrt(N))
print u.reshape(d, d)

# pp = pprint.PrettyPrinter(indent=4)


ans = []
func = lambda x, y : x * (abs(x)-2) / 2
# func = lambda x, y : 1/( np.pi**2/4*(2) ) * cos( pi*(x+1)/2 ) * cos( pi*(y+1)/2 )

const = 0
for i in range(0,len(edges)):
    x = edges[i][0]
    y = edges[i][1]
    ans.append(func(x, y))

ans = np.asarray(ans)
print
# print ans.reshape(d, d)
