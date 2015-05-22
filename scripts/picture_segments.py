
from matplotlib import pyplot as plt
import sys

filename = sys.argv[1]

# could do in pairs but not convinient
def ReadPoints(f):
    X, Y = [], []
    n = int(f.readline())
    for i in range(n):
        x, y = map(int, f.readline().split())
        X.append(x)
        Y.append(y)
    return X, Y

# pairs of indices
def ReadSegments(f):
    S = []
    n = int(f.readline())
    for i in range(n):
        S.append(tuple(map(int, f.readline().split())))
    return S
    
# returns intersection points
def ReadIntersections(f):
    X, Y = [], []
    n = int(f.readline())
    for i in range(n):
        # first two elements are segments
        x, y = map(float, f.readline().split()[2:])
        X.append(x)
        Y.append(y)
    return X, Y


f = open(filename)
X, Y = ReadPoints(f)
S = ReadSegments(f)
XI, YI = ReadIntersections(f)
f.close()

XS, YS = [], []
for i, s in enumerate(S):
    x_0, x_1 = X[s[0]], X[s[1]]
    y_0, y_1 = Y[s[0]], Y[s[1]]
    XS.extend((x_0, x_1, None))
    YS.extend((y_0, y_1, None))
    plt.text((x_0+x_1)/2, (y_0 + y_1)/2, str(i))
    
    
plt.plot(XS, YS, "ro-")
plt.plot(XI, YI, "go")
plt.show()