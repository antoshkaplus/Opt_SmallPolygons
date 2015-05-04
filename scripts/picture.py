
from matplotlib import pyplot as plt
import sys

points_filename = sys.argv[1]
path_filename = sys.argv[2]

points_file = open(points_filename)
N = int(points_file.readline())
# nothing after last line break
lines = points_file.read().split("\n")
coordinates = map(int, lines[:-1])
X = coordinates[::2]
Y = coordinates[1::2]

points_file.close()

path_file = open(path_filename)
path_file.readline()
paths  = map(lambda x: map(int, x.split(" ")), path_file.read().split("\n")[:-1])

for i, (x, y) in enumerate(zip(X, Y)):
    plt.text(x, y, i)

for p in paths:
    OP = [(X[c], Y[c]) for c in p]
    OP.append((X[p[0]], Y[p[0]]))
    x, y = zip(*OP)
    plt.plot(x, y, "ro-")

plt.show()