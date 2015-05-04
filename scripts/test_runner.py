
from subprocess import Popen, PIPE
import sys

input_name = sys.argv[1]
output_name = sys.argv[2]

# input
out = open(input_name, "w")
s = raw_input()
out.write(s + "\n")
Np = int(s)
for i in range(Np):
    s = raw_input()
    out.write(s + "\n")
s = raw_input()
out.write(s + "\n")
out.close()

Popen(["./../bin/small_polygons", input_name, output_name], 
      stdout=PIPE, stderr=PIPE, stdin=PIPE).wait()

input = open(output_name, "r")
tour_count = int(input.readline())
print tour_count
for t in range(tour_count):
    print input.readline()[:-1]
    