
import sys

out = open(sys.argv[1], "w")
s = raw_input()
out.write(s + "\n")
Np = int(s)
for i in range(Np):
    s = raw_input()
    out.write(s + "\n")
s = raw_input()
out.write(s + "\n")
