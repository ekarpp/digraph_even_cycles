#!/usr/bin/env python3

def K(n, path):
	l = [x for x in range(n)]
	grph = ""
	for i in range(n):
		lfilt = [str(x) for x in l if x != i]
		grph += f"{' '.join(lfilt)}\n"
	with open(f"{path}/k{n}", "w") as f:
		f.write(grph)

import sys

n = int(sys.argv[1])
p = sys.argv[2]

K(n, p)
