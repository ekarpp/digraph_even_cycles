#!/usr/bin/env python3

def C(n, path):
    l = [str((x+1)%n) for x in range(n)]
    with open(f"{path}/c{n}", "w") as f:
        f.write(
            "\n".join(l) + "\n"
        )

import sys

n = int(sys.argv[1])
p = sys.argv[2]

C(n, p)
