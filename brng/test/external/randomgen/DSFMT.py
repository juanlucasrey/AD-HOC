from randomgen import DSFMT
from numpy.random import Generator

rg = DSFMT(5489)

state = rg.state["state"]
with open("DSFMT_state.txt", "w") as f:
    p = state["state"]
    for x in range(len(p)):
        f.write(str(p[x]) + "\n")
    f.write(str(state["idx"]) + "\n")

iterations = 100
with open("DSFMT_vals.txt", "w") as f:
    for _ in range(iterations):
        f.write(str(rg.random_raw()) + "\n")

rng2 = Generator(DSFMT(5489))
with open("DSFMT_double_vals.txt", "w") as f:
    for _ in range(iterations):
        f.write(str(rng2.uniform()) + "\n")
