from randomgen import SFMT
rg = SFMT(5489)

state = rg.state["state"]
with open("SFMT_state.txt", "w") as f:
    p = state["state"]
    for x in range(len(p)):
        f.write(str(p[x]) + "\n")
    f.write(str(state["idx"]) + "\n")

iterations = 100
with open("SFMT_vals.txt", "w") as f:
    for _ in range(iterations):
        f.write(str(rg.random_raw()) + "\n")
