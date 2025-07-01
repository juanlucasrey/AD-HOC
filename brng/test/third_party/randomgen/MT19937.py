from randomgen import MT19937
rg = MT19937(5489)

state = rg.state["state"]
with open("mt19937_state.txt", "w") as f:
    p = state["key"]
    for x in range(len(p)):
        f.write(str(p[x]) + "\n")
    f.write(str(state["pos"]) + "\n")


iterations = 100
with open("mt19937_vals.txt", "w") as f:
    for _ in range(iterations):
        f.write(str(rg.random_raw()) + "\n")
