from randomgen import MT64
rg = MT64(5489)

state = rg.state["state"]
with open("mt64_state.txt", "w") as f:
    p = state["key"]
    for x in range(len(p)):
        f.write(str(p[x]) + "\n")
    f.write(str(state["pos"]) + "\n")


iterations = 100
with open("mt64_vals.txt", "w") as f:
    for _ in range(iterations):
        f.write(str(rg.random_raw()) + "\n")
