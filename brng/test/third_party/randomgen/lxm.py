from randomgen import LXM
rg = LXM()

state = rg.state["state"]
with open("lxm_state.txt", "w") as f:
    f.write(str(state["lcg_state"]) + "\n")
    f.write(str(state["b"]) + "\n")
    p = state["x"]
    for x in range(len(p)):
        f.write(str(p[x]) + "\n")


iterations = 100
with open("lxm_vals.txt", "w") as f:
    for _ in range(iterations):
        f.write(str(rg.random_raw()) + "\n")
