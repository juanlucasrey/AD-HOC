from randomgen import SFC64
rg = SFC64(5489)

state = rg.state["state"]
with open("sfc64_state.txt", "w") as f:
    f.write(str(state["a"]) + "\n")
    f.write(str(state["b"]) + "\n")
    f.write(str(state["c"]) + "\n")
    f.write(str(state["w"]) + "\n")
    f.write(str(state["k"]) + "\n")


iterations = 100
with open("sfc64_vals.txt", "w") as f:
    for _ in range(iterations):
        f.write(str(rg.random_raw()) + "\n")
