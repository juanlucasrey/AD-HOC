from randomgen import SPECK128
rg = SPECK128(1234)

state = rg.state["state"]
with open("speck_state.txt", "w") as f:
    p = state["round_key"]
    for x in range(len(p)):
        f.write(str(p[x]) + "\n")
    p = state["ctr"]
    for x in range(len(p)):
        f.write(str(p[x]) + "\n")


iterations = 100
with open("speck_vals.txt", "w") as f:
    for _ in range(iterations):
        f.write(str(rg.random_raw()) + "\n")
