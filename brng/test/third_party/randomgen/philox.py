from randomgen import Philox
rg = Philox(1234)

state = rg.state["state"]
with open("philox_state.txt", "w") as f:
    p = state["key"]
    for x in range(len(p)):
        f.write(str(p[x]) + "\n")
    p = state["counter"]
    for x in range(len(p)):
        f.write(str(p[x]) + "\n")


iterations = 100
with open("philox_vals.txt", "w") as f:
    for _ in range(iterations):
        f.write(str(rg.random_raw()) + "\n")
