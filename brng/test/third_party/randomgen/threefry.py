from randomgen import ThreeFry
rg = ThreeFry(1234)

state = rg.state["state"]
with open("threefry_state.txt", "w") as f:
    p = state["key"]
    for x in range(len(p)):
        f.write(str(p[x]) + "\n")
    p = state["counter"]
    for x in range(len(p)):
        f.write(str(p[x]) + "\n")


iterations = 100
with open("threefry_vals.txt", "w") as f:
    for _ in range(iterations):
        f.write(str(rg.random_raw()) + "\n")
