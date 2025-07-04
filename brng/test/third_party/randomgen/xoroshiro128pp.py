from randomgen import Xoroshiro128
rg = Xoroshiro128(1234, plusplus=True)

state = rg.state
with open("xoroshiro128pp_state.txt", "w") as f:
    p = state["s"]
    for x in range(len(p)):
        f.write(str(p[x]) + "\n")

iterations = 100
with open("xoroshiro128pp_vals.txt", "w") as f:
    for _ in range(iterations):
        f.write(str(rg.random_raw()) + "\n")
