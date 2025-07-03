from randomgen import Xoshiro512
rg = Xoshiro512(1234)

state = rg.state
with open("xoshiro512_state.txt", "w") as f:
    p = state["s"]
    for x in range(len(p)):
        f.write(str(p[x]) + "\n")

iterations = 100
with open("xoshiro512_vals.txt", "w") as f:
    for _ in range(iterations):
        f.write(str(rg.random_raw()) + "\n")
