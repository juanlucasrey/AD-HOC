from randomgen import HC128
rg = HC128(1234)
result = rg.state

state = rg.state["state"]
with open("hc128state.txt", "w") as f:
    p = state["p"]
    for x in range(len(p)):
        f.write(str(p[x]) + "\n")
    q = state["q"]
    for x in range(len(q)):
        f.write(str(q[x]) + "\n")

iterations = 100
with open("hc128vals.txt", "w") as f:
    for _ in range(iterations):
        f.write(str(rg.random_raw()) + "\n")
