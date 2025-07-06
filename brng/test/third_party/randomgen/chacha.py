from randomgen import ChaCha
rg = ChaCha(1234)

state = rg.state["state"]
with open("chacha_state.txt", "w") as f:
    p = state["keysetup"]
    for x in range(len(p)):
        f.write(str(p[x]) + "\n")
    p = state["ctr"]
    for x in range(len(p)):
        f.write(str(p[x]) + "\n")
    p = state["block"]
    for x in range(len(p)):
        f.write(str(p[x]) + "\n")


iterations = 100
with open("chacha_vals.txt", "w") as f:
    for _ in range(iterations):
        f.write(str(rg.random_raw()) + "\n")
