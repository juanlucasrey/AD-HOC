from randomgen import JSF
rg = JSF(1234)

state = rg.state["state"]
with open("jsf_state.txt", "w") as f:
    f.write(str(state["a"]) + "\n")
    f.write(str(state["b"]) + "\n")
    f.write(str(state["c"]) + "\n")
    f.write(str(state["d"]) + "\n")
    f.write(str(state["p"]) + "\n")
    f.write(str(state["q"]) + "\n")
    f.write(str(state["r"]) + "\n")


iterations = 100
with open("jsf_vals.txt", "w") as f:
    for _ in range(iterations):
        f.write(str(rg.random_raw()) + "\n")
