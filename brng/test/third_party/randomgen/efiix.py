from randomgen import EFIIX64
rg = EFIIX64(1234)

state = rg.state["state"]
with open("efiix_state.txt", "w") as f:
    f.write(str(state["a"]) + "\n")
    f.write(str(state["b"]) + "\n")
    f.write(str(state["c"]) + "\n")
    f.write(str(state["i"]) + "\n")
    iteration_table = state["iteration_table"]
    for x in range(len(iteration_table)):
        f.write(str(iteration_table[x]) + "\n")
    indirection_table = state["indirection_table"]
    for x in range(len(indirection_table)):
        f.write(str(indirection_table[x]) + "\n")

iterations = 100
with open("efiix_vals.txt", "w") as f:
    for _ in range(iterations):
        f.write(str(rg.random_raw()) + "\n")
