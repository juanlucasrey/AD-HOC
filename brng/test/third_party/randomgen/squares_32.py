from randomgen import Squares
rg = Squares(1234, variant=32)

state = rg.state["state"]
with open("squares_32_state.txt", "w") as f:
    f.write(str(state["key"]) + "\n")


iterations = 100
with open("squares_32_vals.txt", "w") as f:
    for _ in range(iterations):
        f.write(str(rg.random_raw()) + "\n")
