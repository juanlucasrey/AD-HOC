from randomgen import Squares
rg = Squares(1234)

state = rg.state["state"]
with open("squares_state.txt", "w") as f:
    f.write(str(state["key"]) + "\n")


iterations = 100
with open("squares_vals.txt", "w") as f:
    for _ in range(iterations):
        f.write(str(rg.random_raw()) + "\n")
