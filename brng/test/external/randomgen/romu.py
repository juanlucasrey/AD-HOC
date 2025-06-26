from randomgen import Romu
rg = Romu(123)

state = rg.state["state"]
with open("romu_state.txt", "w") as f:
    f.write(str(state['w']) + "\n")
    f.write(str(state['x']) + "\n")
    f.write(str(state['y']) + "\n")
    f.write(str(state['z']) + "\n")

iterations = 100
with open("romu_vals.txt", "w") as f:
    for _ in range(iterations):
        f.write(str(rg.random_raw()) + "\n")
