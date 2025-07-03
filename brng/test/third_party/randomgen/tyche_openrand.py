from randomgen import Tyche
rg = Tyche(123, original=False)

state = rg.state["state"]
with open("tyche_openrand_state.txt", "w") as f:
    f.write(str(state['a']) + "\n")
    f.write(str(state['b']) + "\n")
    f.write(str(state['c']) + "\n")
    f.write(str(state['d']) + "\n")

iterations = 100
with open("tyche_openrand_vals.txt", "w") as f:
    for _ in range(iterations):
        f.write(str(rg.random_raw()) + "\n")
