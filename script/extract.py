import os, numpy as np, opt_einsum


for folder in list(os.walk("MBPT"))[0][1]:
    contractions = []
    for line in open("MBPT/" + folder + "/CD_output.txt", "r").readlines()[:-1]:
        contractions.append([",".join([contr.strip("{};").split("{{")[-1].replace(", ", "") for contr in line.strip().split("}, ")]), str(1.0 / eval(line.split("{")[1].strip(", ")))])
    for contr in contractions:
        tens = []
        for elem in contr[0].split(","):
            tens.append(np.random.rand(*[7 for i in elem]))
        path = opt_einsum.contract_path(contr[0], *tens)
        contr[0] = contr[0] + ";" + ":".join([("{}-{}/".format(*e[0]) + e[2].replace(",", "+")) for e in path[1].__dict__["contraction_list"]])
    for contr in contractions: contr = ";".join(contr)

    contractions = [";".join(contr).replace(",", "+") for contr in contractions]

    with open("mbpt{}.txt".format(folder.split("-")[1]), "w") as output:
        output.write("TOSTRING(\n"); output.write("\n".join(contractions)); output.write("\n)\n")

