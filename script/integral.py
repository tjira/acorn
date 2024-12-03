import pyscf

molecule = pyscf.M(atom="".join(open("example/molecule/water.xyz").readlines()[2:]), basis="sto-3g")

S, T, V, J = molecule.intor("int1e_ovlp"), molecule.intor("int1e_kin"), molecule.intor("int1e_nuc"), molecule.intor("int2e")

D1, D2, J = " ".join([str(D) for D in S.shape]), " ".join([str(D) for D in J.shape]), J.reshape((J.shape[0] * J.shape[1], J.shape[2] * J.shape[3]))

with open("S_AO.mat", "w") as file:
    file.write(f"{D1}\n" + "\n".join([" ".join(["%20.14f" % S[i, j] for j in range(S.shape[1])]) for i in range(S.shape[0])]) + "\n")

with open("T_AO.mat", "w") as file:
    file.write(f"{D1}\n" + "\n".join([" ".join(["%20.14f" % T[i, j] for j in range(T.shape[1])]) for i in range(T.shape[0])]) + "\n")

with open("V_AO.mat", "w") as file:
    file.write(f"{D1}\n" + "\n".join([" ".join(["%20.14f" % V[i, j] for j in range(V.shape[1])]) for i in range(V.shape[0])]) + "\n")

with open("J_AO.mat", "w") as file:
    file.write(f"{D2}\n" + "\n".join([" ".join(["%20.14f" % J[i, j] for j in range(J.shape[1])]) for i in range(J.shape[0])]) + "\n")
