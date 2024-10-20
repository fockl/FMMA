import pandas as pd
import matplotlib.pyplot as plt

input_csv = pd.read_csv('./error.csv')
Ns = input_csv[input_csv.keys()[0]]
nrnmm_1_error = input_csv[input_csv.keys()[2]]
nrnmm_2_error = input_csv[input_csv.keys()[3]]
nrnmm_3_error = input_csv[input_csv.keys()[4]]

plt.xlabel("N")
plt.ylabel("error")

plt.xscale("log")
plt.yscale("log")

plt.plot(Ns, nrnmm_1_error, label="nrnmm order 1")
plt.plot(Ns, nrnmm_2_error, label="nrnmm order 2")
plt.plot(Ns, nrnmm_3_error, label="nrnmm order 3")
plt.legend()
plt.savefig("error")

