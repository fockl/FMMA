import pandas as pd
import matplotlib.pyplot as plt

def plot_error(name):
  input_csv = pd.read_csv("./" + name + ".csv")
  Ns = input_csv[input_csv.keys()[0]]
  nrnmm_1_error = input_csv[input_csv.keys()[1]]
  nrnmm_2_error = input_csv[input_csv.keys()[2]]
  nrnmm_3_error = input_csv[input_csv.keys()[3]]
  tree_1_error = input_csv[input_csv.keys()[4]]
  tree_2_error = input_csv[input_csv.keys()[5]]
  tree_3_error = input_csv[input_csv.keys()[6]]
  fmm_1_error = input_csv[input_csv.keys()[7]]
  fmm_2_error = input_csv[input_csv.keys()[8]]
  fmm_3_error = input_csv[input_csv.keys()[9]]

  plt.xlabel("N")
  plt.ylabel("relative error")

  plt.xscale("log")
  plt.yscale("log")

  plt.plot(Ns, nrnmm_1_error, label="nrnmm order 1")
  plt.plot(Ns, nrnmm_2_error, label="nrnmm order 2")
  plt.plot(Ns, nrnmm_3_error, label="nrnmm order 3")
  plt.plot(Ns, tree_1_error, label="tree order 1")
  plt.plot(Ns, tree_2_error, label="tree order 2")
  plt.plot(Ns, tree_3_error, label="tree order 3")
  plt.plot(Ns, fmm_1_error, label="fmm order 1")
  plt.plot(Ns, fmm_2_error, label="fmm order 2")
  plt.plot(Ns, fmm_3_error, label="fmm order 3")
  plt.legend()
  plt.savefig(name)

  plt.clf()
  plt.close()

plot_error("error_1")
plot_error("error_2")
