import pandas as pd
import matplotlib.pyplot as plt

def plot_time(name):
  input_csv = pd.read_csv("./" + name + ".csv")
  Ns = input_csv[input_csv.keys()[0]]
  exact_time = input_csv[input_csv.keys()[1]]
  nrnmm_1_time = input_csv[input_csv.keys()[2]]
  nrnmm_1_error = input_csv[input_csv.keys()[3]]
  nrnmm_2_time = input_csv[input_csv.keys()[4]]
  nrnmm_2_error = input_csv[input_csv.keys()[5]]
  nrnmm_3_time = input_csv[input_csv.keys()[6]]
  nrnmm_3_error = input_csv[input_csv.keys()[7]]
  tree_1_time = input_csv[input_csv.keys()[8]]
  tree_1_error = input_csv[input_csv.keys()[9]]
  tree_2_time = input_csv[input_csv.keys()[10]]
  tree_2_error = input_csv[input_csv.keys()[11]]
  tree_3_time = input_csv[input_csv.keys()[12]]
  tree_3_error = input_csv[input_csv.keys()[13]]

  plt.xlabel("time[ms]")
  plt.ylabel("error")

  plt.xscale("log")
  plt.yscale("log")

  plt.plot(nrnmm_1_time, nrnmm_1_error, label="nrnmm order 1")
  plt.plot(nrnmm_2_time, nrnmm_2_error, label="nrnmm order 2")
  plt.plot(nrnmm_3_time, nrnmm_3_error, label="nrnmm order 3")
  plt.plot(tree_1_time, tree_1_error, label="tree order 1")
  plt.plot(tree_2_time, tree_2_error, label="tree order 2")
  plt.plot(tree_3_time, tree_3_error, label="tree order 3")
  plt.legend()
  plt.savefig(name)

  plt.clf()
  plt.close()

plot_time("all_1")
plot_time("all_2")
