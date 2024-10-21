import pandas as pd
import matplotlib.pyplot as plt

def plot_time(name):
  input_csv = pd.read_csv("./" + name + ".csv")
  Ns = input_csv[input_csv.keys()[0]]
  exact_time = input_csv[input_csv.keys()[1]]
  nrnmm_1_time = input_csv[input_csv.keys()[2]]
  nrnmm_2_time = input_csv[input_csv.keys()[3]]
  nrnmm_3_time = input_csv[input_csv.keys()[4]]

  plt.xlabel("N")
  plt.ylabel("time[ms]")

  plt.xscale("log")
  plt.yscale("log")

  plt.plot(Ns, exact_time, label="exact")
  plt.plot(Ns, nrnmm_1_time, label="nrnmm order 1")
  plt.plot(Ns, nrnmm_2_time, label="nrnmm order 2")
  plt.plot(Ns, nrnmm_3_time, label="nrnmm order 3")
  plt.legend()
  plt.savefig(name)

  plt.clf()
  plt.close()

plot_time("time_1")
plot_time("time_2")
