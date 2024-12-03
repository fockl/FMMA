import pyfmma
import numpy as np

def fn(x, y):
  d = np.fabs((x[0]-y[0])*(x[0]-y[0]) + (x[1]-y[1])*(x[1]-y[1]))
  return 1.0/d

fmma = pyfmma.fmmad2()
fmma.set_fn2(fn)
fmma.set_solver_type("fmm")
fmma.set_Depth(3)
fmma.set_poly_ord(3)

N = 10
source = np.zeros((N, 2))
source_weight = np.zeros(N)
target = np.zeros((N, 2))
ans = np.zeros(N)

for i in range(N):
  source[i][0] = (i+0.0)/N
  source[i][1] = (i+0.0)/N
  target[i][0] = (i+0.5)/N
  target[i][1] = (i+0.5)/N
  source_weight[i] = (i+0.0)/N

fmma.solve(target, source_weight, source, ans)

print(ans)