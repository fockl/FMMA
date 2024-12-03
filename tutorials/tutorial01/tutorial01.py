import pyfmma
import numpy as np

def fn(x, y):
  d = np.fabs(x[0]-y[0])
  return 1.0/(d*d)

fmma = pyfmma.fmmad1()
fmma.set_fn2(fn)
fmma.set_solver_type("exact")

N = 10
source = np.zeros((N, 1))
target = np.zeros((N, 1))
ans = np.zeros(N)

for i in range(N):
  source[i][0] = (i+0.0)/N
  target[i][0] = (i+0.5)/N

fmma.solve(target, source, ans)

print(ans)
