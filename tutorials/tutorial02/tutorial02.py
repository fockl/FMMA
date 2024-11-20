import pyfmma
import numpy as np

def fn(x_y):
  d = np.linalg.norm(x_y)
  return 1.0/(d*d)

fmma = pyfmma.fmmad1()
fmma.set_fn(fn)
fmma.set_solver_type("tree")
fmma.set_Depth(3)
fmma.set_poly_ord(3)

N = 10
source = np.zeros((N, 1))
target = np.zeros((N, 1))
ans = np.zeros(N)

for i in range(N):
  source[i][0] = (i+0.0)/N
  target[i][0] = (i+0.5)/N

fmma.solve(target, source, ans)

print(ans)
