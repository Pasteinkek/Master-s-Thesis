from dolfin import *

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--eps', required=True, type=float, default=0.01)
# eps=.2 is good and eps=0.01

args = parser.parse_args()
eps = args.eps

dt = 5.0e-5
t = 0.00
T = 0.05

class LinearBump2D(UserExpression):
    def eval(self, values, x):
        d = sqrt((x[0] - .5)**2 + (x[1] - .5)**2)

        if d + DOLFIN_EPS < .25:
            values[0] = 1
        elif d + DOLFIN_EPS < .25 + eps:
            values[0] = - 2/eps * (d - .25 - eps) - 1
        else:
            values[0] = -1

class PeriodicBoundary(SubDomain):
    # based on https://fenicsproject.org/qa/262/possible-specify-more-than-one-periodic-boundary-condition/
    # Left boundary is "target domain" G
    def inside(self, x, on_boundary):
        return bool((near(x[0], 0) or near(x[1], 0)) and on_boundary)

    def map(self, x, y):
        y[0] = x[0]
        if near(x[0], 1):
            y[0] -= 1
        y[1] = x[1]
        if near(x[1], 1):
            y[1] -= 1


mesh = UnitSquareMesh(100, 100)
pbc = PeriodicBoundary()
V = FunctionSpace(mesh, 'CG', 1, constrained_domain=pbc)

u, v = Function(V), TestFunction(V)
u.rename('u', '')

u_init = LinearBump2D()
u_init = interpolate(u_init, V)

u_pre = Function(V)
u_pre.rename('u', '')
u_pre.interpolate(u_init)

def W(u):
    return (u**2 - 1)**2

def W_prime(u):
    return 4 * u * (u**2 - 1)

F = + 1/dt * u * v * dx \
    - 1/dt * u_pre * v * dx \
    + inner(grad(u), grad(v)) * dx \
    + 1/eps**2 * W_prime(u) * v * dx

file = XDMFFile(f"data_allen_cahn_2d_linear_bump_eps_{eps}.xdmf")
file.parameters["functions_share_mesh"] = True
i = 0
file.write(u_pre, i)

while t < T:
    print(f'Time {t} of {T}.')
    
    i += 1
    t += dt
    solve(F == 0, u)

    if i % 10 == 0:
        file.write(u, i)

    u_pre.assign(u)

file.write(