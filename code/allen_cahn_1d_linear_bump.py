from fenics import *

dt = 5.0e-5
n = 100
eps = 1/n
t = 0.00
T = .1

class LinearBump(UserExpression):
    def eval(self, values, x):
        if x[0] < 0.25 + DOLFIN_EPS:
            values[0] = -1
        elif x[0] < .25 + eps + DOLFIN_EPS:
            values[0] = 2/eps *(x[0]-.25) - 1
        elif x[0] < .75 - eps + DOLFIN_EPS:
            values[0] = 1
        elif x[0] < .75 + DOLFIN_EPS:
            values[0] = 2/eps *(.75 - eps - x[0]) + 1
        else:
            values[0] = -1

class PeriodicBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return bool(x[0] < DOLFIN_EPS and x[0] > -DOLFIN_EPS and on_boundary)

    def map(self, x, y):
        y[0] = x[0] - 1

mesh = UnitIntervalMesh(300)
V = FunctionSpace(mesh, 'CG', 1, constrained_domain=PeriodicBoundary())

u, v = Function(V), TestFunction(V)
u.rename('u', '')

u_init = LinearBump()
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

file = XDMFFile("data_allen_cahn_1d_linear_bump.xdmf")
file.parameters["functions_share_mesh"] = True
i = 0
file.write(u_pre, i)

while t < T:
    print(f'Time {t} of {T}.')
    
    i += 1
    t += dt
    solve(F == 0, u)

    file.write(u, i)
    u_pre.assign(u)
