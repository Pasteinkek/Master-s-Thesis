from fenics import *

dt = 5.0e-3
eps = 0.01
t = 0.1
T = 1

class PeriodicBoundary(SubDomain):

    def inside(self, x, on_boundary):
        return on_boundary

    def map(self, x, y):
        for i in range(2):
            if near(x[i], 0):
                y[i] = 1
            elif near(x[i], 1):
                y[i] = 0
            else:
                raise Exception()

mesh = UnitSquareMesh(100, 100)
pbc = PeriodicBoundary()
V = FunctionSpace(mesh, 'CG', 1, constrained_domain=pbc)

u, v = Function(V), TestFunction(V)
u.rename('u', '')

u_init = Expression("pow(x[0],2)*sin(2*pi*x[0])", degree=2)
u_pre = Function(V)
u_pre.interpolate(u_init)

def W(u):
    return (u**2 - 1)**2

F = - dt * inner(grad(u), grad(v)) * dx \
    - u * v * dx \
    + u_pre * v * dx \
    - dt / eps**2 * diff(W(u), u) * v * dx

file = XDMFFile("data.xdmf")
file.parameters["functions_share_mesh"] = True
i = 0
file.write(u, i)

while t < T:
    print(f'Time {t} of {T}.')
    
    i += 1
    t += dt
    solve(F == 0, u)

    file.write(u, i)
    u_pre.assign(u)
