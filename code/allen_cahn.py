from fenics import *

dt = 0.05
n = 2
eps = 1/n
t = 0.00
T = 0.25

class PeriodicBoundary(SubDomain):
    # based on https://fenicsproject.org/qa/262/possible-specify-more-than-one-periodic-boundary-condition/
    # Left boundary is "target domain" G
    def inside(self, x, on_boundary):
        # return True if on left or bottom boundary AND NOT on one of the two corners (0, 1) and (1, 0)
        return bool((near(x[0], 0) or near(x[1], 0)) and 
                (not ((near(x[0], 0) and near(x[1], 1)) or 
                        (near(x[0], 1) and near(x[1], 0)))) and on_boundary)

    def map(self, x, y):
        y[0] = x[0] + near(x[0], 1) * (-1) 
        y[1] = x[1] + near(x[1], 1) * (-1) 

mesh = UnitSquareMesh(100, 100)
# mesh = UnitIntervalMesh(100)
pbc = PeriodicBoundary()
V = FunctionSpace(mesh, 'CG', 1, constrained_domain=pbc)

u, v = Function(V), TestFunction(V)
u.rename('u', '')

# u_init = Expression("pow(x[0],2)*sin(2*pi*x[0])", degree=2)
# u_init = Expression("sin(pi/eps * (x[0] + x[1]))", degree=3, eps=eps)
u_init = Expression("sin(2*pi/eps * (x[0]))", degree=3, eps=eps)
u_pre = Function(V)
u_pre.rename('u', '')
u_pre.interpolate(u_init)

def W(u):
    return (u**2 - 1)**2

# dt u = \laplace u - 1/eps**2 diff(W(u), u)
# dt u - \laplace u + 1/eps**2 diff(W(u), u) = 0
# (u - u_prev)/dt
F = + 1/dt * u * v * dx \
    - 1/dt * u_pre * v * dx \
    + inner(grad(u), grad(v)) * dx \
    + 1/eps**2 * diff(W(u), u) * v * dx

file = XDMFFile("data.xdmf")
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
