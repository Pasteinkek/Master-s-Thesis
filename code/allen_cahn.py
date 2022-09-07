from dolfin import *

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--eps', required=True, type=float)

args = parser.parse_args()
eps = args.eps

dt = 5.0e-5
# dt = 5.0e-4
# n = 100
# eps = 1/n
t = 0.00
# T = 0.05
T = 100 * dt
# T = 1*dt


# class LinearBump(UserExpression):
#     def eval(self, values, x):
#         if x[0] < 0.25 + DOLFIN_EPS:
#             values[0] = -1
#         elif x[0] < .25 + eps + DOLFIN_EPS:
#             values[0] = 2/eps *(x[0]-.25) - 1
#         elif x[0] < .75 - eps + DOLFIN_EPS:
#             values[0] = 1
#         elif x[0] < .75 + DOLFIN_EPS:
#             values[0] = 2/eps *(.75 - eps - x[0]) + 1
#         else:
#             values[0] = -1

# def distance(a, b):
#     return sqrt((a[0]-b[0])**2 + (a[1]-b[1])**2)

# class Dumbel2D(UserExpression):
#     def eval(self, values, x):
#         r = 3/16
#         a = 3/8

#         if x[0] < .5 + DOLFIN_EPS:
#             d = distance(x, (a, .5))
#             if d < r + DOLFIN_EPS:
#                 values[0] = 1
#             elif d < r + eps + DOLFIN_EPS:
#                 values[0] = 2/eps * (r - d) + 1
#             else:
#                 values[0] = -1
#         else:
#             d = distance(x, (1-a, .5))
#             if d < r + DOLFIN_EPS:
#                 values[0] = 1
#             elif d < r + eps + DOLFIN_EPS:
#                 values[0] = 2/eps * (r - d) + 1
#             else:
#                 values[0] = -1

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

(x0, xf) = (-1.0, 1.0) # boundaries

# class PeriodicBoundary(SubDomain):

#     # Left boundary is "target domain" G
#     def inside(self, x, on_boundary):
#         return bool(x[0] < DOLFIN_EPS and x[0] > -DOLFIN_EPS and on_boundary)

#     # Map right boundary (H) to left boundary (G)
#     def map(self, x, y):
#         y[0] = x[0] - (xf - x0)

mesh = UnitSquareMesh(100, 100)
# mesh = UnitIntervalMesh(300)
# mesh = IntervalMesh(100, x0, xf)
# mesh = RectangleMesh(Point(-1, -1), Point(1, 1), 100, 100)

pbc = PeriodicBoundary()
V = FunctionSpace(mesh, 'CG', 1, constrained_domain=pbc)

u, v = Function(V), TestFunction(V)
u.rename('u', '')

# u_init = Expression("pow(x[0],2)*sin(2*pi*x[0])", degree=2)
u_init = Expression("sin(pi/eps * (x[0] + x[1]))", degree=3, eps=eps)
# u_init = Expression("sin(2*pi/eps * (x[0]))", degree=3, eps=eps)
# u_init = Expression("pow(x[0],2)*sin(2*pi*x[0])", degree=2)
# u_init = Dumbel2D()
u_init = interpolate(u_init, V)
# def boundary(x, on_boundary):
#     return on_boundary
# dbc = [DirichletBC(V, u_init, boundary)]
u_pre = Function(V)
u_pre.rename('u', '')
u_pre.interpolate(u_init)

def W(u):
    return (u**2 - 1)**2

def W_prime(u):
    return 4 * u * (u**2 - 1)

# dt u = \laplace u - 1/eps**2 diff(W(u), u)
# dt u - \laplace u + 1/eps**2 diff(W(u), u) = 0
# (u - u_prev)/dt
F = + 1/dt * u * v * dx \
    - 1/dt * u_pre * v * dx \
    + inner(grad(u), grad(v)) * dx \
    + 1/eps**2 * W_prime(u) * v * dx
    # + 1/eps**2 * diff(W(u), u) * v * dx

# J = derivative(F, u, TrialFunction(V))
# print(diff(W(u), u))
# exit()
file = XDMFFile("data.xdmf")
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

file.write(u, i)