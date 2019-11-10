from numpy import *             # Grab all of the NumPy functions
from matplotlib.pyplot import * # Grab MATLAB plotting functions
from control.matlab import *    # MATLAB-like functions

m = 4;                         # mass of aircraft
J = 0.0475;                    # inertia around pitch axis
r = 0.25;                      # distance to center of force
g = 9.8;                       # gravitational constant
c = 0.05;                      # damping factor (estimated)

# State space dynamics
xe = [0, 0, 0, 0, 0, 0];        # equilibrium point of interest
ue = [0, m*g];                  # (note these are lists, not matrices)

# Dynamics matrix (use matrix type so that * works for multiplication)
A = matrix(
   [[ 0,    0,    0,    1,    0,    0],
    [ 0,    0,    0,    0,    1,    0],
    [ 0,    0,    0,    0,    0,    1],
    [ 0, 0, (-ue[0]*sin(xe[2]) - ue[1]*cos(xe[2]))/m, -c/m, 0, 0],
    [ 0, 0, (ue[0]*cos(xe[2]) - ue[1]*sin(xe[2]))/m, 0, -c/m, 0],
    [ 0,    0,    0,    0,    0,    0 ]])

# Input matrix
B = matrix(
   [[0, 0], [0, 0], [0, 0],
    [cos(xe[2])/m, -sin(xe[2])/m],
    [sin(xe[2])/m,  cos(xe[2])/m],
    [r/J, 0]])

# Output matrix
C = matrix([[1, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0]])
D = matrix([[0, 0], [0, 0]])

Qx1 = diag([1, 1, 1, 1, 1, 1]);
Qu1a = diag([1, 1]);
(K, X, E) = lqr(A, B, Qx1, Qu1a); K1a = matrix(K);

xd = matrix([[1], [0], [0], [0], [0], [0]]);
yd = matrix([[0], [1], [0], [0], [0], [0]]);

# Indices for the parts of the state that we want
lat = (0,2,3,5);
alt = (1,4);

# Decoupled dynamics
Ax = (A[lat, :])[:, lat];       #! not sure why I have to do it this way
Bx = B[lat, 0]; Cx = C[0, lat]; Dx = D[0, 0];

Ay = (A[alt, :])[:, alt];       #! not sure why I have to do it this way
By = B[alt, 1]; Cy = C[1, alt]; Dy = D[1, 1];

# Step response for the first input
H1ax = ss(Ax - Bx*K1a[0,lat], Bx*K1a[0,lat]*xd[lat,:], Cx, Dx);
(Tx, Yx) = step(H1ax, T=linspace(0,10,100));

# Step response for the second input
H1ay = ss(Ay - By*K1a[1,alt], By*K1a[1,alt]*yd[alt,:], Cy, Dy);
(Ty, Yy) = step(H1ay, T=linspace(0,10,100));

plot(Tx, Yx[0,:].T, '-', Ty, Yy[0,:].T, '--'); hold(True);
plot([0, 10], [1, 1], 'k-'); hold(True);
ylabel('position');
legend(('x', 'y'), loc='lower right');

# Look at different input weightings
Qu1a = diag([1, 1]); (K1a, X, E) = lqr(A, B, Qx1, Qu1a);
H1ax = ss(Ax - Bx*K1a[0,lat], Bx*K1a[0,lat]*xd[lat,:], Cx, Dx);

Qu1b = (40**2)*diag([1, 1]); (K1b, X, E) = lqr(A, B, Qx1, Qu1b);
H1bx = ss(Ax - Bx*K1b[0,lat], Bx*K1b[0,lat]*xd[lat,:],Cx, Dx);

Qu1c = (200**2)*diag([1, 1]); (K1c, X, E) = lqr(A, B, Qx1, Qu1c);
H1cx = ss(Ax - Bx*K1c[0,lat], Bx*K1c[0,lat]*xd[lat,:],Cx, Dx);

[T1, Y1] = step(H1ax, T=linspace(0,10,100));
[T2, Y2] = step(H1bx, T=linspace(0,10,100));
[T3, Y3] = step(H1cx, T=linspace(0,10,100));

plot(T1, Y1[0,:].T, 'b-'); hold(True);
plot(T2, Y2[0,:].T, 'b-'); hold(True);
plot(T3, Y3[0,:].T, 'b-'); hold(True);
plot([0 ,10], [1, 1], 'k-'); hold(True);

axis([0, 10, -0.1, 1.4]);
# arcarrow([1.3, 0.8], [5, 0.45], -6);
text(5.3, 0.4, 'rho');
