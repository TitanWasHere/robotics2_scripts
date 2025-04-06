syms q1 q2 real
syms q1d q2d real
syms q1dd q2dd real
syms l1 l2 real
syms a1 a2 real
syms alpha1 alpha2 real
syms m1 m2 real
syms rcx1 rcy1 rcz1 real
syms rcx2 rcy2 rcz2 real
syms Ic1xx Ic1xy Ic1xz Ic1yx Ic1yy Ic1yz Ic1zx Ic1zy Ic1zz real
syms Ic2xx Ic2xy Ic2xz Ic2yx Ic2yy Ic2yz Ic2zx Ic2zy Ic2zz real
syms g real


g = [0; 0; -9.81]; % gravity vector TO ADJUST, if + then my axis in direction of g
r1 = [rcx1; rcy1; rcz1];
r2 = [rcx2; rcy2; rcz2];
I1 = [Ic1xx Ic1xy Ic1xz; Ic1yx Ic1yy Ic1yz; Ic1zx Ic1zy Ic1zz];
I2 = [Ic2xx Ic2xy Ic2xz; Ic2yx Ic2yy Ic2yz; Ic2zx Ic2zy Ic2zz];
I = [I1 I2];
r = [r1 r2];
m = [m1; m2];

DH = [alpha1 a1 l1 q1; alpha2 a2 l2 q2];
q = [q1; q2];
qd = [q1d; q2d];
qdd = [q1dd; q2dd];

[DK, A] = DHMatrix(DH);
p_e = DK(1:3, 4);

J = simplify(jacobian(p_e, q));

[v, w, vc, T] = movingFramesAlgorithm(DH, qd, r, ["r", "r"]);

T_tot = 0;
for i = 1:size(DH,1)
    T_tot = T_tot + KonigTheorem(m(i), vc{i}, w{i}, I{i});
end

T_tot = simplify(T_tot);

M = inertia_from_kinetic(T, qd);
M = simplify(M);

if is_inertia_matrix(M)
    disp('The matrix is a valid inertia matrix.');
else
    error('The matrix is not a valid inertia matrix.');
end

[c, C] = inertia_matrix_to_coriolis(M, q, dq);

c = simplify(c);
C = simplify(C);

% compute the potential energy
U = 0;
syms U01 U02 real
U0 = [U01 U02];
for i = 1:size(DH,1)
    U = U + U0(i) - m(i) * g' * (A{i}(1:3, 1:3) * r(i));
end
U = simplify(U);

% Now we can compute the gravity vector
g = g_from_potential_energy(U, q);
g = simplify(g);

tau = M * qdd + C * qd + g;