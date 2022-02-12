%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fourth Order Butterworth Filter Design
% Fraction and Integer Bit Selection
% Embedded Scientific Computing 4450:410
% The University of Akron
% Wilson Woods
% 12.8.2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script calculates integer and fraction bits for the fourth-order
% Butterworth filter described in the report "Fourth-Order Butterworth
% Filter Design" by Wilson Woods. Methodolgy for calculating integer and
% fraction bits for individual signals is outlined in the report. The
% signals themselves are determined based on the parallel-form state-
% space representation of the filter. Refer to the block diagram contained
% in the report for signal locations and relationships.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Fraction Bit Calculations
%
% The number of fraction bits required for a given signal are found
% through the following steps:
%
% 1.)   Determine a state-space representation that maps the signal's
%       quantization error to the filter output using a modified B
%       vector shown below
%
% 2.)   Convert the state-space represenation into a transfer function,
%       and find the impulse response of that transfer function
%
% 3.)   Take the infinite sum of the absolute value of the impulse
%       response to determine its L-1 norm
%
% 4.)   Calculate fraction bits from the L-1 norm
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

format short

A = [3765/4096 647/4096 0 0; -647/4096 3765/4096 0 0; 0 0 861/1024 123/2048; 0 0 -123/2048 861/1024];
B = [-2533/481; -1027/317; -3250/697; -13945/1024];
C = [107/1024 73/2048 -53/512 -7/512];
D = 1/4096;

%%%%%%%%%%%%%%%%%%%%%%%% Fraction Bit Calculation %%%%%%%%%%%%%%%%%%%%%%%%%
B1 = [1; 0; 0; 0];
B2 = [0; 1; 0; 0];
B3 = [0; 0; 1; 0];
B4 = [0; 0; 0; 1];

% transfer function from quantization error at input to the output
% G_deltaU->Y
[num0, den0] = ss2tf(A, B, C, D);
[r0, p0, k0] = residue(num0, den0);

% the impulse (time) response of transfer function
N = 500;
gu = zeros(1, N);
d = zeros(1, N);
d(1) = 1;
u = ones(1, N);
u(1) = 0;
for k = 1 : N
    gu(k) = k0 * d(k) + ((2 * abs(r0(1)) * (abs(p0(1))^(k - 2)) * cos(angle(r0(1)) + (k - 2) * angle(p0(1)))) + (2 * abs(r0(3)) * (abs(p0(3))^(k - 2)) * cos(angle(r0(3)) + (k - 2) * angle(p0(3))))) * u(k);
end

% the l1 norm of the impulse reponse
gu_l1_norm = sum(abs(gu))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% transfer function for error in state update equation 1 to output
% G_e1->Y
[num1, den1] = ss2tf(A, B1, C, D);
[r1, p1, k1] = residue(num1, den1);

g1 = zeros(1, N);
for k = 1 : N
    g1(k) = k1 * d(k) + ((2 * abs(r1(1)) * (abs(p1(1))^(k - 2)) * cos(angle(r1(1)) + (k - 2) * angle(p1(1)))) + (2 * abs(r1(3)) * (abs(p1(3))^(k - 2)) * cos(angle(r1(3)) + (k - 2) * angle(p1(3))))) * u(k);
end

% the l1 norm of the impulse reponse
g1_l1_norm = sum(abs(g1))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% transfer function for error in state update equation 2 to output
% G_e2->Y
[num2, den2] = ss2tf(A, B2, C, D);
[r2, p2, k2] = residue(num2, den2);

g2 = zeros(1, N);
for k = 1 : N
    g2(k) = k2 * d(k) + ((2 * abs(r2(1)) * (abs(p2(1))^(k - 2)) * cos(angle(r2(1)) + (k - 2) * angle(p2(1)))) + (2 * abs(r2(3)) * (abs(p2(3))^(k - 2)) * cos(angle(r2(3)) + (k - 2) * angle(p2(3))))) * u(k);
end

% the l1 norm of the impulse reponse
g2_l1_norm = sum(abs(g2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% transfer function for error in state update equation 3 to output
% G_e3->Y
[num3, den3] = ss2tf(A, B3, C, D);
[r3, p3, k3] = residue(num3, den3);

g3 = zeros(1, N);
for k = 1 : N
    g3(k) = k3 * d(k) + ((2 * abs(r3(1)) * (abs(p3(1))^(k - 2)) * cos(angle(r3(1)) + (k - 2) * angle(p3(1)))) + (2 * abs(r3(3)) * (abs(p3(3))^(k - 2)) * cos(angle(r3(3)) + (k - 2) * angle(p3(3))))) * u(k);
end

% the l1 norm of the impulse reponse
g3_l1_norm = sum(abs(g3));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% transfer function for error in state update equation 4 to output
% G_e3->Y
[num4, den4] = ss2tf(A, B4, C, D);
[r4, p4, k4] = residue(num4, den4);

g4 = zeros(1, N);
for k = 1 : N
    g4(k) = k4 * d(k) + ((2 * abs(r4(1)) * (abs(p4(1))^(k - 2)) * cos(angle(r4(1)) + (k - 2) * angle(p4(1)))) + (2 * abs(r4(3)) * (abs(p4(3))^(k - 2)) * cos(angle(r4(3)) + (k - 2) * angle(p4(3))))) * u(k);
end

% the l1 norm of the impulse reponse
g4_l1_norm = sum(abs(g4));

% we have 2^-7 total allowable error
% simple solution: divide allowable error evenly among 6 signals
M = 2^-7 / 6;

% use L-1 norms to calculate fraction bits (yields a fractional number)
fu_ = -log(M / gu_l1_norm) / log(2);
fy_ = log(5 / M) / log(2);
f1_ = -log(M / (3 * g1_l1_norm)) / log(2);
f2_ = -log(M / (3 * g2_l1_norm)) / log(2);
f3_ = -log(M / (3 * g3_l1_norm)) / log(2);
f4_ = -log(M / (3 * g4_l1_norm)) / log(2);

% take ceiling to give integer value of minimum bits required
fu = ceil(fu_);
fy = ceil(fy_);
f1 = ceil(f1_);
f2 = ceil(f2_);
f3 = ceil(f3_);
f4 = ceil(f4_);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Integer Bit Calculation %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% L-infinity norm of input u
u_l_inf = 2

r1b = u_l_inf * abs(B(1))
r2b = u_l_inf * abs(B(2))
r3b = u_l_inf * abs(B(3))
r4b = u_l_inf * abs(B(4))

% L-infinity norm of output y
y_l_inf = gu_l1_norm * u_l_inf

% c vectors to isolate our signals
C1 = [1 0 0 0];
C2 = [0 1 0 0];
C3 = [0 0 1 0];
C4 = [0 0 0 1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Signals: r1, r11, r21, r1c, r1s

% ss2tf(): state-space representation to transfer function 
[num11, den11] = ss2tf(A, B, C1, D)

% poles, zeros, and residue to assemble transfer function
[r11, p11, k11] = residue(num11, den11)

gu_x1 = zeros(1, N);

% impulse response
for k = 1 : N
    gu_x1(k) = k11 * d(k) + (2 * abs(r11(1)) * (abs(p11(1))^(k - 2)) * cos(angle(r11(1)) + (k - 2) * angle(p11(1)))) * u(k);
end

% the l1 norm of the impulse reponse
gu_x1_l1_norm = sum(abs(gu_x1))
 
% L-infinity norms
x1_l_inf = gu_x1_l1_norm * u_l_inf
r11_l_inf = x1_l_inf * abs(A(1, 1))
r21_l_inf = x1_l_inf * abs(A(2, 1))
r1c_l_inf = x1_l_inf * abs(C(1))

% state space representation to find r1s signal integer bits
A_r1s = [3765/4096 647/4096; -647/4096 3765/4096];
B_r1s = [-2533/481; -1027/317];
C_r1s = [0 647/4096];
D_r1s = -2533/481;

[num_r1s, den_r1s] = ss2tf(A_r1s, B_r1s, C_r1s, D_r1s)
[r_r1s, p_r1s, k_r1s] = residue(num_r1s, den_r1s)

g_r1s = zeros(1, N);

% impulse response
for k = 1 : N
    g_r1s(k) = k_r1s * d(k) + (2 * abs(r_r1s(1)) * (abs(p_r1s(1))^(k - 2)) * cos(angle(r_r1s(1)) + (k - 2) * angle(p_r1s(1)))) * u(k);
end

g_r1s_l1_norm = sum(abs(g_r1s))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Signals: r2, r22, r12, r2c, r2s
[num22, den22] = ss2tf(A, B, C2, D)
[r22, p22, k22] = residue(num22, den22)

gu_x2 = zeros(1, N);

% impulse response
for k = 1 : N
    gu_x2(k) = k22 * d(k) + (2 * abs(r22(1)) * (abs(p22(1))^(k - 2)) * cos(angle(r22(1)) + (k - 2) * angle(p22(1)))) * u(k);
end

% the l1 norm of the impulse reponse
gu_x2_l1_norm = sum(abs(gu_x2))

x2_l_inf = gu_x2_l1_norm * u_l_inf
r22_l_inf = x2_l_inf * abs(A(2, 2))
r12_l_inf = x2_l_inf * abs(A(1, 2))
r2c_l_inf = x2_l_inf * abs(C(2))

% state space representation for r2s signal integer bits
A_r2s = A_r1s;
B_r2s = B_r1s;
C_r2s = [ -647/4096 0 ];
D_r2s = -1027/317;

[num_r2s, den_r2s] = ss2tf(A_r2s, B_r2s, C_r2s, D_r2s)
[r_r2s, p_r2s, k_r2s] = residue(num_r2s, den_r2s)

g_r2s = zeros(1, N);

% impulse response
for k = 1 : N
    g_r2s(k) = k_r2s * d(k) + (2 * abs(r_r2s(1)) * (abs(p_r2s(1))^(k - 2)) * cos(angle(r_r2s(1)) + (k - 2) * angle(p_r2s(1)))) * u(k);
end

g_r2s_l1_norm = sum(abs(g_r2s))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Signals: r3, r33, r43, r3c, r3s
[num33, den33] = ss2tf(A, B, C3, D)
[r33, p33, k33] = residue(num33, den33)

gu_x3 = zeros(1, N);

% impulse response
for k = 1 : N
    gu_x3(k) = k33 * d(k) + (2 * abs(r33(3)) * (abs(p33(3))^(k - 2)) * cos(angle(r33(3)) + (k - 2) * angle(p33(3)))) * u(k);
end

% the l1 norm of the impulse reponse
gu_x3_l1_norm = sum(abs(gu_x3))

x3_l_inf = gu_x3_l1_norm * u_l_inf
r33_l_inf = x3_l_inf * abs(A(3, 3))
r43_l_inf = x3_l_inf * abs(A(4, 3))
r3c_l_inf = x3_l_inf * abs(C(3))

A_r3s = [ 861/1024 123/2048 ; -123/2048 861/1024 ];
B_r3s = [ -3250/697 ; -13945/1024 ];
C_r3s = [ 0 123/2048 ];
D_r3s = -3250/697;

[num_r3s, den_r3s] = ss2tf(A_r3s, B_r3s, C_r3s, D_r3s)
[r_r3s, p_r3s, k_r3s] = residue(num_r3s, den_r3s)

g_r3s = zeros(1, N);

% impulse response
for k = 1 : N
    g_r3s(k) = k_r3s * d(k) + (2 * abs(r_r3s(1)) * (abs(p_r3s(1))^(k - 2)) * cos(angle(r_r3s(1)) + (k - 2) * angle(p_r3s(1)))) * u(k);
end

g_r3s_l1_norm = sum(abs(g_r3s))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Signals: r4, r44, r34, r4c, r4s
[num44, den44] = ss2tf(A, B, C4, D)
[r44, p44, k44] = residue(num44, den44)

gu_x4 = zeros(1, N);

% impulse response
for k = 1 : N
    gu_x4(k) = k44 * d(k) + (2 * abs(r44(3)) * (abs(p44(3))^(k - 2)) * cos(angle(r44(3)) + (k - 2) * angle(p44(3)))) * u(k);
end

% the l1 norm of the impulse reponse
gu_x4_l1_norm = sum(abs(gu_x4))


x4_l_inf = gu_x4_l1_norm * u_l_inf
r44_l_inf = x4_l_inf * abs(A(4, 4))
r34_l_inf = x4_l_inf * abs(A(3, 4))
r4c_l_inf = x4_l_inf * abs(C(4))

A_r4s = A_r3s;
B_r4s = B_r3s;
C_r4s = [ -123/2048 0 ];
D_r4s = -13945/1024;

[num_r4s, den_r4s] = ss2tf(A_r4s, B_r4s, C_r4s, D_r4s)
[r_r4s, p_r4s, k_r4s] = residue(num_r4s, den_r4s)

g_r4s = zeros(1, N);

% impulse response
for k = 1 : N
    g_r4s(k) = k_r4s * d(k) + (2 * abs(r_r4s(1)) * (abs(p_r4s(1))^(k - 2)) * cos(angle(r_r4s(1)) + (k - 2) * angle(p_r4s(1)))) * u(k);
end

g_r4s_l1_norm = sum(abs(g_r4s))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Signal: rp1
A_rp1 = A_r1s;
B_rp1 = B_r1s;
C_rp1 = [ 107/1024 73/2048 ];
D_rp1 = 0;

[num_rp1, den_rp1] = ss2tf(A_rp1, B_rp1, C_rp1, D_rp1)
[r_rp1, p_rp1, k_rp1] = residue(num_rp1, den_rp1)
k_rp1 = 0;
g_rp1 = zeros(1, N);

% impulse response
for k = 1 : N
    g_rp1(k) = k_rp1 * d(k) + (2 * abs(r_rp1(1)) * (abs(p_rp1(1))^(k - 2)) * cos(angle(r_rp1(1)) + (k - 2) * angle(p_rp1(1)))) * u(k);
end

g_rp1_l1_norm = sum(abs(g_rp1))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Signal: rp3
A_rp3 = A_r3s;
B_rp3 = B_r1s;
C_rp3 = [ 107/1024 73/2048 ];
D_rp3 = 0;

[num_rp3, den_rp3] = ss2tf(A_rp3, B_rp3, C_rp3, D_rp3)
[r_rp3, p_rp3, k_rp3] = residue(num_rp3, den_rp3)
k_rp3 = 0;
g_rp3 = zeros(1, N);

% impulse response
for k = 1 : N
    g_rp3(k) = k_rp3 * d(k) + (2 * abs(r_rp3(1)) * (abs(p_rp3(1))^(k - 2)) * cos(angle(r_rp3(1)) + (k - 2) * angle(p_rp3(1)))) * u(k);
end

g_rp3_l1_norm = sum(abs(g_rp3))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end of file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
