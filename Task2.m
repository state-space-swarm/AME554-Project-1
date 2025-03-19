% This script will compute the solutions to Task 2 of AME554 Project 1

% Task 2 directs us to complete problems 3.38 and 3.39c in the text

% Problem 3.38
%  A spacecraft has 4 attitude sensors, the first sensor is the most accurate
% Each sensor gives a unit vector, described in the body frame as well as the inertial frame

% Each vector may not actually be a unit vector, normalize as needed
% Compute 3 different estimates for the attitude using the TRIAD algorithm
% Determine the least accurate measurement from the attitude estimates

% Problem 3.39c
% Using the first 3 sensor measurements, determine the attitude using the OLAE method

%% Problem Statement

v1_b = [0.8273 0.5541 -0.0920]';
v2_b = [-0.8285 0.5522 -0.0955]';
v3_b = [0.2155 0.5522 0.8022]';
v4_b = [0.5570 -0.7442 -0.2884]';

v1_n = [-0.1517 -0.9669 0.2050]';
v2_n = [-0.8393 0.4494 -0.3044]';
v3_n = [-0.0886 -0.5856 -0.8000]';
v4_n = [0.8814 -0.0303 0.5202]';

v_body = [v1_b v2_b v3_b v4_b];
v_inertial = [v1_n v2_n v3_n v4_n];

%% Triad Algorithm

% Initialize array of attitude matrices
BNs = zeros(3,3,3);

for i = 2:4

    BN = triad_attitude(v_body(:,1), v_body(:,i), v_inertial(:,1), v_inertial(:,i));
    BNs(:,:,i-1) = BN;

end

% Compute the attitude error for each pair of attitude matrices

% Initialize array of attitude errors
errors = zeros(1,3);

for i = 1:3

    j = mod(i-2,3)+1;
    k = mod(i,3)+1;

    error = prangle(BNs(:,:,j),BNs(:,:,k));
    errors(i) = error;

end

% The error is the lowest when measurement 4 is omitted, so measurement 4 is the least accurate

%% OLAE Algorithm

% Construct the sum and difference arrays

S = v_body + v_inertial;
D = v_body - v_inertial;

S = S(:, 1:3);
D = D(:, 1:3);
S = S(:);
D = D(:);

[S, D] = olaeSD(S, D);

% Constrtuct the Weights matrix
W = [2 1 1];
W = olae_weights(W, D);

% Compute the CRP attitude description
q = (S'*W*S)\(S'*W*D);

BN_olae = rod2dcm(q');

%% Control, Drive to Identity

I1 = 1500;
I2 = 1000;
I3 = 800;

% Construct the inertia matrix
I = diag([I1 I2 I3]);

t0 = 0;
tf = 120;

% Equations of Motion
% I_c*w_dot = -w_skew*I_c*w + L_c
% x_dot = A*x + B*u

% This stuff is a start...
% % Initial conditions
% w0 = [0.1; 0.1; 0.1]; % Initial angular velocity
% q0 = [0.1; 0.1; 0.1]; % Initial CRP

% % Desired final conditions
% wf = [0; 0; 0]; % Final angular velocity
% qf = [0; 0; 0]; % Final CRP (identity attitude)

% % Time vector
% tspan = [t0 tf];

% % State vector [w; q]
% x0 = [w0; q0];

% % Control input (assuming zero external torque for simplicity)
% u = @(t) [0; 0; 0];

% % Equations of motion
% A = [zeros(3), eye(3); zeros(3), -inv(I)];
% B = [zeros(3); inv(I)];

% % Define the ODE function
% odefun = @(t, x) A*x + B*u(t);

% % Solve the ODE
% [t, x] = ode45(odefun, tspan, x0);

% % Extract the results
% w = x(:, 1:3);
% q = x(:, 4:6);

% % Plot the results
% figure;
% subplot(2,1,1);
% plot(t, w);
% title('Angular Velocity');
% xlabel('Time (s)');
% ylabel('Angular Velocity (rad/s)');
% legend('w1', 'w2', 'w3');

% subplot(2,1,2);
% plot(t, q);
% title('Classical Rodriguez Parameters');
% xlabel('Time (s)');
% ylabel('CRP');
% legend('q1', 'q2', 'q3');


function TT = triad_Tframe(v1,v2)

    v1 = v1/norm(v1);
    v2 = v2/norm(v2);

    t1 = v1;
    t2 = cross(v1,v2);
    t2 = t2/norm(t2);
    t3 = cross(t1,t2);

    TT = [t1 t2 t3];

end

function BN = triad_attitude(v1b, v2b, v1n, v2n)

    BT = triad_Tframe(v1b,v2b);
    NT = triad_Tframe(v1n,v2n);

    BN = BT*NT';

end

function PA = prangle(M1, M2)

    B_true = M1;
    B_guess = M2;

    BB = B_true*B_guess';
    PA = acosd((trace(BB)-1)/2);

end

function [S, D] = olaeSD(ss, dd)

    slen = length(ss)/3;
    S = zeros(slen*3,3);
    for i = 1:slen
        S((i-1)*3+1:i*3,:) = skew(ss((i-1)*3+1:i*3)');
    end

    D = dd;

end

function W = olae_weights(ww, dd)

    slen = length(dd)/3;
    W = zeros(slen*3);
    for i = 1:slen
        W((i-1)*3+1:i*3,(i-1)*3+1:i*3) = ww(i)*eye(3);
    end

end