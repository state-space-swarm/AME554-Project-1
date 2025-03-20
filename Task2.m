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

% Optimal Rest to Rest Maneuver
% Drive to identity using 3 sequential single axis maneuvers

clear; close all; clc;

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

% Construct the Weights matrix
W = [2 1 1];
W = olae_weights(W, D);

% Compute the CRP attitude description
q = (S'*W*S)\(S'*W*D);

BN_olae = rod2dcm(q');
[EA1, EA2, EA3] = dcm2angle(BN_olae); % Z angle, Y angle, X angle
EA = fliplr([EA1 EA2 EA3]); % X angle, Y angle, Z angle

%% Control, Drive to Identity

I1 = 1500;
I2 = 1000;
I3 = 800;

% Construct the inertia matrix
I = diag([I1 I2 I3]);
I_vec = [I1 I2 I3];

tspace = [0 40 80 120];

% Initialize maneuver array

u_coeffs = zeros(3,4);

for i = 1:3

    [a1, a2, a3, a4] = opt1axis(EA(i), 0, tspace(i), tspace(i+1));
    u_coeffs(i,:) = [a1 a2 a3 a4];

end

% Compute EA time history based on a 1-2-3 maneuver, undoing a 3-2-1 set of EA's

eom_OL = @(t, state) dOL(tspace, u_coeffs, t, state);

tspan = [0 120];
state0 = [EA1 EA2 EA3 0 0 0]';

tol = 1e-8;
options = odeset('RelTol', tol, 'AbsTol', tol);

[t, state] = ode45(eom_OL, tspan, state0, options);

% Plot results

figure;
subplot(2,1,1);
plot(t, state(:,1:3));
title('Euler Angles Time History');
xlabel('Time (s)');
ylabel('Angle (rad)');
legend('Yaw', 'Pitch', 'Roll');

% Compute the nominal control effort u(t) over the time span
u_nominal_values = zeros(length(t), 3);
for idx = 1:length(t)
    u_nominal_values(idx, :) = u_nominal(tspace, u_coeffs, t(idx));
end

subplot(2,1,2);
plot(t, u_nominal_values);
title('Nominal Control Effort Time History');
xlabel('Time (s)');
ylabel('Control Effort');
legend('u1', 'u2', 'u3');

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

function [a1, a2, a3, a4] = opt1axis(theta0, thetaf, t0, tf)

    dt = tf-t0;

    a1 = theta0;
    a2 = 0;
    a3 = 3*(thetaf-theta0)/dt^2;
    a4 = -2*(thetaf-theta0)/dt^3;

end

function ut = u_nominal(tspace, u_coeffs, time)

    if time <= tspace(2)
        ut = [1; 0; 0]*polyval([6 2].*fliplr(u_coeffs(1,3:4)), time);
    elseif time <= tspace(3)
        ut = [0; 1; 0]*polyval([6 2].*fliplr(u_coeffs(2,3:4)), time - tspace(2));
    else
        ut = [0; 0; 1]*polyval([6 2].*fliplr(u_coeffs(3,3:4)), time - tspace(3));
    end

end

function dstate = dOL(tspace, u_coeffs, time, state)

    % This may be wonky, but the state is defined here as:
    % state = [yaw pitch roll body3rate body2rate body1rate]'

    % The only control input are the body axis rates

    ut = u_nominal(tspace, u_coeffs, time);

    A = zeros(6,6);
    B = zeros(6,3);

    yaw = state(1);
    pitch = state(2);
    roll = state(3);
    body3_rate = state(4);
    body2_rate = state(5);
    body1_rate = state(6);

    A(1:3,4:6) = (1/cos(pitch))*[0 sin(roll) cos(roll);
                                 0 cos(roll)*cos(pitch) -sin(roll)*cos(pitch);
                                 cos(pitch) sin(roll)*sin(pitch) cos(roll)*sin(pitch)];
    B(4:6,:) = eye(3);

    dstate = A*state + B*ut;

end