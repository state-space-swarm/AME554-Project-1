% This script determines the optimal tspace vector that minimizes the total control effort
% of the optimal rest to rest maneuver. It graphs the surface where the x and y coordinates
% are the middle values in the tspace vector and the z coordinate is the cost function.

clear; close all; clc;

% Define the initial and final times
t0 = 0;
tf = 120;

% Define the range for the middle values in the tspace vector
t1_range = linspace(0, 120, 20);
t2_range = linspace(0, 120, 20);

% Initialize the cost function matrix
J_matrix = zeros(length(t1_range), length(t2_range));

% Define the inertia matrix
I1 = 1500;
I2 = 1000;
I3 = 800;
I = diag([I1 I2 I3]);

% Define the initial Euler angles
EA = [-2.8846 -0.3157 -1.1180];
EA1 = EA(1);
EA2 = EA(2);
EA3 = EA(3);

% Loop over the range of middle values
for i = 1:length(t1_range)
    for j = 1:length(t2_range)
        tspace = [t0 t1_range(i) t2_range(j) tf];
        
        % Initialize maneuver array
        u_coeffs = zeros(3,4);
        
        for k = 1:3
            [a1, a2, a3, a4] = opt1axis(EA(k), 0, tspace(k), tspace(k+1));
            u_coeffs(k,:) = [a1 a2 a3 a4];
        end
        
        % Compute the nominal control effort u(t) over the time span
        tspan = linspace(t0, tf, 100);
        u_nominal_values = zeros(length(tspan), 3);
        for idx = 1:length(tspan)
            u_nominal_values(idx, :) = u_nominal(tspace, u_coeffs, tspan(idx), I);
        end
        
        % Compute the total control effort
        J = sum(trapz(tspan, u_nominal_values.^2));
        J_matrix(i, j) = J;
    end
end

% Plot the surface with log scale for the z coordinate
[T1, T2] = meshgrid(t1_range, t2_range);
figure;
surf(T1, T2, log(J_matrix'));
xlabel('t1 (s)');
ylabel('t2 (s)');
zlabel('Log of Total Control Effort');
title('Optimal tspace Vector Minimization');
colorbar;

function [a1, a2, a3, a4] = opt1axis(theta0, thetaf, t0, tf)

    dt = tf-t0;

    a1 = theta0;
    a2 = 0;
    a3 = 3*(thetaf-theta0)/dt^2;
    a4 = -2*(thetaf-theta0)/dt^3;

end

function ut = u_nominal(tspace, u_coeffs, time, Imat)

    if time <= tspace(2)
        ut = Imat*[1; 0; 0]*polyval([6 2].*fliplr(u_coeffs(1,3:4)), time);
    elseif time <= tspace(3)
        ut = Imat*[0; 1; 0]*polyval([6 2].*fliplr(u_coeffs(2,3:4)), time - tspace(2));
    elseif time <= tspace(4)
        ut = Imat*[0; 0; 1]*polyval([6 2].*fliplr(u_coeffs(3,3:4)), time - tspace(3));
    else
        ut = [0; 0; 0];
    end

end