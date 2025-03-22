
compute_optimal_gains()

function optimal_gains = compute_optimal_gains()

    % Define the initial state and time span
    % Perturb the first 3 elements of the initial state vector randomly
    rng(0);
    state0 = [0 0 0 0 0 0]' + [randn(1, 3) * 0.05, 0, 0, 0]';
    tol = 1e-9;
    options = odeset('RelTol', tol, 'AbsTol', tol);

    % Define the inertia matrix
    I1 = 1500;
    I2 = 1000;
    I3 = 800;
    I = diag([I1 I2 I3]);

    % Define the time space and coefficients for the optimal maneuver
    tspan = [0 120];
    tspace = [0 40 80 120];
    u_coeffs = zeros(3,4);
    for i = 1:3
        [a1, a2, a3, a4] = opt1axis(0, 0, tspace(i), tspace(i+1));
        u_coeffs(i,:) = [a1 a2 a3 a4];
    end

    % Define the objective function to minimize
    objective_function = @(gains) compute_error(gains, tspan, state0, tspace, u_coeffs, I, options);

    % Set initial guess for the gains
    % initial_guess = [1e-5 1e-5];
    initial_guess = 4e-2;

    % Use fminsearch to find the optimal gains
    optimal_gains = fminsearch(objective_function, initial_guess);

end

function error = compute_error(gains, tspan, state0, tspace, u_coeffs, I, options)

    % Extract the gains
    wn = gains(1);
    kp = wn^2 * I;
    kd = 2 * I * wn;

    % Define the closed-loop equations of motion
    eom_CL = @(t, state) dCL(tspace, u_coeffs, t, state, kp, kd, I);

    % Solve the ODE
    [~, state_CL] = ode45(eom_CL, tspan, state0, options);

    % Compute the error as the norm of the difference between the final state and the identity
    error = dot(state_CL(end, 1:6), state_CL(end, 1:6));

end

function [a1, a2, a3, a4] = opt1axis(theta0, thetaf, t0, tf)

    dt = tf - t0;
    a1 = theta0;
    a2 = 0;
    a3 = 3 * (thetaf - theta0) / dt^2;
    a4 = -2 * (thetaf - theta0) / dt^3;

end

function dstate = dCL(tspace, u_coeffs, time, state, kp, kd, I)

    % Define the control input
    u_ref = u_nominal(tspace, u_coeffs, time, I);
    angle_ref = theta_nominal(tspace, u_coeffs, time);
    dangle_ref = theta_dot_nominal(tspace, u_coeffs, time);

    % Extract the state variables
    yaw = state(1);
    pitch = state(2);
    roll = state(3);
    body_rates = state(4:6);

    % Define the state-space matrices
    A = zeros(6,6);
    B = zeros(6,3);
    A(1:3,4:6) = (1/cos(pitch)) * [0 sin(roll) cos(roll);
                                   0 cos(roll)*cos(pitch) -sin(roll)*cos(pitch);
                                   cos(pitch) sin(roll)*sin(pitch) cos(roll)*sin(pitch)];
    B(4:6,:) = inv(I);

    % Compute the control input
    ut = u_ref - kp * (state(1:3) - angle_ref) - kd * (state(4:6) - dangle_ref);

    % Compute the state derivative
    dstate = A * state + B * ut;

end

function ut = u_nominal(tspace, u_coeffs, time, I)

    if time <= tspace(2)
        ut = I * [1; 0; 0] * polyval([6 2] .* fliplr(u_coeffs(1,3:4)), time);
    elseif time <= tspace(3)
        ut = I * [0; 1; 0] * polyval([6 2] .* fliplr(u_coeffs(2,3:4)), time - tspace(2));
    elseif time <= tspace(4)
        ut = I * [0; 0; 1] * polyval([6 2] .* fliplr(u_coeffs(3,3:4)), time - tspace(3));
    else
        ut = [0; 0; 0];
    end

end

function theta_ref = theta_nominal(tspace, u_coeffs, time)

    if time <= tspace(2)
        theta_ref = [0; 0; 1] .* polyval(fliplr(u_coeffs(1,:)), time) ...
                  + [0; 1; 0] .* polyval(fliplr(u_coeffs(2,:)), 0) ...
                  + [1; 0; 0] .* polyval(fliplr(u_coeffs(3,:)), 0);
    elseif time <= tspace(3)
        theta_ref = [0; 1; 0] .* polyval(fliplr(u_coeffs(2,:)), time - tspace(2)) ...
                  + [1; 0; 0] .* polyval(fliplr(u_coeffs(3,:)), 0) ...
                  + [0; 0; 1] .* polyval(fliplr(u_coeffs(1,:)), tspace(2));
    elseif time <= tspace(4)
        theta_ref = [1; 0; 0] .* polyval(fliplr(u_coeffs(3,:)), time - tspace(3)) ...
                  + [0; 1; 0] .* polyval(fliplr(u_coeffs(2,:)), tspace(3) - tspace(2)) ...
                  + [0; 0; 1] .* polyval(fliplr(u_coeffs(1,:)), tspace(2));
    else
        theta_ref = [0; 0; 0];
    end

end

function theta_dot_ref = theta_dot_nominal(tspace, u_coeffs, time)

    if time <= tspace(2)
        theta_dot_ref = [1; 0; 0] * polyval(polyder(fliplr(u_coeffs(1,:))), time) ...
                      + [0; 1; 0] * polyval(polyder(fliplr(u_coeffs(2,:))), 0) ...
                      + [0; 0; 1] * polyval(polyder(fliplr(u_coeffs(3,:))), 0);
    elseif time <= tspace(3)
        theta_dot_ref = [0; 1; 0] * polyval(polyder(fliplr(u_coeffs(2,:))), time - tspace(2)) ...
                      + [0; 0; 1] * polyval(polyder(fliplr(u_coeffs(3,:))), 0) ...
                      + [1; 0; 0] * polyval(polyder(fliplr(u_coeffs(1,:))), tspace(2));
    elseif time <= tspace(4)
        theta_dot_ref = [0; 0; 1] * polyval(polyder(fliplr(u_coeffs(3,:))), time - tspace(3)) ...
                      + [0; 1; 0] * polyval(polyder(fliplr(u_coeffs(2,:))), tspace(3) - tspace(2)) ...
                      + [1; 0; 0] * polyval(polyder(fliplr(u_coeffs(1,:))), tspace(2));
    else
        theta_dot_ref = [0; 0; 0];
    end

end