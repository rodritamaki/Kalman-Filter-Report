% This script performs the second experiment carried out for the report. The
% user is allowed to add several Q values to observe the implications of
% this parameter in the execution.
%
% Author:
%- Rodrigo Tamaki

clc;
clear;
close;

% Input Parameters
true_x = randn;                                                            % True value of the constant
Q_values = [1e-7, 1e-5, 1e-3];                                             % Process noise covariance values to test
R = 0.01;                                                                  % Fixed measurement noise covariance
num_meas = 50;                                                             % Number of measurements

% Generate random measurements
z = true_x + sqrt(R) * randn(1, num_meas);                                 % Simulated measurements with fixed R

% Arrays to store results
all_estimates_Q = zeros(numel(Q_values), num_meas);                        % Store estimates for each Q
all_covariances_Q = zeros(numel(Q_values), num_meas);                      % Store covariances for each Q

% Kalman Filter for different Q values
for q_idx = 1:numel(Q_values)
    Q = Q_values(q_idx);                                                   % Current process noise covariance
    
    % Initial conditions
    x_hat = 0;                                                             % Initial estimate of x
    P = 1;                                                                 % Initial error covariance estimate
    
    % Arrays for storing the estimated values
    x_estimates = zeros(1, num_meas);
    P_values = zeros(1, num_meas);

    % Kalman Filter recursive algorithm
    for k = 1:num_meas
        % Time Update Equations (Predict)
        x_hat_prior = x_hat;                                               % Predicted state (constant process)
        P_prior = P + Q;                                                   % Predicted error covariance
        
        % Measurement Update Equations (Correct)
        K = P_prior / (P_prior + R);                                       % Kalman Gain
        x_hat = x_hat_prior + K * (z(k) - x_hat_prior);                    % Updated a posteriori estimate
        P = (1 - K) * P_prior;                                             % Updated error covariance
        
        % Results storage
        x_estimates(k) = x_hat;
        P_values(k) = P;
    end
    
    % Store results for this Q
    all_estimates_Q(q_idx, :) = x_estimates;
    all_covariances_Q(q_idx, :) = P_values;
end

%% Graphical Representation
% Simulation results for different Q values
figure;
for q_idx = 1:numel(Q_values)
    subplot(1, numel(Q_values), q_idx);
    plot(1:num_meas, z, 'kx', 'DisplayName', 'Measurements');
    hold on;
    plot(1:num_meas, all_estimates_Q(q_idx, :), 'r-', 'LineWidth', 1.5, ...
         'DisplayName', 'Kalman Estimate');
    yline(true_x, 'b--', 'DisplayName', 'True Value');
    legend;
    title(sprintf('Kalman Filter Estimate (Q = %.1e)', Q_values(q_idx)), ...
          'Interpreter', 'latex');
    xlabel('Iteration', 'Interpreter', 'latex');
    ylabel('Voltage', 'Interpreter', 'latex');
    grid on;
end

% Variance Convergence for Different Q values
figure;
for q_idx = 1:numel(Q_values)
    plot(1:num_meas, all_covariances_Q(q_idx, :), 'LineWidth', 1.5, ...
         'DisplayName', sprintf('Q = %.1e', Q_values(q_idx)));
    hold on;
end
title('Variance Convergence for Different Q values', 'Interpreter', 'latex');
xlabel('Iteration', 'Interpreter', 'latex');
ylabel('Covariance (Voltage$^2$)', 'Interpreter', 'latex');
legend show;
grid on;
