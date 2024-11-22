% This script performs the third experiment carried out for the report. The
% user is allowed to add several initialP values to observe the implications of
% this parameter in the execution.
%
% Author:
%- Rodrigo Tamaki

clc;
clear;
close;

% Input Parameters
true_x = randn;                                                            % True value of the constant
R = 0.01;                                                                  % Fixed measurement noise covariance
Q = 1e-5;                                                                  % Fixed process noise covariance
P_values_init = [0, 0.01, 1, 100];                                         % Different initial P values to test
num_meas = 50;                                                             % Number of measurements

% Generate random measurements
z = true_x + sqrt(R) * randn(1, num_meas);                                 % Simulated measurements with fixed R

% Arrays to store results
all_estimates_P = zeros(numel(P_values_init), num_meas);                   % Store estimates for each initial P
all_covariances_P = zeros(numel(P_values_init), num_meas);                 % Store covariances for each initial P

% Kalman Filter for different initial P values
for p_idx = 1:numel(P_values_init)
    P_init = P_values_init(p_idx);                                         % Current initial P value
    
    % Initial conditions
    x_hat = 0;                                                             % Initial estimate of x
    P = P_init;                                                            % Initial error covariance estimate
    
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
    
    % Store results for this P_init
    all_estimates_P(p_idx, :) = x_estimates;
    all_covariances_P(p_idx, :) = P_values;
end

%% Graphical Representation
% Simulation results for different initial P values
figure;
for p_idx = 1:numel(P_values_init)
    subplot(2, 2, p_idx);
    plot(1:num_meas, z, 'kx', 'DisplayName', 'Measurements');
    hold on;
    plot(1:num_meas, all_estimates_P(p_idx, :), 'r-', 'LineWidth', 1.5, ...
         'DisplayName', 'Kalman Estimate');
    yline(true_x, 'b--', 'DisplayName', 'True Value');
    legend;
    title(sprintf('Kalman Filter Estimate (Po = %.2f)', P_values_init(p_idx)), ...
          'Interpreter', 'latex');
    xlabel('Iteration', 'Interpreter', 'latex');
    ylabel('Voltage', 'Interpreter', 'latex');
    grid on;
end

% Variance Convergence for Different Initial P values
figure;
for p_idx = 1:numel(P_values_init)
    plot(1:num_meas, all_covariances_P(p_idx, :), 'LineWidth', 1.5, ...
         'DisplayName', sprintf('Po = %.2f', P_values_init(p_idx)));
    hold on;
end
title('Variance Convergence for Different Initial P Values', 'Interpreter', 'latex');
xlabel('Iteration', 'Interpreter', 'latex');
ylabel('Covariance (Voltage$^2$)', 'Interpreter', 'latex');
legend show;
grid on;