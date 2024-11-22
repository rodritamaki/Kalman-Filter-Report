% This script performs the first experiment carried out for the report. The
% user is allowed to add several R values to observe the implications of
% this parameter in the execution.
%
% Author:
%- Rodrigo Tamaki

clc;
clear;
close;

% Input Parameters
true_x = randn;                                                            % True value of the constant
Q = 1e-5;                                                                  % Process noise covariance
R_values = [0.001, 0.01, 0.1, 1];                                          % Different measurement noise covariance values
num_meas = 50;                                                             % Number of measurements

% Generate random measurements (same across all R for consistency)
z = true_x + sqrt(R_values(2)) * randn(1, num_meas);                       

% Arrays to store results
all_estimates = zeros(numel(R_values), num_meas);                          % Store estimates for each R
all_covariances = zeros(numel(R_values), num_meas);                        % Store covariances for each R

% Kalman Filter for different R values
for r_idx = 1:numel(R_values)
    R = R_values(r_idx);                                                   % Current measurement noise covariance
    sigma_R = sqrt(R);                                                     % Standard deviation of measurement noise
    
    % Initial conditions
    x_hat = 0;                                                             % Initial estimate of x 
    P = 1;                                                                 % Initial erro covariance estimate
    
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
    
    % Store results for this R
    all_estimates(r_idx, :) = x_estimates;
    all_covariances(r_idx, :) = P_values;
end

%% Graphical Representation 
% Simulation results for Different R values
figure;
for r_idx = 1:numel(R_values)
    subplot(2, 2, r_idx);
    plot(1:num_meas, z, 'kx', 'DisplayName', 'Measurements');
    hold on;
    plot(1:num_meas, all_estimates(r_idx, :), 'r-', 'LineWidth', 1.5, ...
         'DisplayName', 'Kalman Estimate');
    yline(true_x, 'b--', 'DisplayName', 'True Value');
    legend;
    title(sprintf('Kalman Filter Estimate (R = %.3f)', R_values(r_idx)), ...
          'Interpreter', 'latex');
    xlabel('Iteration', 'Interpreter', 'latex');
    ylabel('Voltage', 'Interpreter', 'latex');
    grid on;
end

% Variance Convergence for Different R values
figure;
for r_idx = 1:numel(R_values)
    plot(1:num_meas, all_covariances(r_idx, :), 'LineWidth', 1.5, ...
         'DisplayName', sprintf('R = %.3f', R_values(r_idx)));
    hold on;
end
title('Variance Convergence for Different R values', 'Interpreter', 'latex');
xlabel('Iteration', 'Interpreter', 'latex');
ylabel('Covariance (Voltage$^2$)', 'Interpreter', 'latex');
legend show;
grid on;