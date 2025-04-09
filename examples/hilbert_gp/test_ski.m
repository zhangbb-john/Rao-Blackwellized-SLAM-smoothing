% Clear workspace
clear; clc;

% Generate synthetic data
n = 100; % Number of data points
x = linspace(-5, 5, n)'; % Input points
y = sin(x) + 0.1 * randn(n, 1); % Noisy observations

% Define the kernel function (squared exponential / RBF)
kernel = @(x1, x2, l) exp(-0.5 * (x1 - x2).^2 / l^2);

% Kernel hyperparameters
l = 1.0; % Length scale

% SKI: Create a grid
m = 20; % Number of grid points
x_grid = linspace(-5, 5, m)'; % Grid points

% Compute the kernel matrix on the grid
K_grid = zeros(m, m);
for i = 1:m
    for j = 1:m
        K_grid(i, j) = kernel(x_grid(i), x_grid(j), l);
    end
end

% Interpolate the kernel matrix to the data points
K_ski = zeros(n, n);
for i = 1:n
    for j = 1:n
        % Find the nearest grid points for interpolation
        idx_i = interp1(x_grid, 1:m, x(i), 'nearest', 'extrap');
        idx_j = interp1(x_grid, 1:m, x(j), 'nearest', 'extrap');
        K_ski(i, j) = K_grid(idx_i, idx_j);
    end
end

% Add a small noise term to the diagonal for numerical stability
noise = 1e-5;
K_ski = K_ski + noise * eye(n);

% Perform regression (predict mean)
alpha = K_ski \ y; % Solve for weights
y_pred = K_ski * alpha; % Predictions

% Plot the results
figure;
plot(x, y, 'k.', 'MarkerSize', 10); hold on; % Original data
plot(x, y_pred, 'r-', 'LineWidth', 2); % SKI predictions
plot(x, sin(x), 'b--', 'LineWidth', 2); % True function
legend('Noisy Data', 'SKI Prediction', 'True Function');
title('SKI Approximation');
xlabel('x');
ylabel('y');
grid on;



% Example: Scalable Kernel Interpolation (SKI) using Random Fourier Features for GP

% Generate synthetic data (1D input, output)
n = 100;  % Number of data points
x = linspace(0, 10, n)';  % Input data points
y = sin(x) + 0.2 * randn(size(x));  % Output (with some noise)

% Define the RBF kernel parameters
length_scale = 1.0;  % Length scale for RBF kernel
sigma_f = 1.0;       % Amplitude of the kernel

% Random Fourier Features approximation of the RBF kernel
num_features = 50;   % Number of random features
omega = randn(num_features, 1) / length_scale;  % Random frequencies
b = 2 * pi * rand(num_features, 1);  % Random phase shifts

% Generate random features (Fourier features)
phi_x = sqrt(2 / num_features) * cos(x * omega' + b');  % Random feature mapping

% Compute kernel matrix K using random Fourier features
K_approx = phi_x * phi_x';  % Approximated kernel matrix using random features

% Compute the weights for GP regression
alpha = pinv(K_approx + 1e-6 * eye(n)) * y;  % Regularize the matrix inversion

% Predict using the Gaussian Process model
y_pred = phi_x' * alpha;  % Prediction based on the approximated kernel

% Visualize the results
figure;
hold on;
plot(x, y, 'bo', 'MarkerFaceColor', 'b', 'DisplayName', 'True Data');
plot(x, y_pred, 'r-', 'LineWidth', 2, 'DisplayName', 'GP Prediction (SKI)');
xlabel('Input (x)');
ylabel('Output (y)');
title('Gaussian Process Regression with Scalable Kernel Interpolation (SKI)');
legend;
hold off;
