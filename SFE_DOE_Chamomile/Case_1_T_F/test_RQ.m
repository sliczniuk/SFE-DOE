% Ensure CasADi is installed and available in your MATLAB environment
addpath('\\home.org.aalto.fi\sliczno1\data\Documents\casadi-3.6.3-windows64-matlab2018b');
% Ensure CasADi is installed and available in your MATLAB environment
import casadi.*

% Number of data points
N = 100;    % Number of points in the dataset
alpha_init = 1.0;  % Initial value of Rational quadratic kernel parameter

% Generate sine wave data
x_data = linspace(0, 2*pi, N)';
y_data = sin(x_data) + 0.1*randn(N, 1);  % Adding a bit of noise

% Select every third data point for the kernel (downsampling)
indices = 1:3:N;
x_selected = x_data(indices);  % Use every third point for kernel calculations

% Rational Quadratic Kernel function
function K = RationalQuadratic_kernel(x, x_prime, l, alpha)
    % x and x_prime are input vectors
    % l is the length-scale
    % alpha is the relative weighting parameter
    K = (1 + ((x - x_prime).^2) / (2 * alpha * l^2)).^(-alpha);
end

% Setting up CasADi OptiStack for optimization
opti = Opti();

% Define Rational Quadratic Kernel parameters for optimization
l = opti.variable(1);           % Length-scale parameter
alpha = opti.variable(1);        % Alpha parameter for kernel flexibility
kernel_weights = opti.variable(length(indices), 1);  % Weights for each selected data point

% Initialize the Rational Quadratic kernel matrix for the full dataset
K_full = MX(zeros(N, length(indices)));  % Only uses the selected data points

% Fill the kernel matrix using the Rational Quadratic Kernel
for i = 1:N
    for j = 1:length(indices)
        K_full(i, j) = RationalQuadratic_kernel(x_data(i), x_selected(j), l, alpha);
    end
end

% Construct the model: linear combination of kernel weights and kernel matrix
model = K_full * kernel_weights;

% Define the loss function: least squares error between the model and actual data
loss = sum((model - y_data).^2);

% Set the objective to minimize the loss
opti.minimize(loss);

% Set initial guesses for the kernel parameters and weights
opti.set_initial(l, 1);
opti.set_initial(alpha, alpha_init);
opti.set_initial(kernel_weights, ones(length(indices), 1));

% Solver options
p_opts = struct('print_time', false, 'ipopt', struct('print_level', 0));
opti.solver('ipopt', p_opts);  % Use IPOPT solver for optimization

% Solve the optimization problem
sol = opti.solve();

% Extract optimized parameters
l_opt = sol.value(l);
alpha_opt = sol.value(alpha);
kernel_weights_opt = sol.value(kernel_weights);

% Recompute the fitted model with optimized parameters
K_opt = zeros(N, length(indices));
for i = 1:N
    for j = 1:length(indices)
        K_opt(i, j) = RationalQuadratic_kernel(x_data(i), x_selected(j), l_opt, alpha_opt);
    end
end
model_opt = K_opt * kernel_weights_opt;

% Plot the original data and the fitted model
figure;
plot(x_data, y_data, 'bo', 'DisplayName', 'Data'); hold on;
plot(x_data, model_opt, 'r-', 'DisplayName', 'Fitted Model');
title('Rational Quadratic Kernel Regression (Every Third Point)');
legend;

% Display the optimized parameters
disp(['Optimized Length-Scale (l): ', num2str(l_opt)]);
disp(['Optimized Alpha: ', num2str(alpha_opt)]);
