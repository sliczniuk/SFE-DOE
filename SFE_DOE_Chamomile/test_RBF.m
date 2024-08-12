startup;
delete(gcp('nocreate'));

%addpath('C:\Dev\casadi-3.6.3-windows64-matlab2018b');
addpath('\\home.org.aalto.fi\sliczno1\data\Documents\casadi-3.6.3-windows64-matlab2018b');
import casadi.*

%%
% Define the dataset
[x_data, y_data] = meshgrid(linspace(0, 2*pi, 100), linspace(0, 2*pi, 100));
z_data = sin(x_data) + cos(y_data);
x_data = x_data(:);
y_data = y_data(:);
z_data = z_data(:);

% Number of RBF kernels specified by the user
N = 3; % For example, you can change this to any number of RBF kernels

% Create an Opti instance
opti = Opti();

% Parameters for N RBFs
cx = opti.variable(N, 1); % Centers of the RBFs in x
cy = opti.variable(N, 1); % Centers of the RBFs in y
w = opti.variable(N, 1);  % Weights of the RBFs
sx = opti.variable(N, 1); % Widths of the RBFs in x (standard deviations)
sy = opti.variable(N, 1); % Widths of the RBFs in y (standard deviations)
b = opti.variable();      % Bias term

% RBF function
rbf = @(x, y, cx, cy, sx, sy) exp(-((x - cx).^2) / (2 * sx^2) - ((y - cy).^2) / (2 * sy^2));

% Model prediction
z_pred = @(x, y) b;
for i = 1:N
    z_pred = @(x, y) z_pred(x, y) + w(i) * rbf(x, y, cx(i), cy(i), sx(i), sy(i));
end

% Objective: minimize the mean squared error
mse = sum((z_pred(x_data, y_data) - z_data).^2) / length(z_data);
opti.minimize(mse);

% Set initial guesses
initial_cx = linspace(0, 2*pi, N)';
initial_cy = linspace(0, 2*pi, N)';
initial_w = ones(N, 1);
initial_sx = ones(N, 1);
initial_sy = ones(N, 1);

opti.set_initial(cx, initial_cx);
opti.set_initial(cy, initial_cy);
opti.set_initial(w, initial_w);
opti.set_initial(sx, initial_sx);
opti.set_initial(sy, initial_sy);
opti.set_initial(b, 0);

% Solver options
opts = struct;
opts.ipopt.print_level = 5;
opti.solver('ipopt', opts);

% Solve the problem
sol = opti.solve();

% Extract the optimal parameters
cx_opt = sol.value(cx);
cy_opt = sol.value(cy);
w_opt = sol.value(w);
sx_opt = sol.value(sx);
sy_opt = sol.value(sy);
b_opt = sol.value(b);

% Display the optimal parameters
disp('Optimal parameters:');
disp(['cx = ', num2str(cx_opt')]);
disp(['cy = ', num2str(cy_opt')]);
disp(['w = ', num2str(w_opt')]);
disp(['sx = ', num2str(sx_opt')]);
disp(['sy = ', num2str(sy_opt')]);
disp(['b = ', num2str(b_opt)]);

% Plot the results
x_plot = linspace(0, 2*pi, 100);
y_plot = linspace(0, 2*pi, 100);
[X_plot, Y_plot] = meshgrid(x_plot, y_plot);
Z_plot = arrayfun(@(x, y) sol.value(z_pred(x, y)), X_plot, Y_plot);

%%
figure;
surf(x_plot, y_plot, Z_plot, 'linestyle','none');
hold on;
surf(x_plot, y_plot, sin(X_plot) + cos(Y_plot), 'FaceAlpha', 0.5, 'linestyle','none');
legend('RBF Approximation', 'Original Function');
xlabel('x');
ylabel('y');
zlabel('z');
title(['RBF Approximation of sin(x) + cos(y) with ', num2str(N), ' Kernels']);
grid on;