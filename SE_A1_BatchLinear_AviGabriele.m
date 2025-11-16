clear 
clear all
clear clc
%% --- 1. Load Data and Define Parameters ---
load('dataset1.mat');

T = 0.1;           % Sampling period [s] (assumed uniform)
K = length(t);     % Total number of full timesteps

% Variances
sigma_r2 = r_var;   % Observation noise variance (n_k) [m^2]
sigma_v2 = v_var;   % Speed noise variance [m^2/s^2]

% Process noise variance (w_k)
sigma_q2 = T^2 * sigma_v2; 

% Pre-calculate inverse variances (Information Weights)
lambda_r = 1 / sigma_r2; % Observation information weight
lambda_q = 1 / sigma_q2; % Motion information weight

%% Minimization

% Transformed observation vector
y = l - r;

% Odometry input
f = T * v;

% Define the deltas to test
deltas = [1000, 100, 10, 1];

for delta = deltas
    
    % Total number of sparse states
    M = floor(K / delta);

    % New process noise variance
    sigma_q_tilde2 = delta * sigma_q2;
    lambda_q_tilde = 1 / sigma_q_tilde2;

    % Sparse index
    idx_sparse = (1:M) * delta;

    % Sparse vectors
    y_sparse = y(idx_sparse);
    x_true_sparse = x_true(idx_sparse);
    t_sparse = t(idx_sparse);

    % Integrated odometry
    for m = 1 : M
        idx_start = (m - 1) * delta + 1;
        idx_end = m * delta;
        tilde_f(m) = sum(f(idx_start:idx_end));
    end

    % Matrices implementation
    H_sparse = sparse(M, M);
    H_sparse_q = sparse(M, M);
    H_sparse_r = sparse(M, M);
    b_sparse = sparse(M, 1);
    b_sparse_r = sparse(M, 1);
    b_sparse_q = sparse(M, 1);

    % H_sparse matrix
    diag_val = lambda_r + 2*lambda_q_tilde;
    for m = 2 : (M-1)
        H_sparse(m, m) = diag_val;

        H_sparse_r(m, m) = lambda_r;

        H_sparse_q(m, m) = 2*lambda_q_tilde;
    end

    bound_val = lambda_r + lambda_q_tilde;
    H_sparse(1, 1) = bound_val;
    H_sparse(M, M) = bound_val;

    H_sparse_r(1, 1) = lambda_r;
    H_sparse_r(M, M) = lambda_r;

    H_sparse_q(1, 1) = lambda_q_tilde;
    H_sparse_q(M, M) = lambda_q_tilde;

    out_val = -lambda_q_tilde;
    for m = 1 : (M-1)
        H_sparse(m, m+1) = out_val;
        H_sparse(m+1, m) = out_val;

        H_sparse_r(m, m+1) = 0;
        H_sparse_r(m+1, m) = 0;

        H_sparse_q(m, m+1) = out_val;
        H_sparse_q(m+1, m) = out_val;
    end

    % b_sparse vector
    b_sparse = b_sparse + lambda_r*y_sparse;
    b_sparse_r = b_sparse_r + lambda_r*y_sparse;
    b_sparse_q = b_sparse_q;

    prior = 0;
    
    b_sparse_q(1) = b_sparse_q(1) + lambda_q_tilde*(tilde_f(2)) + prior;
    b_sparse(1) = b_sparse(1) + lambda_q_tilde*(tilde_f(2));
    b_sparse_r(1) = 0;

    for m = 2:(M - 1)
        b_sparse(m) = b_sparse(m) + ...
                            lambda_q_tilde * (tilde_f(m) - tilde_f(m + 1));

        b_sparse_r(m) = b_sparse_r(m);

        b_sparse_q(m) = b_sparse_q(m) + lambda_q_tilde * (tilde_f(m) - tilde_f(m + 1)) + prior;
    end
    b_sparse(M) = b_sparse(M) + lambda_q_tilde*(tilde_f(M));
    b_sparse_q(M) = b_sparse_q(M) + lambda_q_tilde*(tilde_f(M)) + prior;

    % Solution

    % Optimal state estimate: x* = H^-1 * b'
    x_star_sparse = H_sparse \ b_sparse;
    x_star_sparse_r = H_sparse_r \ b_sparse_r;
    x_star_sparse_q = H_sparse_q \ b_sparse_q;

    % Covariance matrix
    sigma_x_sparse = full(inv(H_sparse));
    sigma_x2_sparse = diag(sigma_x_sparse);
    sigma_x_sparse = sqrt(sigma_x2_sparse);

    sigma_x_sparse_r = full(inv(H_sparse_r));
    sigma_x2_sparse_r = diag(sigma_x_sparse_r);
    sigma_x_sparse_r = sqrt(sigma_x2_sparse_r);

    sigma_x_sparse_q = full(inv(H_sparse_q));
    sigma_x2_sparse_q = diag(sigma_x_sparse_q);
    sigma_x_sparse_q = sqrt(sigma_x2_sparse_q);

    % Error and Uncertainty Envelope
    error = x_star_sparse - x_true_sparse;
    uncertainty_envelope = 3 * sigma_x_sparse;

    error_r = x_star_sparse_r - x_true_sparse;
    uncertainty_envelope_r = 3 * sigma_x_sparse_r;

    error_q = x_star_sparse_q - x_true_sparse;
    uncertainty_envelope_q = 3 * sigma_x_sparse_q;

    % --- 6. Generate Figures (Enhanced for Report Quality) ---

    % (a) Plot Error and Uncertainty Envelope
    % NOTE: Assuming 'uncertainty_envelope' is already calculated as 3 * sqrt(diag(P_sparse))
    figure('Name', sprintf('Error and Uncertainty (delta=%d)', delta), 'Position', [100 100 800 500]);
    
    % Plot the actual position error (Solid Black Line)
    plot(t_sparse, error, 'k-', 'LineWidth', 2);
    hold on;
    
    % Plot the theoretical uncertainty envelope (Dashed Red Lines)
    plot(t_sparse, uncertainty_envelope, 'r--', 'LineWidth', 1.5);
    plot(t_sparse, -uncertainty_envelope, 'r--', 'LineWidth', 1.5);
    
    % NOTE: The redundant '_r' uncertainty lines from the original code are removed for clarity.
    hold off;
    
    title(sprintf('Position Error vs. Time (\\delta=%d)', delta), 'Interpreter', 'tex', 'FontSize', 14);
    xlabel('Time [s]', 'FontSize', 12);
    ylabel('Position Error and Uncertainty [m]', 'FontSize', 12);
    legend('Error $x_{k\delta}^{*} - x_{k\delta}$', '$\pm 3\sigma_{x_{k\delta}}$', ...
           'Interpreter', 'latex', 'Location', 'SouthEast', 'FontSize', 10);
    grid on;
    box on;
    set(gca, 'FontSize', 10);

    % (b) Plot Histogram of Errors
    figure('Name', sprintf('Error Histogram (delta=%d)', delta), 'Position', [900 100 800 500]);
    
    % Calculate error statistics for the plot title
    mean_error = mean(error);
    std_error = std(error);
    
    % Plot the actual error distribution (Histogram)
    histogram(error, 'Normalization', 'pdf', 'FaceColor', [0.7 0.7 0.9], 'EdgeColor', 'k');
    hold on;
    
    % Theoretical Gaussian PDF calculated using the actual error statistics
    x_range = linspace(min(error)*1.1, max(error)*1.1, 200);
    pdf_theoretical = normpdf(x_range, mean_error, std_error);
    
    % Plot the fitted Gaussian curve (Solid Red Line)
    plot(x_range, pdf_theoretical, 'r-', 'LineWidth', 2.5);
    
    % NOTE: The redundant '_r' and '_q' Gaussian curves from the original code are removed for clarity.
    hold off;
    
    % 1. Create the title string using single backslashes for LaTeX symbols and newline
    title_text = sprintf('Histogram of Position Errors (delta=%d) \nMean: %.4f [m], StdDev: %.4f [m]', ...
                     delta, mean_error, std_error);

    % 2. Apply the title and explicitly set the Interpreter to 'latex'
    title(title_text, 'Interpreter', 'tex', 'FontSize', 14);
    
    % Axes Labels
    xlabel('Error $x_{k\delta}^{*} - x_{k\delta}$ [m]', 'Interpreter', 'latex', 'FontSize', 12); 
    ylabel('Probability Density Function', 'FontSize', 12);
    
    % Legend
    legend('Actual Error Distribution', 'Fitted Gaussian Distribution', ...
           'Location', 'NorthWest', 'FontSize', 10);
    grid on;
    box on;
    set(gca, 'FontSize', 10);
end
