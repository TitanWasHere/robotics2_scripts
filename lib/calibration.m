function [est_params, delta_params] = calibration(Q, P, fkine_sym, q_sym, param_sym, nominal_params)
%CALIBRATION Solve kinematic parameter calibration via linear least squares
%   [est_params, delta_params] = calibration(Q, P, fkine_sym, q_sym, param_sym, nominal_params)
%   inputs:
%     Q             MxN matrix of joint configurations (each row is one q)
%     P             MxD matrix of measured end-effector positions (each row is one p)
%     fkine_sym     Dx1 symbolic direct kinematic map f(q_sym, param_sym)
%     q_sym        Nx1 symbolic vector of joint variables
%     param_sym    Kx1 symbolic vector of kinematic parameters to calibrate
%     nominal_params Kx1 numeric vector of nominal parameter values
%   outputs:
%     est_params   Kx1 numeric vector of estimated parameters (nominal + delta)
%     delta_params Kx1 numeric vector of parameter corrections

% Check dimensions
[M, N] = size(Q);
[Mp, D] = size(P);
if M~=Mp
    error('Number of Q and P samples must match');
end
K = numel(param_sym);

% Symbolic Jacobian wrt parameters
Jp_sym = jacobian(fkine_sym, param_sym); % DxK symbolic

% Substitute nominal parameters into kinematics and Jacobian
fkine_nom_sym = subs(fkine_sym, param_sym, nominal_params(:));
Jp_nom_sym    = subs(Jp_sym,    param_sym, nominal_params(:));

% Preallocate stacked matrices
Phi     = zeros(M*D, K);
DeltaP  = zeros(M*D, 1);

% Loop over measurements
for i = 1:M
    qi = Q(i, :)';      % Nx1
    pi = P(i, :)';      % Dx1
    % Evaluate nominal forward kinematics at qi
    pi_nom = double(subs(fkine_nom_sym, q_sym, qi));
    % Residual
    ri = pi - pi_nom;
    % Evaluate parameter Jacobian at qi
    Jpi = double(subs(Jp_nom_sym, q_sym, qi));
    % Stack into regression
    idx = (i-1)*D + (1:D);
    Phi(idx, :)   = Jpi;
    DeltaP(idx)   = ri;
end

%% Example of usage
%--- Example usage of calibration.m for a 2R planar robot ---
%clear; clc;

% 1) Define symbolic model
%syms l1 l2 q1 q2
%param_sym  = [l1; l2];            % parameters to calibrate
%q_sym      = [q1; q2];            % joint variables
%fkine_sym  = [ ...
%    l1*cos(q1) + l2*cos(q1+q2);   % x-position
%    l1*sin(q1) + l2*sin(q1+q2)    % y-position
%];

% 2) Nominal link-lengths
%nominal_params = [1; 1];  % \hat l1 = \hat l2 = 1 m

% 3) Measured data (from the midterm exercise)
%Q = [ ...
%     0,        0;
%     pi/2,     0;
%     pi/4,   -pi/4;
%     0,     pi/4
%];                               % M×2 matrix of joint samples
%P = [ ...
%     2.0000, 0.0000;
%     0.0000, 2.0000;
%     1.6925, 0.7425;
%     1.7218, 0.6718
%];                               % M×2 matrix of measured EE positions

% 4) Run the calibration
%[est_params, delta_params] = calibration(Q, P, fkine_sym, q_sym, param_sym, nominal_params);

% 5) Display results
%fprintf('Estimated l1 = %.4f m,  l2 = %.4f m\n', est_params(1), est_params(2));
%fprintf('Corrections   Δl1 = %+0.4f m, Δl2 = %+0.4f m\n', delta_params(1), delta_params(2));


% Solve for parameter corrections (least squares)
%delta_params = pinv(Phi) * DeltaP;
% Updated estimates
%est_params   = nominal_params(:) + delta_params;
%end
