% LMPC for SRB attitude (QUATERNION-BASED) - Final Verified Version
clear; clc;

deg2rad = pi/180;

% Inertia (example)
% I = diag([0.04565442, 0.03779987, 0.01585311]);

%%%  SATELLITE INERTIA %%%
I = [1.10, 0.05, 0.00;
     0.05, 1.90, -0.01;
     0.00, -0.01, 1.12];

Iinv = inv(I);

% Simulation parameters
dt = 0.05;

% LMPC Parameters
N = 5; % Prediction Horizon

% --- TUNING FOR QUATERNION CONTROLLER ---
% Q penalizes: [q_vec_error(3), omega_error(3)]
Q_diag = [250, 180, 100, 10, 15, 5]; % Aggressive on attitude, moderate on rates
R_val  = 0.01;                  % Low penalty on control effort

% Input Constraints
tau_max = 10;
umax = tau_max*ones(3,1);
umin = -umax;

% --- INITIAL AND TARGET STATES (QUATERNION) ---
% State x: [q0, q1, q2, q3, wx, wy, wz]' = [w, x, y, z, rates]'
x = [1; 0; 0; 0; 0; 0; 0]; % Start at identity rotation, zero velocity

% Initial attitude: [12, 31, -40] degrees (ZYX)
eul_init = [12; 31; -40] * deg2rad;
q_init = eul2quat(eul_init', 'ZYX')';

% Initial body rates: [10, -20, 12] deg/s converted to rad/s
omega_init = [10; -20; 12] * deg2rad;

x = [q_init; omega_init];

% Target state: Convert 40;-30;90 degrees euler target to quaternion
eul_target = [40;-30;90]*deg2rad;
q_target = eul2quat(eul_target', 'ZYX')'; % MATLAB's eul2quat output is [w,x,y,z]
omega_target=[0; 0; 0];
x_target_quat = [q_target; omega_target];

tol = 1e-6;
max_steps = 700;

% Initialize linearization input
u0 = zeros(3,1);

% Storage for history
X_hist = x';
U_hist = [];
EUL_hist = quat2eul(x(1:4)','ZYX')/deg2rad;

quadprog_options = optimoptions('quadprog', 'Display', 'off');

for step = 1:max_steps
    % --- ERROR CALCULATION (QUATERNION) ---
    q_current = x(1:4)'; % as row vector
    q_target_inv = [q_target(1), -q_target(2:4)'];
    q_error = quatmultiply(q_target_inv, q_current);
    
    err = norm(q_error(2:4)); % Error is the vector part of the error quaternion
    if err < tol && norm(x(5:7)) < tol
        disp(['Converged at step ', num2str(step)]);
        break;
    end

    % ---- 1. Linearize and Discretize ----
    f = @(x,u) dynamics_quat(x,u,Iinv);
    A = numerical_jacobian_x(f,x,u0); 
    B = numerical_jacobian_u(f,x,u0);
    sysd = c2d(ss(A,B,[],[]), dt, 'zoh');
    Ad = sysd.A;
    Bd = sysd.B;
    
    % ---- 2. Formulate the LMPC Problem ----
    nx = size(Ad, 1); % Should be 7
    nu = size(Bd, 2); % Should be 3
    
    % Build Prediction Matrices: ? and ?
    Phi = zeros(N*nx, nx);
    Gamma = zeros(N*nx, N*nu);
    
    % Populate Phi and the first block-column of Gamma
    for i = 1:N
        Phi((i-1)*nx+1:i*nx, :) = Ad^i;
        Gamma((i-1)*nx+1:i*nx, 1:nu) = (Ad^(i-1)) * Bd;
    end
    
    % --- CORRECTED AND VERIFIED GAMMA CONSTRUCTION ---
    % Populate the rest of Gamma by copying and shifting the first column
    for i = 2:N
        col_start = (i-1)*nu + 1;
        col_end = i*nu;
        row_start_dest = (i-1)*nx + 1;
        num_rows_to_copy = (N-i+1)*nx;
        
        Gamma(row_start_dest:end, col_start:col_end) = Gamma(1:num_rows_to_copy, 1:nu);
    end
    % --- END CORRECTION ---

    % --- Build Cost Matrices based on error dynamics ---
    Q_full_state = diag([0, Q_diag(1:3), Q_diag(4:6)]); % 7x7 matrix, no penalty on q0
    R_mat = R_val*eye(nu);
    
    Q_bar = kron(eye(N), Q_full_state);
    R_bar = kron(eye(N), R_mat);

    H = 2 * (Gamma' * Q_bar * Gamma + R_bar);
    H=(H+H')/2;
    X_ref = repmat(x_target_quat, N, 1);
    error_term = Phi * x - X_ref;
    f = 2 * (error_term' * Q_bar * Gamma)';
    
    U_min = repmat(umin, N, 1);
    U_max = repmat(umax, N, 1);

    % ---- 3. Solve the Quadratic Program ----
    [U_opt, ~, exitflag] = quadprog(H, f, [], [], [], [], U_min, U_max, [], quadprog_options);
    
    if exitflag ~= 1
        warning('quadprog failed or found no feasible solution at step %d', step);
        u = zeros(3,1); % Failsafe
    else
        u = U_opt(1:nu);
    end
    
    u = min(max(u,umin),umax);

    % ---- 4. Integrate and Normalize ----
    xdot = dynamics_quat(x,u,Iinv);
    x = x + dt*xdot;
    
    x(1:4) = x(1:4) / norm(x(1:4)); % CRITICAL: Normalize quaternion

    u0 = u; % Update linearization input for next step

    % ---- Log Data ----
    X_hist = [X_hist; x'];
    U_hist = [U_hist; u'];
    EUL_hist = [EUL_hist; quat2eul(x(1:4)','ZYX')/deg2rad];

end

% ---- Plot Results ----
t = (0:size(X_hist,1)-1)*dt;
figure;
subplot(2,1,1)
plot(t, EUL_hist); grid on;
title('Attitude (Euler Angles from Quaternion State)');
ylabel('Angles (deg)');
legend('\psi (Z)','\theta (Y)','\phi (X)');

subplot(2,1,2)
plot(t, X_hist(:,5:7)); grid on;
title('Angular Velocity');
ylabel('Angular rates (rad/s)');
xlabel('Time (s)');
legend('w_x','w_y','w_z');

figure;
plot(t(1:end-1), U_hist); grid on;
title('Control Input (Torque)');
ylabel('Torque (N*m)');
xlabel('Time (s)');
legend('\tau_x','\tau_y','\tau_z');

% ---- NEW DYNAMICS & JACOBIAN FUNCTIONS ----
function xdot = dynamics_quat(x,u,Iinv)
    q = x(1:4); % [w, x, y, z]
    w = x(5:7); % [wx, wy, wz]
    
    % Quaternion kinematics
    Omega = [0, -w(1), -w(2), -w(3);
             w(1), 0, w(3), -w(2);
             w(2), -w(3), 0, w(1);
             w(3), w(2), -w(1), 0];
    
         q_dot = 0.5 * Omega * q;
    
    omega_dot = Iinv*u;
    
    xdot = [q_dot;
           omega_dot];
    
end

function A = numerical_jacobian_x(f,x,u)
    nx = length(x); fx = f(x,u);
    A = zeros(nx,nx);
    eps = 1e-6;
    for i=1:nx
        x2 = x; x2(i)=x2(i)+eps;
        A(:,i)=(f(x2,u)-fx)/eps;
    end
end

function B = numerical_jacobian_u(f,x,u)
   nu = length(u); fx = f(x,u);
   nx = length(fx);
    B = zeros(nx,nu);
  eps = 1e-6;
   for i=1:nu
        u2 = u; u2(i)=u2(i)+eps;
        B(:,i)=(f(x,u2)-fx)/eps;
   end
end



