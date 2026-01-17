%% Ustawienia środowiska
clear; clc; close all;

% =========================================================
% STYL WYKRESÓW
% =========================================================
set(0, 'DefaultLineLineWidth', 2);
set(0, 'DefaultAxesFontSize', 11);
set(0, 'DefaultAxesXGrid', 'on'); 
set(0, 'DefaultAxesYGrid', 'on');
set(0, 'DefaultAxesGridLineStyle', ':');

%% 0. MODEL OBIEKTU
% G(s) = (-s + 6) / ((s+1)(s+7)) = (-s+6)/(s^2+8s+7)

A = [0.0, -7.0;
     1.0, -8.0];

B = [6.0;
    -1.0];

C = [0.0, 1.0]; % y = x2

% ZMIANA CZASU SYMULACJI NA 10 SEKUND
t_span = [0, 10];
t_eval = linspace(0, 10, 1000); 
x0 = [0; 0];

%% 1. ODPOWIEDŹ SKOKOWA – OTWARTA PĘTLA
u_open = 1.0;
open_loop_ode = @(t, x) A*x + B*u_open;

[t_open, y_open] = ode45(open_loop_ode, t_eval, x0);

figure;
plot(t_open, y_open(:, 2)); 
xlabel('Czas [s]'); ylabel('y(t)');
title('Odpowiedź skokowa – otwarta pętla');

%% 2. REGULATOR P – PORÓWNANIE
kp_values = [1, 3, 6];
styles = {'-', '--', ':'}; 

figure; hold on;
for i = 1:length(kp_values)
    kp = kp_values(i);
    style = styles{i};
    
    p_control_ode = @(t, x) p_dynamics(t, x, A, B, C, kp);
    [t_p, y_p] = ode45(p_control_ode, t_eval, x0);
    
    plot(t_p, y_p(:, 2), style, 'DisplayName', sprintf('Kp = %d', kp));
end

yline(1.0, 'k--', 'LineWidth', 1, 'DisplayName', 'Wartość zadana');
xlabel('Czas [s]'); ylabel('y(t)');
title('Regulator P – porównanie');
legend('Location', 'best'); hold off;

%% 3. LQR Z REFERENCJĄ
big_mat = [A, B; 
           C, 0];
rhs = [0; 0; 1];
res = big_mat \ rhs; 

x_ref = res(1:2);
u_ref = res(3);

Q = diag([10, 10]);
R = 1;

% Funkcja care z Control System Toolbox
[P, ~, ~] = care(A, B, Q, R); 
K = inv(R) * B' * P;

lqr_ode = @(t, x) A*x + B * (-K*(x - x_ref) + u_ref);
[t_lqr, y_lqr] = ode45(lqr_ode, t_eval, x0);

figure;
plot(t_lqr, y_lqr(:, 2), 'DisplayName', 'LQR');
hold on;
yline(1.0, 'r--', 'DisplayName', 'Wartość zadana');
xlabel('Czas [s]'); ylabel('y(t)');
title('Sterowanie LQR');
legend; hold off;

%% 4. LQR + OBSERWATOR LUENBERGERA
eig_cl = eig(A - B*K);
observer_poles = 3 * eig_cl;

% Funkcja place z Control System Toolbox
L = place(A', C', observer_poles)'; 

z0 = [0; 0; 0.2; 0.2]; 
obs_ode = @(t, z) observer_dynamics(t, z, A, B, C, K, L, x_ref, u_ref);

[t_obs, y_obs] = ode45(obs_ode, t_eval, z0);

figure;
plot(t_obs, y_obs(:, 2), 'DisplayName', 'x2 rzeczywiste');
hold on;
plot(t_obs, y_obs(:, 4), '--', 'DisplayName', 'x2 estymowane');
xlabel('Czas [s]'); ylabel('x2');
title('LQR z obserwatorem stanu');
legend; hold off;

disp('Gotowe – symulacja dla t=10s.');

%% --- FUNKCJE POMOCNICZE (na końcu pliku) ---

function dx = p_dynamics(~, x, A, B, C, kp)
    y = C * x;
    u = kp * (1.0 - y);
    dx = A * x + B * u;
end

function dz = observer_dynamics(~, z, A, B, C, K, L, x_ref, u_ref)
    x = z(1:2);
    x_hat = z(3:4);
    
    y = C * x;
    y_hat = C * x_hat;
    
    u = -K * (x_hat - x_ref) + u_ref;
    
    dx = A * x + B * u;
    dx_hat = A * x_hat + B * u + L * (y - y_hat);
    
    dz = [dx; dx_hat];
end