clc
clear
close all

tspan = (0:0.01:5);
reference = sin(2*pi*(1/5)*tspan);

% Define mass matrix M and stiffness matrix K
m = [1, 1, 1, 1, 1]; % Masses
L = 10;
n = 5;
EA = 1;
k = n*EA/L;

K = [3*k -k 0 0 0;...
    -k 2*k -k 0 0;...
    0 -k 2*k -k 0;...
    0 0 -k 2*k -k;...
    0 0 0 -k 3*k;];

% Mass matrix
M = diag(m);

% Calculate mode shapes and frequencies
[modes, frequencies, eigenvalues] = mode_shapes(M, K);
modes = modes*1;

figure
plot(tspan,reference,'r')
hold on
for i = 1:5
    subplot(5,1,i)
    if i ==1
        plot(modes(:,i),'o-')
%     elseif i ==5
%         plot(modes(:,i),'o-')
    else
        plot(modes(:,1),'o-')


    end
    title(['Mode Shape ', num2str(i)])
    xlabel('Mass Index')
    ylabel('Displacement')
    grid on
    grid minor
end



% modes = [zeros(1,5);modes;zeros(1,5)]
% Plot mode shapes
figure;
for i = 1:5
    subplot(5, 1, i)
    plot(modes(:, i), 'o-')
    hold on
    plot(tspan,reference,'r')
    title(['Mode Shape ', num2str(i)])
    xlabel('Mass Index')
    ylabel('Displacement')
end


%% Question 5
clc;clear;close all
%phi = (x/L)^i+1

%Linear Case
%Constants
EI = 1.4*10^4;
rho = 1.2;
L = 10;
n = 3;
omega = 1;

%mass/stiffness matrix
m = zeros(3,3);
k = zeros(3,3);
for i = 1:n
    for j = 1:n
        m(i,j) = rho*L/(i+j+3);
        k(i,j) = (EI/(rho*L^4))*i*j*(i+1)*(j+1)/(i+j-1);
    end
end

tspan = (0:0.001:10);
cond0 = [0.1; 0; 0; 0; 0; 1];  % Ensure initial condition is a column vector

options = odeset('RelTol',1e-12,'AbsTol',1e-12);
[~, q] = ode45(@(time, q) linbeam(time, q, m, k, omega), tspan, cond0, options);

subplot(3,1,1)
plot(tspan, q(:,1))
grid on
grid minor
xlabel('Timespan (s)')
ylabel('\delta y')
title('Deflection of Beam')

subplot(3,1,2)
plot(tspan, q(:,2))
grid on
grid minor
xlabel('Timespan (s)')
ylabel('\delta u')
title(' Axial deflection of Beam')

subplot(3,1,3)
plot(tspan, q(:,3))
grid on
grid minor
xlabel('Timespan (s)')
ylabel('\theta')
title('Angle of rotation')

%% Functions
function [modes, frequencies, eigenvalues] = mode_shapes(M, K)
[eigenvectors, eigenvalues] = eig(K, M);
frequencies = sqrt(diag(eigenvalues));
modes = eigenvectors;

% Normalize mode shapes
modes = modes ./ max(abs(modes));
end

function solutions = linbeam(~, q, m, k, omega)
y = q(1);
u = q(2);
theta = q(3);

dy = q(4);
du = q(5);
dtheta = q(6);

dq = [dy; du; dtheta];
q = [y; u; theta];

%     ddq = inv(m) * ((omega^2 * m - k) * q);
ddq =     (((omega^2 * m - k)/m) * q);


solutions = [dq; ddq];
end
