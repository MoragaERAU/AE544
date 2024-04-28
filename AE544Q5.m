clc
clear
close all

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
        k(i,j) =   
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
function solutions = linbeam(~, q, m, k, omega)
y = q(1);
u = q(2);
theta = q(3);

dy = q(4);
du = q(5);
dtheta = q(6);

dq = [dy; du; dtheta];
q = [y; u; theta];

%   ddq = inv(m) * ((omega^2 * m - k) * q);
 ddq =     (((omega^2 * m - k)/m) * q);


solutions = [dq; ddq];
end
