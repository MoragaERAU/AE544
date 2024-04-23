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

%% Functions
function [modes, frequencies, eigenvalues] = mode_shapes(M, K)
[eigenvectors, eigenvalues] = eig(K, M);
frequencies = sqrt(diag(eigenvalues));
modes = eigenvectors;

% Normalize mode shapes
modes = modes ./ max(abs(modes));
end
