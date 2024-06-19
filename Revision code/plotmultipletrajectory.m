clc
clear
close all
A=[-5 1;3 -5];
B=[.5;1];
C=eye(2);

W1 = load('W1.csv');
W2 = load('W2.csv');
W3 = load('W3.csv');
W4 = load('W4.csv');

%% plots
% Define the hyperplane parameters
%sector bound for tanh
sec_bound = 1;
coefficients_z1 = -(sec_bound^4)*abs(W4)*abs(W3)*abs(W2)*abs(W1)*C;
coefficients_z2 = (sec_bound^4)*abs(W4)*abs(W3)*abs(W2)*abs(W1)*C;
% Create a meshgrid of x and y values
[X, Y] = meshgrid(0:0.1:50, 0:0.1:50);

% Calculate the corresponding z values for both hyperplanes
z1 = coefficients_z1 * [X(:)'; Y(:)'];
z2 = coefficients_z2 * [X(:)'; Y(:)'];


% Reshape the z values to match the size of x and y
z1 = reshape(z1, size(X));
z2 = reshape(z2, size(X));


% Plot the hyperplanes
figure(1)
surf(X, Y, z1, 'FaceAlpha', 0.5, 'FaceColor', 'b', 'EdgeColor', 'none');
hold on;
surf(X, Y, z2, 'FaceAlpha', 0.5, 'FaceColor', 'r', 'EdgeColor', 'none');
hold on;

for h = 1:50
global output_history
output_history = zeros(1+size(B,2)+size(A,1),1);

% Set the initial condition
% disp('The sample initial condition')
x0 = 100*abs(rand(size(A,1),1));

%% Putting NN in the loop

% Define the time span for simulation
tspan = [0 3];

% Define the system dynamics as a function
f = @(t, x) A * x + B * NNcontrol(W1,W2,W3,W4,x,t);

% Solve the system using ode45
[t, x] = ode45(f, tspan, x0);

%% Plots

% 
%Checking if the input lies in the sector bound, by plotting the input trajectory as a function of x_1 and x_2
figure(1)
plot3(output_history(3,2:end),output_history(4,2:end),output_history(2,2:end),'k','LineWidth', 2)
set(gca, 'FontSize', 15);
legend('\Gamma_1x','\Gamma_2x','\pi(x)','Location', 'best','FontSize',20)
xlabel('x_1','FontSize',20);
ylabel('x_2','FontSize',20);
zlabel('Output','FontSize',20);

% 
% %Plot the States
% figure(2)
% hold on
% plot(t, x(:,1), 'LineWidth', 2);
% set(gca, 'FontSize', 25);
% ylabel('x_1','FontSize',30);
% xlabel('Time','FontSize',30);
% % title('Simulation of LTI System with NN controller','FontSize',15);
% % legend('x_1','Location', 'best','FontSize',20);
% grid on;
% 
% %Plot the States
% figure(3)
% hold on
% plot(t, x(:,2), 'LineWidth', 2);
% set(gca, 'FontSize', 25);
% ylabel('x_2','FontSize',30);
% xlabel('Time','FontSize',30);
% % title('Simulation of LTI System with NN controller','FontSize',15);
% % legend('x_2','Location', 'best','FontSize',20);
% grid on;

%Plot the States
figure(4)
hold on
plot(t, x(:,1),'b',t,x(:,2),'r', 'LineWidth', 2);
ylabel('State Variables','FontSize',20);
xlabel('Time','FontSize',20);
set(gca, 'FontSize', 15);
% title('Simulation of Continuous-Time System with NN controller','FontSize',15);
legend('x_1','x_2','Location', 'best','FontSize',20);
grid on;
end

%% define the NN controller output
function control_input = NNcontrol(W1,W2,W3,W4,xminus,t)
global output_history
v1 = W1*xminus;
w1 = tanh(v1);
v2 = W2*w1;
w2 = tanh(v2);
v3 = W3*w2;
w3 = tanh(v3);
v4 = W4*w3;
control_input = v4;

output_data = [t;v4;xminus];
output_history = [output_history output_data];
end
