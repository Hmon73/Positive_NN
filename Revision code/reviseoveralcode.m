clc
clear
close all
A=[-5 1;3 -5];
B=[.5;1];
C=eye(2);
syms s11 s12 s21 s22 p
Sigma1 = [s11 s12];
Sigma2 = [s21 s22];


%% Check condition 1
disp('The following matrix should be Metzler')
cond1 = A+B*[s11 s12]*C

%% Check condition 2
disp('The following matrix should be Hurwitz')
cond2 = A+B*[s21 s22]*C
disp('The eigenvalues are:')
eig(cond2)
disp('The Characteristic Polynomial is:')
charpoly(cond2,p)

gh = -2:.1:5;
jh = 10-gh/2;
figure
plot(gh,jh)
hold on
jh = (22-7*gh/2)*2/13;
plot(gh,jh)
hold off

%% lqr gain
% Q = 5*eye(2);
Q = diag([0.3, 3])
R = .3;
[K_lqr,S,P] = lqr(A,B,Q,R);
K_lqr = -K_lqr;
disp('The lqr gain (-K) is:');
disp(K_lqr)


K = K_lqr

% K =[3, 1.5]

disp('This matrix should be Metzler (This is for selected gain)');
ans = A+B*K
disp('The eigenvalues should be negative (This is for selected gain)');
eig(ans)
% 
% %% Generating Data to train NN
% 
% % Generating the input sample to NN
% 
% % Generating sample vectors between 0 and 0.01
% num_vectors = 500;
% vector4 = rand(num_vectors, 2)*.01;
% 
% % Generating sample vectors between 0 and 0.1
% num_vectors = 500;
% vector3 = rand(num_vectors, 2)*0.1;
% 
% % Generating sample vectors between 0 and 1
% num_vectors = 500;
% vector2 = rand(num_vectors, 2)*1;
% 
% % Generating sample vectors between 0 and 10
% num_vectors = 500;
% vector1 = rand(num_vectors, 2)*10;
% 
% % Gathering all the random samples. it's a 2000 by 2 vector
% vector = [vector1;vector2;vector3;vector4];
% 
% % Shuffling the rows of the sample data to get a random sample between 0 and 10
% originalMatrix = vector;
% 
% % Generate a random permutation of row indices
% randomIndices = randperm(size(originalMatrix, 1));
% 
% % Use the permutation to shuffle the rows of the matrix
% shuffledMatrix = originalMatrix(randomIndices, :);
% 
% 
% % Generating the Target data by u = K * sample input
% u = (K * shuffledMatrix')';
% 
% 
% % Save the data vector and matrix to text files
% dlmwrite('K_xs.txt', shuffledMatrix);
% dlmwrite('K_us.txt', u);


%% reading the weights
W1 = load('W1.csv');
W2 = load('W2.csv');
W3 = load('W3.csv');
W4 = load('W4.csv');

%% Putting NN in the loop

global output_history
output_history = zeros(1+size(B,2)+size(A,1),1);

%Simulation of the system
% c = max (|a|,|b|)
sec_bound = 1;

% Define the hyperplane parameters
% coefficients_z1 = [-1, -2];
% coefficients_z2 = [2 3.1];
% coefficients_z3 = K;
coefficients_z4 = (sec_bound^4)*abs(W4)*abs(W3)*abs(W2)*abs(W1)*C;
coefficients_z5 = -(sec_bound^4)*abs(W4)*abs(W3)*abs(W2)*abs(W1)*C;

%random initial condition for system
x0 = 100*abs(rand(size(A,1),1))
% Create a meshgrid of x and y values
[x, y] = meshgrid(-abs(x0(1)):0.1:abs(x0(1)), -abs(x0(2)):0.1:abs(x0(2)));

% Calculate the corresponding z values for both hyperplanes
% z1 = coefficients_z1 * [x(:)'; y(:)'];
% z2 = coefficients_z2 * [x(:)'; y(:)'];
% z3 = coefficients_z3 * [x(:)'; y(:)'];
z4 = coefficients_z4 * [x(:)'; y(:)'];
z5 = coefficients_z5 * [x(:)'; y(:)'];

% Reshape the z values to match the size of x and y
% z1 = reshape(z1, size(x));
% z2 = reshape(z2, size(x));
% z3 = reshape(z3, size(x));
z4 = reshape(z4, size(x));
z5 = reshape(z5, size(x));

% Plot the hyperplanes
figure;
% surf(x, y, z1, 'FaceAlpha', 0.5, 'FaceColor', 'b', 'EdgeColor', 'none');
% hold on;
% surf(x, y, z2, 'FaceAlpha', 0.5, 'FaceColor', 'r', 'EdgeColor', 'none');
% hold on;
surf(x, y, z4, 'FaceAlpha', 0.5, 'FaceColor', 'r', 'EdgeColor', 'none');
hold on;
surf(x, y, z5, 'FaceAlpha', 0.5, 'FaceColor', 'b', 'EdgeColor', 'none');
hold on;

xlabel('X_1');
ylabel('X_2');
zlabel('output');
% title('Hyperplane Plot');
grid on;
hold on
% legend('z1 = [-1, -2] \cdot [x; y]', 'z2 = [1, 4] \cdot [x; y]');



% Define the time span for simulation
tspan = [0 5];

% Define the system dynamics as a function
f = @(t, x) A * x + B * NNcontrol(W1,W2,W3,W4,x,t);

% Solve the system using ode45
[t, x] = ode45(f, tspan, x0);


%Checking if the input lies in the sector bound, by plotting the input trajectory as a function of x_1 and x_2
plot3(output_history(3,2:end),output_history(4,2:end),output_history(2,2:end),'k','LineWidth', 2)
legend('\Sigma_1','\Sigma_2','Input trajectory','Location', 'best')

%Plot the States
figure
hold on
plot(t, x, 'LineWidth', 2);
ylabel('State Variables','FontSize',20);
xlabel('Time','FontSize',20);
title('Simulation of Continuous-Time System with NN controller','FontSize',15);
legend('x_1','x_2','Location', 'best','FontSize',20);
grid on;

%% %define the NN controller output
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


