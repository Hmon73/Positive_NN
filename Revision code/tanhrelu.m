clc
clear
close all

% ReLu
x = -3:0.1:3;
x1 = -3:.1:0;
x2 = 0:.1:3;
y1 = tanh(x);
y2 = x;
y3 = zeros(size(x));
y4 = zeros(size(x1));
y5 = x2;

figure;
% plot(x, y1, 'LineWidth', 3);  % Plot y1 (tanh) with thicker line
plot(x1,y4+.02, 'b' ,x2,y5,'b', 'LineWidth', 3)
% hold on
% plot(x2,y5,'b', 'LineWidth', 3)
hold on;  % Hold the current plot


% Display x and y axes as arrows
quiver(-3, 0, 6.5, 0, 'MaxHeadSize', 0.1,'LineWidth', 2.5, 'Color', 'k'); % x-axis arrow
quiver(0, -3, 0, 6.5, 'MaxHeadSize', 0.1, 'LineWidth', 2.5, 'Color', 'k'); % y-axis arrow
plot(x, y2+.05, 'r--', 'LineWidth', 3);  % Plot y2 (linear) with thicker line
plot(x, y3+.07, 'r--', 'LineWidth', 3);  % Plot y3 (horizontal line) with thicker line
% Label x and y axes
xlabel('x', 'FontSize', 120); % Larger font size for x label
ylabel('f(x)', 'FontSize', 120); % Larger font size for y label

% Set the axis limits to fit all data
xlim([-3 3]);
ylim([-3 3]);

% Increase font size of tick labels on both axes
set(gca, 'FontSize', 30); 

% Add grid
grid on;
title('ReLu function', 'FontSize', 30);
% Add legend
% legend('tanh(x)', 'y = x', 'y = 0');

% Release hold
hold off;

% %% Tanh

figure;
plot(x, y1, 'LineWidth', 3);  % Plot y1 (tanh) with thicker line
hold on;  % Hold the current plot


% Display x and y axes as arrows
quiver(-3, 0, 6.5, 0, 'MaxHeadSize', 0.1,'LineWidth', 2.5, 'Color', 'k'); % x-axis arrow
quiver(0, -3, 0, 6.5, 'MaxHeadSize', 0.1, 'LineWidth', 2.5, 'Color', 'k'); % y-axis arrow
plot(x, y2, 'r--', 'LineWidth', 3);  % Plot y2 (linear) with thicker line
plot(x, y3, 'r--', 'LineWidth', 3);  % Plot y3 (horizontal line) with thicker line
% Label x and y axes
xlabel('x', 'FontSize', 120); % Larger font size for x label
ylabel('f(x)', 'FontSize', 120); % Larger font size for y label

% Set the axis limits to fit all data
xlim([-3 3]);
ylim([-3 3]);

% Increase font size of tick labels on both axes
set(gca, 'FontSize', 30); 

% Add grid
grid on;
title('tanh function', 'FontSize', 30);
% Add legend
% legend('tanh(x)', 'y = x', 'y = 0');

% Release hold
hold off;