% This code plots the sector bounds and the outputs of the NN for random
% inputs.
%% Uploading the weights
clc
clear
close all

W1 = load('W1.csv');
W2 = load('W2.csv');
W3 = load('W3.csv');
W4 = load('W4.csv');
%% generating random positive vectors z and calculating upper and lower bounds
m = 100;
for o = 1:m
    z = abs(rand(size(W1,2),1));
    NNoutput(o) = NN4(W1,W2,W3,W4,z);
    Sigma1(o) = [-3 -2]*z;
    Sigma2(o) = [3,1.76]*z;
    UpperSectorBound(o) = abs(W4)*abs(W3)*abs(W2)*abs(W1)*z;
    LowerSectorBound(o) = -abs(W4)*abs(W3)*abs(W2)*abs(W1)*z;
end

%% plotting
figure
hold on

plot(1:m, Sigma1,'LineWidth',2, 'color',[0 0.4470 0.7410],'DisplayName','\Sigma_1z');
plot(1:m, Sigma2,'LineWidth',2,'Color',[0.6350 0.0780 0.1840],'DisplayName','\Sigma_2z');
fill_between_lines(Sigma1, Sigma2, [0.8 0.8 0.8]); % Fills the area between Sigma1 and Sigma2 with gray
set(gca, 'FontSize', 15);
plot(1:m, LowerSectorBound,'b','LineWidth',2,'DisplayName','\Gamma_1');
plot(1:m, UpperSectorBound,'r','LineWidth',2,'DisplayName','\Gamma_2');
plot(1:m, NNoutput,'k','LineWidth',2,'DisplayName','\pi(z)');
grid on
legend('\Sigma_1z','\Sigma_2z','','\Gamma_1z', '\Gamma_2z','\pi(z)','FontSize',20,'Orientation','horizontal')
xlabel('Random Positive Input','FontSize',20)
ylabel('Output','FontSize',20)


%%
% This code is for a 3 layer NN

%% Uploading the weights
% clc
% clear
% close all

W1 = load('W1m.csv');
W2 = load('W2m.csv');
W3 = load('W3m.csv');

%% generating random positive vectors z and calculating upper and lower bounds
m = 100;
for o = 1:m
    z = abs(rand(size(W1,2),1));
    NNoutput(o) = NN3(W1,W2,W3,z);
    Sigma1(o) = [-3 -2]*z;
    Sigma2(o) = [3,1.76]*z;
    UpperSectorBound(o) = abs(W3)*abs(W2)*abs(W1)*z;
    LowerSectorBound(o) = -abs(W3)*abs(W2)*abs(W1)*z;
end

%% plotting
figure
hold on

plot(1:m, Sigma1,'LineWidth',2, 'color',[0 0.4470 0.7410],'DisplayName','\Sigma_1z');
plot(1:m, Sigma2,'LineWidth',2,'Color',[0.6350 0.0780 0.1840],'DisplayName','\Sigma_2z');
fill_between_lines(Sigma1, Sigma2, [0.8 0.8 0.8]); % Fills the area between Sigma1 and Sigma2 with gray
set(gca, 'FontSize', 15);
plot(1:m, LowerSectorBound,'b','LineWidth',2,'DisplayName','\Gamma_1');
plot(1:m, UpperSectorBound,'r','LineWidth',2,'DisplayName','\Gamma_2');
plot(1:m, NNoutput,'k','LineWidth',2,'DisplayName','\pi(z)');
grid on
legend('\Sigma_1z','\Sigma_2z','','\Gamma_1z', '\Gamma_2z','\pi(z)','FontSize',20,'Orientation','horizontal')
xlabel('Random Positive Input','FontSize',20)
ylabel('Output','FontSize',20)


%% fill between lines function
function fill_between_lines(line1, line2, color)
    x = 1:length(line1);
    fill([x, fliplr(x)], [line1, fliplr(line2)], color,'FaceAlpha',0.5);
end
%% function to calculate NN output
function NNoutput = NN4(W1,W2,W3,W4,z)
v_1 = W1*z;
w_1 = tanh(v_1);
v_2 = W2*w_1;
w_2 = tanh(v_2);
v_3 = W3*w_2;
w_3 = tanh(v_3);
v_4 = W4*w_3;
NNoutput=v_4;
end

%% function to calculate NN output
function NNoutput = NN3(W1,W2,W3,z)
v_1 = W1*z;
w_1 = tanh(v_1);
v_2 = W2*w_1;
w_2 = tanh(v_2);
v_3 = W3*w_2;
NNoutput=v_3;
end


% %% 2ta W
% 
% %% Uploading the weights
% clc
% clear
% close all
% 
% W1 = load('W1.csv');
% W2 = load('W2.csv');
% W3 = load('W3.csv');
% W4 = load('W4.csv');
% 
% % W1 = randn(size(W1));
% % W2 = randn(size(W2));
% % W3 = randn(size(W3));
% 
% %% generating random positive vectors z and calculating upper and lower bounds
% m = 100;
% for o = 1:m
%     z = abs(rand(size(W1,2),1));
%     NNoutput(o) = NN(W1,W2,z);
%     Sigma1(o) = [-3 -2]*z;
%     Sigma2(o) = [3,1.76]*z;
% %     UpperSectorBound(o) = abs(W2)*abs(W1)*z;
% %     LowerSectorBound(o) = -abs(W2)*abs(W1)*z;
% % end
% % 
% % %% plotting
% % figure
% % hold on
% % 
% % plot(1:m, Sigma1,'LineWidth',2, 'color',[0 0.4470 0.7410],'DisplayName','\Sigma_1z');
% % plot(1:m, Sigma2,'LineWidth',2,'Color',[0.6350 0.0780 0.1840],'DisplayName','\Sigma_2z');
% % fill_between_lines(Sigma1, Sigma2, [0.8 0.8 0.8]); % Fills the area between Sigma1 and Sigma2 with gray
% % set(gca, 'FontSize', 15);
% % plot(1:m, LowerSectorBound,'b','LineWidth',2,'DisplayName','\Gamma_1');
% % plot(1:m, UpperSectorBound,'r','LineWidth',2,'DisplayName','\Gamma_2');
% % plot(1:m, NNoutput,'k','LineWidth',2,'DisplayName','\pi(z)');
% % grid on
% % legend('\Gamma_1z','\Gamma_2z','','\Sigma_1z', '\Sigma_2z','\pi(z)','FontSize',20,'Orientation','horizontal')
% % xlabel('Random Positive Input','FontSize',20)
% % ylabel('Output','FontSize',20)
% % figure
% % hold on
% % plot(1:m, LowerSectorBound,'b','LineWidth',2);
% % plot(1:m, UpperSectorBound,'r','LineWidth',2);
% % plot(1:m, Sigma1,'LineWidth',2, 'color',[0 0.4470 0.7410]);
% % plot(1:m, Sigma2,'LineWidth',2,'Color',[0.6350 0.0780 0.1840]);
% % plot(1:m, NNoutput,'k','LineWidth',2);
% % set(gca, 'FontSize', 15);
% % grid on
% % legend('\Gamma_1z','\Gamma_2z','\pi(z)','FontSize',20)
% % xlabel('Random Positive Input','FontSize',20)
% % ylabel('Output','FontSize',20)
% % legend('Location','best')
% 
% % figure
% % data = [LowerSectorBound; UpperSectorBound; NNoutput]';
% % % Plotting
% % bar(data, 'stacked');
% 
% function fill_between_lines(line1, line2, color)
%     x = 1:length(line1);
%     fill([x, fliplr(x)], [line1, fliplr(line2)], color,'FaceAlpha',0.5);
% end
% 
% %% function to calculate NN output
% function NNoutput = NN(W1,W2,z)
% v_1 = W1*z;
% w_1 = tanh(v_1);
% v_2 = W2*w_1;
% w_2 = tanh(v_2);
% v_3 = W3*w_2;
% w_3 = tanh(v_3);
% v_4 = W4*w_3;
% NNoutput=v_2;
% end
% 
% 
% 
% 
% 
%  % Uploading the weights
% clc
% clear
% close all
% 
% W1 = load('W1.csv');
% W2 = load('W2.csv');
% W3 = load('W3.csv');
% W4 = load('W4.csv');
% 
% 
% %% generating random positive vectors z and calculating upper and lower bounds
% m = 100;
% % % Define z1 and z2
% % a1 = 0:0.1:2;
% % a2 = 0:0.1:2;
% % 
% % % Create a meshgrid of z1 and z2
% % [A1, A2] = meshgrid(a1, a2);
% % 
% % % Stack Z1 and Z2 into a single matrix z
% % A = [A1(:), A2(:)]';
% % m=size(A,2);
% for o = 1:m
%     z = abs(rand(size(W1,2),1));
%     % z = A(:,o)
%     NNoutput(o) = NN(W1,W2,W3,z);
%     Sigma1(o) = [-3 -2]*z;
%     Sigma2(o) = [3,1.76]*z;
%     UpperSectorBound(o) = abs(W3)*abs(W2)*abs(W1)*z;
%     LowerSectorBound(o) = -abs(W3)*abs(W2)*abs(W1)*z;
% end
% 
% %% plotting
% figure
% hold on
% 
% % plot(1:m, Sigma1,'LineWidth',2, 'color',[0 0.4470 0.7410],'DisplayName','\Sigma_1z');
% % plot(1:m, Sigma2,'LineWidth',2,'Color',[0.6350 0.0780 0.1840],'DisplayName','\Sigma_2z');
% % fill_between_lines(Sigma1, Sigma2, [0.8 0.8 0.8]); % Fills the area between Sigma1 and Sigma2 with gray
% % set(gca, 'FontSize', 15);
% plot(1:m, LowerSectorBound,'b','LineWidth',2,'DisplayName','\Gamma_1');
% plot(1:m, UpperSectorBound,'r','LineWidth',2,'DisplayName','\Gamma_2');
% plot(1:m, NNoutput,'k','LineWidth',2,'DisplayName','\pi(z)');
% grid on
% legend('\Gamma_1z','\Gamma_2z','','\Sigma_1z', '\Sigma_2z','\pi(z)','FontSize',20,'Orientation','horizontal')
% xlabel('Random Positive Input','FontSize',20)
% ylabel('Output','FontSize',20)
% % figure
% % hold on
% % plot(1:m, LowerSectorBound,'b','LineWidth',2);
% % plot(1:m, UpperSectorBound,'r','LineWidth',2);
% % plot(1:m, Sigma1,'LineWidth',2, 'color',[0 0.4470 0.7410]);
% % plot(1:m, Sigma2,'LineWidth',2,'Color',[0.6350 0.0780 0.1840]);
% % plot(1:m, NNoutput,'k','LineWidth',2);
% % set(gca, 'FontSize', 15);
% % grid on
% % legend('\Gamma_1z','\Gamma_2z','\pi(z)','FontSize',20)
% % xlabel('Random Positive Input','FontSize',20)
% % ylabel('Output','FontSize',20)
% % legend('Location','best')
% 
% % figure
% % data = [LowerSectorBound; UpperSectorBound; NNoutput]';
% % % Plotting
% % bar(data, 'stacked');
% 
% function fill_between_lines(line1, line2, color)
%     x = 1:length(line1);
%     fill([x, fliplr(x)], [line1, fliplr(line2)], color,'FaceAlpha',0.5);
% end
% 
% 
% %% function to calculate NN output
% function NNoutput = NN(W1,W2,W3,W4,z)
% v_1 = W1*z;
% w_1 = tanh(v_1);
% v_2 = W2*w_1;
% w_2 = tanh(v_2);
% v_3 = W3*w_2;
% w_3 = tanh(v_3);
% v_4 = W4*w_3;
% NNoutput=v_4;
% end
