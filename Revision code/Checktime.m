clc
clear
A = [-5,1;3,-5];
B = [0.5;1];
C = eye(2);
W1 = load('W1final.csv');
W2 = load('W2final.csv');
W3 = load('W3final.csv');
% W1 = load('W1.csv');
% W2 = load('W2.csv');
% W3 = load('W3.csv');
% W4 = load('W4.csv');
Ws = abs(W3)*abs(W2)*abs(W1);
% Ws = abs(W4)*abs(W3)*abs(W2)*abs(W1)

normws = norm(W3)*norm(W2)*norm(W1);
% normws = norm(W4)*norm(W3)*norm(W2)*norm(W1);
for q=1:1000
tic
n = size(A,1);
Metzler = A-B*Ws*C;
Hurwitz = eigs(A+B*Ws*C);

for i = 1:n
    if Hurwitz(i) > 0
       disp("No")
       break
    end
    for j = 1:n
        if j ~= i
           if Metzler(i,j) < 0
               disp("Non")
               break
           end
        end
    end
end
a=toc;
end
mean(a)