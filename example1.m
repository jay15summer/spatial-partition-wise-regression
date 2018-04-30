clear; close all; 
% 1. data simulation
[x, y] = meshgrid(1:0.1:10,1:0.1:10);
x = reshape(x, size(x,1)*size(x,2),1);
y = reshape(y, size(y,1)*size(y,2),1);
S = [x y];% locations
X1 = normrnd(2, 0.5, length(S), 1);% independent variable 1
X2 = normrnd(1, 0.1, length(S), 1);% independent variable 2
X = [X1 X2]; 
Y = zeros(length(S), 1); % initialize dependent variable
% simulate partition-wise regression with varying coefficients (varying
% sharply with crisp boundary)
for i = 1: length(S)
    if S(i, 1) < 3 && S(i, 2) < 5
        Y(i) = 1 + 2*X1(i) + X2(i) + normrnd(0, 0.2);
    elseif S(i, 1) < 6 && S(i, 2) < 8 
        Y(i) = 1.5 + 3*X1(i) + 2*X2(i) + normrnd(0, 0.2);
    else
        Y(i) = 2 + 3.5*X(i) + 1.5*X2(i) + normrnd(0, 0.2);
    end
end
% 2. main function for spatial partitioning
h = 9; v = 9; T = 0.1; 
[partition_all, partition_slt] = spatial_partition_reg(S, X, Y, h, v, T);
