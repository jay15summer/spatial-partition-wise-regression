clear; close all;
% 1. data simulation
[s1, s2] = meshgrid(1:0.1:10,1:0.1:10);
s1 = reshape(s1, size(s1,1)*size(s1,2),1);
s2 = reshape(s2, size(s2,1)*size(s2,2),1);
S = [s1 s2];% locations
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
