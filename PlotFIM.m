%% 2D plot of FIM for a range of c (Recurrent excitation) and Tau (Noise)

% Initialization
clear; close all; clc

% The range of parameters
c_bar = 1.0:0.2:1.4;
Tau = 0.05:0.1:0.25;

% Get FIM
FIM_matrix = zeros(size(c_bar, 2), size(Tau, 2));
for i = 1:size(c_bar, 2)
    for j = 1:size(Tau, 2)
        FIM_matrix(i, j) = FIMfunction(c_bar(i), Tau(j));
        j
    end
end

   
