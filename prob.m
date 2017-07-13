function [p] = prob(X,tau,s,c,Tau,dt)
% probability function for 

% dX_i/dt = s/tau - 1/tau X_i + c/tau Sum_{j!=i} r_j + 1/tau epsilon
% r_j = tanh(X_j)

% input: solution path for all X_i, i=1,...,M, 
% ouput: prob(X_i(t+1)|X_i(t)) for all i=1,...,M and t=1,...,T-1

    [index, time] = size(X);
    p = zeros(index, time-1);
    for i=1:index
        for j=1:time-1
            p(i,j+1) = tau/(sqrt(2*pi)*dt*Tau) ...
                *exp( -tau^2/(2*(dt)^2*Tau^2)...
                *( X(i,j+1) -...
                ( X(i,j) + (s/tau-1/tau*X(i,j)+c/tau*(sum(tanh(X(:,j)))-tanh(X(i,j))) )*dt ) )^2);
        end
    end
end