
%% Part 1: get solution path for SDEs

% Euler-Maruyama method on SDE

% dX_i/dt = s/tau - 1/tau X_i + c/tau Sum_{j!=i} r_j + 1/tau epsilon
% r_j = tanh(X_j)

randn('state', 100);
% parameters
tau = 10; M = 500; c = 1.1/(M-1); s = 0.03; Tau = 0.16;
% initiate X
Xzeros = zeros(1,M);
% time scale, total running time 1000ms
T = 1000; N = 1000; dt = T/N; 
% solution repeat fequency, the more repeat, the more stable of solution.
rep = 1;
% initiate random variables dt*N(0, 0.16^2)
dW = dt*0.16*randn(M,N,rep);


Xem = zeros(M,N,rep);
for repeat=1:rep
    Xtemp = zeros(1,N);
    Xem(:,1,repeat) = Xzeros';
    for time=1:N-1
        for index=1:M
            Xtemp(time+1) = Xem(index,time,repeat) + (s/tau - 1/tau*Xem(index,time,repeat))+ ...
            c/tau*(M*mean(tanh(Xem(:,time,repeat)))-tanh(Xem(index,time,repeat)))*dt + ...
            1/tau*dW(index,time,repeat);
            Xem(index,time+1,repeat) = Xtemp(time+1);
        end
    end
end
% average over all repeats, now its just one repeat.
X = mean(Xem, 3);

% (optional) plot solution path for X_1.
% plot([1:dt:T], X(1,:), 'r--'), hold on

%% Part 2: compute Markovian path fisher information estimator

% interested in parameter 's'
% initiate FIM
FIM = 0;
p = prob(X,tau,s,c,Tau,dt);
for j=1:N-1
    sum = 0;
    for i=1:M
        sum = sum + tau^2/((dt)^2*Tau^2)...
            * (X(i,j+1)- (X(i,j) + (s/tau-1/tau*X(i,j)...
            +c/tau*(M*mean(tanh(X(:,j)))-tanh(X(i,j))) )*dt)  )*dt/tau;
    end
    %FIM = FIM + prod(p(:,j))*sum^2;
    FIM = FIM + sum^2;
end
FIM = FIM/(N-1);
