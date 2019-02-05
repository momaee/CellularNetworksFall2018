clear
clc
aStruct = load('AdjMat');
A = aStruct.A;
specRadA = max(eig(A));
deltaMax = 0.8;
deltaMin = 3.9*0.2;
epsilonBar = 0.001;
DeltaTilda = max(epsilonBar,deltaMax);
epidTresh = (1-deltaMax)/specRadA;
betaMax = 4*epidTresh;
betaMin = 0.3*betaMax;
N = length(A);
rho = 4;

%%
numofIterations = 110;
phi = zeros(N,N,numofIterations);
phi(:,:,1) = 0.2*ones(N,N);
beta = zeros(N,numofIterations);
delta = zeros(N,numofIterations);
u = zeros(N,N,numofIterations);
for k = 1:numofIterations-1
    for i = 1:N
        phi(:,i,k+1) = phi(:,i,k) + rho * sum(  A(i,:).*(u(:,i,k)-u(:,:,k)), 2 ) ; 
        [beta(i,k+1), delta(i,k+1), u(:,i,k+1)] = updateLocalVariablesFmincon(i, phi(:,i,k+1), u(:,:,k), betaMin, betaMax, deltaMin, deltaMax, A(i,:), rho, DeltaTilda, epsilonBar);
        %[beta(i,k+1), delta(i,k+1), u(:,i,k+1)] = updateLocalVariables(i, phi(:,i,k+1), u(:,:,k), betaMin, betaMax, deltaMin, deltaMax, A(i,:), rho);
    end
end

%%
S = zeros(numofIterations,1);
for t = 1:numofIterations
    for i=1:N
        S(t) = S(t) + sum(sum((u(:,i,t)-u(:,:,t)).^2) .* A(i,:));
    end
end
figure
plot(S)
max(eig(    diag(beta(:,end))*A - diag(delta(:,end))    ))