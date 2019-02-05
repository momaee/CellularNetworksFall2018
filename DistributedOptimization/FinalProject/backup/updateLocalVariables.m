function [ b, d, u ] = updateLocalVariables( i, phi_i, u_k , betaMin, betaMax, deltaMin, deltaMax, Arow_i, rho)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
n = length(phi_i);
cvx_solver gurobi
cvx_begin gp
    variable b nonnegative
    variable d nonnegative
    variable u(n) 
    f = (pow_p(b,-1)-betaMax^(-1))/(betaMin^(-1)-betaMax^(-1));
    g = (pow_p(1-d,-1)-(1-deltaMin)^(-1))/((1-deltaMax)^(-1)-(1-deltaMin)^(-1));
    h = rho * sum(sum(pow_pos((repmat(u,1,n) - (u_k(:,i)+u_k)/2 ), 2)).*Arow_i) + (phi_i')*u;
    minimize f+g+h
    subject to
        b*(Arow_i*u)+d*u(i) <= u(i);
        prod(u) == 1;
        deltaMin <= d <= deltaMax;
        betaMin <= b <= betaMax;
cvx_end

end

