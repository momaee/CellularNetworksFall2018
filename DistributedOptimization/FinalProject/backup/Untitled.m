clc
cvx_begin
    variable u
    minimize u^2-u
    subject to
       u >= 0
cvx_end
