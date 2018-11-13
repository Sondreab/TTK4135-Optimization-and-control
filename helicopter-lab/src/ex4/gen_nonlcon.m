function [c,ceq] = gen_nonlcon(z)

global N mx

alp = 0.2;
bet = 20;
lambda_t = 2*pi/3;

c = alp*exp(-bet*(z(1:mx:N*mx)-lambda_t).^2)-z(5:mx:N*mx);
ceq = [];

end
