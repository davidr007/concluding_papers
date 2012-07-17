function eps = calc_eps(stda, meana)
% computes the total fractional variance for a normal distribution with 
% mean = meana and standard deviation = stda. 

% first create a large set of normally distributed with zero mean and
% standard deviation = stda
tmp = stda*randn(1,1000000);
eps = tot_fracvar(tmp,meana);