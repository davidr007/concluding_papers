function ans = iseven(number)
% tells whether a number is even 
%
% INPUTS:
% number    [scalar - double] number to be tested
%
% OUTPUTS:
% ans       [scalar - ] 
%               1 => number is even
%               0 => number is odd
%               -1 => number is neither odd or even

tmp = number/2;
tmp2 = floor(tmp);
diff = tmp-tmp2;

if diff == 0 % number is even
    ans = 1;
elseif diff == 0.5  % number is odd
    ans = 0;
else
    ans = -1;
end

