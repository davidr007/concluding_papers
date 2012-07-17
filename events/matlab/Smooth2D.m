function A = smooth2D(B)
% Smooth2D smooths a matrix with a 3 x 3 running window. 
%
% David Robinson
% APRIL 2003

% Orientation as perceived by program
% Note that 5 represents the point of interest
% 1   4   7
% 2   5   8
% 3   6   9


[n,m] = size(B);
A = zeros(n,m);

% All points on top row


% All central points
A(2:n-1,2:m-1)= ...
    B(1:n-2,1:m-2) + ...       % 1   i.e. one shift up and one shift to the left    
    B(2:n-1,1:m-2) + ...       % 2   i.e. one shift to the left
    B(3:n,1:m-2)   + ...       % 3   i.e. one shift down and one shift to the left  
    B(1:n-2,2:m-1) + ...       % 4   i.e. one shift up
    B(2:n-1,2:m-1) + ...       % 5   i.e. central position     
    B(3:n,2:m-1)   + ...       % 6   i.e. one shift down
    B(1:n-2,3:m)   + ...       % 7   i.e. one shift up and one shift to the right
    B(2:n-1,3:m)   + ...       % 8   i.e. one shift to the right
    B(3:n,3:m);                % 9   i.e. one shift down and one shift to the right
    
 % Points in top centre
 A(1,2:m-1)= ...
    B(1,1:m-2) + ...           % 1   i.e. one shift up (take copy of actual) and one shift to the left    
    B(1,1:m-2) + ...           % 2   i.e. one shift to the left
    B(2,1:m-2) + ...           % 3   i.e. one shift down and one shift to the left  
    B(1,2:m-1) + ...           % 4   i.e. one shift up (take copy of actual)
    B(1,2:m-1) + ...           % 5   i.e. central position     
    B(2,2:m-1) + ...           % 6   i.e. one shift down
    B(1,3:m)   + ...           % 7   i.e. one shift up (take copy of actual) and one shift to the right
    B(1,3:m)   + ...           % 8   i.e. one shift to the right
    B(2,3:m);                  % 9   i.e. one shift down and one shift to the right

% Points in Bottom centre
A(n,2:m-1)= ...
    B(n-1,1:m-2) + ...         % 1   i.e. one shift up and one shift to the left    
    B(n,1:m-2) + ...           % 2   i.e. one shift to the left
    B(n,1:m-2)   + ...         % 3   i.e. one shift down (take copy of actual) and one shift to the left  
    B(n-1,2:m-1) + ...         % 4   i.e. one shift up
    B(n,2:m-1) + ...           % 5   i.e. central position     
    B(n,2:m-1)   + ...         % 6   i.e. one shift down (take copy of actual)
    B(n-1,3:m)   + ...         % 7   i.e. one shift up and one shift to the right
    B(n,3:m)   + ...           % 8   i.e. one shift to the right
    B(n,3:m);                  % 9   i.e. one shift down (take copy of actual) and one shift to the right

% Points in the left centre
A(2:n-1,1)= ...
    B(1:n-2,1) + ...           % 1   i.e. one shift up and one shift to the left (take copy of actual)   
    B(2:n-1,1) + ...           % 2   i.e. one shift to the left  (take copy of actual)
    B(3:n,1)   + ...           % 3   i.e. one shift down and one shift to the left  (take copy of actual)
    B(1:n-2,1) + ...           % 4   i.e. one shift up
    B(2:n-1,1) + ...           % 5   i.e. central position     
    B(3:n,1)   + ...           % 6   i.e. one shift down
    B(1:n-2,2)   + ...         % 7   i.e. one shift up and one shift to the right
    B(2:n-1,2)   + ...         % 8   i.e. one shift to the right
    B(3:n,2);                  % 9   i.e. one shift down and one shift to the right

% Points in the right centre
A(2:n-1,m)= ...
    B(1:n-2,m-1) + ...         % 1   i.e. one shift up and one shift to the left    
    B(2:n-1,m-1) + ...         % 2   i.e. one shift to the left
    B(3:n,m-1)   + ...         % 3   i.e. one shift down and one shift to the left  
    B(1:n-2,m) + ...           % 4   i.e. one shift up
    B(2:n-1,m) + ...           % 5   i.e. central position     
    B(3:n,m)   + ...           % 6   i.e. one shift down
    B(1:n-2,m)   + ...         % 7   i.e. one shift up and one shift to the right (take copy of actual)
    B(2:n-1,m)   + ...         % 8   i.e. one shift to the right (take copy of actual)
    B(3:n,m);                  % 9   i.e. one shift down and one shift to the right (take copy of actual)

% Top left point
A(1,1) = 4*B(1,1)+2*B(1,2)+2*B(2,1)+B(2,2);

% Bottom left point
A(n,1) = 4*B(n,1)+2*B(n-1,1)+2*B(n,2)+B(n-1,2);

% Top right point
A(1,m) = 4* B(1,m)+2*B(2,m)+2*B(1,m-1)+B(2,m-1);

% Bottom right point
A(n,m) = 4*B(n,m) + +2*B(n,m-1) + 2*B(n-1,m) + B(n-1,m-1);

A = A./9;  % compute the value as an average