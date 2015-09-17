% Finite-Difference Method for Linear Problems
% To approximate the solution of the boundary value problem
% y''=p(x) y' + q(x) y + r(x), for a6x 6 b, with y(a) = ? and y(b) = ?.
% INPUT: function handles p, q, r, boundary values alpha and beta, two
% boundary points a and b, number of interior points N where we want to
% approximate the values of the solution
% OUTPUT: approximate solution [y x] where x is the vector consisting of xj
% and y is the vector consisting of approximation to y(xj).
% ALL VECTORS ARE SUPPOSED TO BE COLUMN VECTORS.
function [y x] = LinearFiniteDifference (p, q, r, x0, xend, alpha, beta, N)
% x = []; % x is an empty matrix.
% y = [];
% Check for input error
if (xend < x0)
error(' b must be larger than a'); % print out error message then quit')
end
if( N < 1 )
    % then I intended N to be h
    h = N;
    x = (x0+h):h:(xend-h);
    x = x.';
    N = length( x );
else
    h = (xend-x0)/(N+1);
    x = x0 + (1:N)*h; % xj =x0+j h, j=1,2, ..., N
    x = x';
end
% Forming the matrix A and b to solve for Ay = b
b = -h^2 * r(x);    % make sure to input  p,q,r correctly, otherwise b will be incorrect
b(1) = b(1) + (1 + h/2 * p(x(1))) * alpha;
b(N) = b(N) + (1 - h/2 * p(x(N))) * beta;
A = zeros(N);
for i = 1:N
A(i, i) = 2 + h^2 * q(x(i));
if(i+1 < N)
    A(i, i+1) = -1 + h/2 * p(x(i));
    A(i+1, i) = -1 - h/2 * p(x(i+1));
end
A(N,N-1) = -1 - h/2 * p(x(N));
A(N-1,N) = -1 + h/2 * p(x(N-1));
end
% define y as the solution
y = A\b;
y = [alpha; y ; beta];
x = [x0 ; x ; xend];
end % end of function LinearFiniteDifference