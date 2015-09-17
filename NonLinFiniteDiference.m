% BRIAN SARRACINO-AGUILERA :-)
function [y x] = NonLinFiniteDiference(f, dfdyp, dfdy, xvalues, boundaryC, N, max, tol)
% w is length N. Then becomes length N+2 after inclding BCs
% f, dfdy, dfdyp: derivatives of the RHS of y''
% N: number of steps to take between a<=x<=b
% guess_w0: initial guess of solution for Newtons method
x0 = xvalues(1);
xend = xvalues(2);
alpha = boundaryC(1);
beta = boundaryC(2);
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

%%%%% initializing approximation w and wp to y and yp
w = alpha + (1:N)*(beta - alpha)/(xend - x0)*h;     % linear approximation
w = w.';
wp = zeros(N,1);    % wp: wprime using centered difference O(h^2)
F = zeros(N,1);
J = zeros(N,N);
%%%%%
count = 1;
while( count < max)
    % Build F(x) vector
    w = [alpha; w ; beta]; % w is length N+2
    for j=1:N
        wp(j) = ( w(j+2)-w(j) ) / (2*h);
        F(j) = -( w(j+2) - 2*w(j+1) + w(j) ) + h^2*f(x(j),w(j+1), wp(j) );
    end
    % Build Jacobian(x)

    for i = 1:N
    J(i, i) = 2 + h^2 * dfdy(x(i),w(i),wp(i) );
    if(i+1 < N)
        J(i, i+1) = -1 + h/2 * dfdyp(x(i),w(i),wp(i));
        J(i+1, i) = -1 - h/2 * dfdyp(x(i+1),w(i),wp(i));
    end
    J(N,N-1) = -1 - h/2 * dfdyp(x(N),w(N),wp(N));
    J(N-1,N) = -1 + h/2 * dfdyp(x(N-1),w(N),wp(N));
    end

    v = -J\F;
    % take boundaryC off w
    w(1) = [];
    w(end) = [];
    % wold = w
    w = w + v;
    
    if( norm(v) < tol ) %abs(wold - w) <= tol )
        fprintf('Beat the tolerance, success :-) \n')
        break
    end
    count = count + 1;
end
%
if( count == max)
    fprintf('Reached max iterations :-( \n')
end
y = w;
y = [alpha;y;beta];
x = [x0 ; x ; xend];

end % end of the function