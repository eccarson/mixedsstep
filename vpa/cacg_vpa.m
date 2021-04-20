% cacg_vpa.m
% Run the s-step CG algorithm (CACG) to solve Ax=b
%
% Input:
%   A : square, sparse matrix with dimension n
%   b : right hand side of system to solve, Ax=b; vector of dimension n
%   s : number of inner-loop iterations per outer loop; the "s" in "s-step
%   methods"
%   x0 : initial guess for solution, vector of dimension n
%   maxits : maximum number of iterations to complete before returning
%   basis_type: string denoting which basis to use. Acceptable values are
%   'monomial', 'newton', or 'chebyshev'. If something besides these strings
%   entered, will default to using monomial basis.
%
%   Output:
%   results struct stores:
%       r_exact_norm : 2-norm of true residual computed in each iteration
%       (results.r_exact_norm)
%       r_comp_norm : 2-norm of computed residual computed in each iteration
%       (results.r_comp_norm)
%       x : approximate solution computed in each iteration
%       (results.x)
%
% Last edited by: Erin Carson, 2021
%

function results = cacg_vpa(A, b, s, x0, basis_info, options)

% Size of matrix
N = size(A,1);

% Set initial values for vectors
r0 = b - A*x0;
p0 = r0;
x(:,1)  = x0;
r(:,1)  = r0;
p(:,1)  = p0;

% Set outer loop iteration count to 0
k = 0;

% Set total number of iterations to 0
its = 0;

% Initialize initial true and computed residuals and approximate solution
results.r_exact_norm(1) = norm(b-A*x0);
results.r_comp_norm(1) = norm(r0);
results.error_A_norm(1) = 1;
results.x = x0;

% Measure quantities needed for the error bounds
sigma = norm(full(A));
theta = norm(full(abs(A)),2)/sigma;

% Compute/set basis parameters
[alp,bet,gam, T] = basisparams(s, A, basis_info);

% Store the basis parameters used for output
basis_info.alp = alp;
basis_info.bet = bet;
basis_info.gam = gam;


Tt  = sparse([T, zeros(s+1,s+1); zeros(s,s+1), T(1:s,1:(s-1)),zeros(s,1)]);
basis_info.T = Tt;

% Initialize quantity that stores maximum basis condition number
gammax = 0;

% Begin the iterations
while its < options.xlim
    
    % Compute Krylov basis with starting vector p
    P = computeBasis(A,p(:,s*k+1),s+1,basis_info);
    
    % Compute Krylov basis with starting vector r
    R = computeBasis(A,r(:,s*k+1),s,basis_info);
     
    % Initialize CG coordinate vectors for current outer loop
    p_c = [[1;zeros(2*s,1)],zeros(2*s+1,s)];
    r_c = [[zeros(s+1,1);1;zeros(s-1,1)],zeros(2*s+1,s)];
    x_c = zeros(2*s+1,s+1);
    
    % Compute Gram matrix
    G = [P,R]'*[P,R];
    
    if (k>0)
        gammax = max(cond(G), gammax);
    end
    
    % Begin s inner iterations
    for j = 1:s
        
        if (its>=options.xlim)
            break;
        end
        
        % Increase iteration count
        its = its + 1;
        
        % Compute scalar alpha using Gram matrix for inner products
        alpha(its) = (r_c(:,j)'*G*r_c(:,j)) / (p_c(:,j)'*G*Tt*p_c(:,j));
        
        % Update x coordinate vector
        x_c(:,j+1) = x_c(:,j) + alpha(its)*p_c(:,j);
        
        % Perform basis change to compute x vector in standard basis (note we wouldn't need to do this
        % in the inner loop in practice)
        x(:,s*k+j+1) = [P,R]*x_c(:,j+1) + x(:,s*k+1);
        
        % Update r coordinate vector
        r_c(:,j+1) = r_c(:,j) - alpha(its)*Tt*p_c(:,j);
        
        % Perform basis change to compute r vector in standard basis (note we wouldn't need to do this
        % in the inner loop in practice)
        r(:,s*k+j+1)  = [P,R]*r_c(:,j+1);
        
        % Compute scalar beta using Gram matrix for inner products
        beta(its) = (r_c(:,j+1)'*G*r_c(:,j+1))/ (r_c(:,j)'*G*r_c(:,j));
        
        % Update p coordinate vector
        p_c(:,j+1) = r_c(:,j+1) + beta(its)*p_c(:,j);
        
        % Perform basis change to compute p vector in standard basis (note we wouldn't need to do this
        % in the inner loop in practice)
        p(:,s*k+j+1)  = [P,R]*p_c(:,j+1);
        
        % Compute and store true residual norm (note we wouldn't do this in
        % the inner loop in practice)
        results.r_exact_norm(its+1) = norm(b-A*x(:,s*k+j+1));
        
        % Compute and store computed residual norm (note we wouldn't do this, at least in this way, in
        % the inner loop in practice)
        results.r_comp_norm(its+1) = norm(r(:,s*k+j+1));
        
        % Compute the relative error in the A-norm 
        err = vpa(x(:,its+1) - options.truesol,34);
        num = sqrt(double(vpa(err',34)*vpa(A,34)*vpa(err,34)));
        denom = sqrt(double(vpa(options.truesol',34)*vpa(A,34)*vpa(options.truesol,34)));
        results.error_A_norm(its+1) = double(num/denom);results.error_A_norm(its+1) = double(num/denom);
   
        
        % Store current approximate solution
        results.x = x(:,s*k+j+1);
        
    end
    
    % Increment k, outer iteration index
    k = k+1;
    
end

results.gammax = double(gammax);



