% mixed_ca_lanczos_bounds.m
% This file runs mixed precision s-step Lanczos algorithm and measures
% quantities related to the normality, loss of orthogonality, column error,
% and deviation of AV and T, as well as the bounds on these quantities. All
% computations are performed in double precision except for computing and
% applying the Gram matrix G, which is done in quadruple precision (using
% the Advanpix toolbox). 
%
% Input:
%   A : a square, symmetric coefficient matrix
%   v : the starting vector for the Lanczos method
%   s : the parameter s in s-step (s>0)
%   basis_info : a struct containing a quantity 'type' which defines the
%   polynomial basis to be used. Options are 'monomial', 'newton', or
%   'chebyshev'
%   options : a struct containing a quantity 'xlim', which gives the number
%   of iterations to perform, and a quantity 'name', which is used for
%   naming the output file
%
%
% Last edited by: Erin Carson, 2021
%
function [results, basis_info] = mixed_ca_lanczos_bounds(A, v, s, basis_info, options)

% Measure quantities needed for the error bounds
sigma = norm(full(A));
theta = norm(full(abs(A)),2)/sigma;

% Compute/set basis parameters
[alp,bet,gam,~] = basisparams (s+1, A, basis_info);

% Store the basis parameters used for output
basis_info.alp = alp;
basis_info.bet = bet;
basis_info.gam = gam;

% Construct and store the "change of basis matrix"
Tt = diag(alp(1:s),0) + diag(bet(1:s-1),1) + diag(gam(1:s-1),-1);
Tt = [Tt;zeros(1,s-1),gam(s)];
Tt = [Tt,zeros(s+1,1)];
T  = sparse ( [Tt, zeros(s+1,s+1); zeros(s+1,s+1), Tt] );
basis_info.T = T;


%Set tau parameter for bounds
tau = norm(full(abs(T)))/sigma;

% Get size of matrix and maximum number of nonzeros per row
n = size(A,1);
N = max(sum(A~=0,2));

% Set initial values for vectors (ensure unit starting vector)
v = v./norm(v);
v0 = v;
u0 = A*v0;
v(:,1) = v0;
u(:,1) = u0;

% Initialize quantities
beta(1) = 0;
k = 0;
m = 0;
Gamma = 0;


% Begin the iterations!
while m < options.xlim

    % Compute s-step Krylov bases
    V = computeBasis(A,v(:,m+1),s+1,basis_info);
    U = computeBasis(A,u(:,m+1),s+1,basis_info);
    
    % Record basis rank if basis is rank deficient
    if (rank([V,U])<2*s+2)
        basisrank = rank([V,U]);
    end
    
    % Set coordinate vectors
    vcoeff = [[1;zeros(2*s+1,1)],zeros(2*s+2,s)];
    ucoeff = [[zeros(s+1,1);1;zeros(s,1)],zeros(2*s+2,s)];
    
    % Compute Gram matrix in quadruple precision
    G = mp([V,U],34)'*mp([V,U],34);

    % Set Gamma_k parameter for bound
    if( k>0)
        Gammat = (norm(pinv([V,U]),2)*norm(abs([V,U]),2));
    else
        Gammat = (norm(pinv([V]),2)*norm(abs([V]),2));
    end
    if(Gammat > Gamma)
        Gamma = Gammat;
    end
    
    % Set quantities in bounds
    e0 = 2*eps*(6*s+11)*Gamma;
    e1 = eps*( (N+2*s+5)*theta + (4*s+9)*tau + 10*s+16)*Gamma;
    
    %Inner iterations
    for j = 1:s
        
        % Increment global iteration count
        m = m + 1;
        
        % Record Gamma for bounds
        results.Gamval(m) = Gamma;
        
        % Update iteration vectors
        alpha(m) = (vcoeff(:,j)'*double(mp(G,34)*mp(ucoeff(:,j),34))); % apply G to vector in quadruple precision
        wcoeff = ucoeff(:,j) - alpha(m)*vcoeff(:,j);
        w = ([V,U]*wcoeff);
        beta(m+1) = sqrt(wcoeff'* double(mp(G,34)*mp(wcoeff,34))); % apply G to vector in quadruple precision
        vcoeff(:,j+1) = wcoeff/beta(m+1);
        v(:,m+1) = ([V,U]*vcoeff(:,j+1));
        ucoeff(:,j+1) = T*vcoeff(:,j+1) -beta(m+1)*vcoeff(:,j);
        u(:,m+1) = ([V,U]*ucoeff(:,j+1));
        
        % Store measured quantities and their bounds
        if(m>1)
            results.orthactual(m) = mp(beta(m+1),64)*abs(mp(v(:,m),64)'*mp(v(:,m+1),64));
            results.orthbound(m) =  e0*sigma;
            results.normactual(m) = abs(mp(mp(v(:,m+1),64)'*mp(v(:,m+1),64)-1,64));
            results.normbound(m) = e0/2;
            results.colactual(m) = norm(mp( mp(A,64)*mp(v(:,m),64)-(mp(beta(m+1),64)*mp(v(:,m+1),64) + mp(beta(m),64)*mp(v(:,m-1),64) + mp(alpha(m),64)*mp(v(:,m),64)),64));
            results.colbound(m) = e1*sigma;
            results.diffactual(m) = abs(mp(beta(m+1),64)^2 + mp(alpha(m),64)^2 + mp(beta(m),64)^2 - norm(mp(A,64)*mp(v(:,m),64))^2);
            results.diffbound(m) = 2*m*(3*e0+2*e1)*sigma^2;
            
        end
        
    end
    
    % Increment outer iteration count
    k = k+1;
end

% Store tridiagonal T matrix
results.T = diag(alpha(1:end),0)+diag(beta(2:end-1),1)+ diag(beta(2:end-1),-1);

end
