% mixed_ca_lanczos_bounds_tight.m
% This file runs mixed precision s-step Lanczos algorithm and measures
% quantities related to the normality, loss of orthogonality, column error,
% and deviation of AV and T, as well as the (tight) bounds on these quantities. All
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
function [results, basis_info] = mixed_ca_lanczos_bounds_tight(A, v, s, basis_info, options)

% Measure quantities needed for the error bounds
sigma = norm(full(A));
theta = norm(full(abs(A)),2)/sigma;

Gamma_diffbound(1) = 0;

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
    
    Y = [V,U];
    
    
    % Record basis rank if basis is rank deficient
    if (rank(Y)<2*s+2)
        basisrank = rank([V,U]);
    end
    
    % Set coordinate vectors
    vcoeff = [[1;zeros(2*s+1,1)],zeros(2*s+2,s)];
    ucoeff = [[zeros(s+1,1);1;zeros(s,1)],zeros(2*s+2,s)];
    
    % Compute Gram matrix in quadruple precision
    G = mp(Y,34)'*mp(Y,34);
    Gabs = abs(mp(Y,34))'*abs(mp(Y,34));
    
    % Set Gamma_k parameter for bound
    if( k>0)
        Gammat = (norm(pinv(Y),2)*norm(abs(Y),2));
    else
        Gammat = (norm(pinv([V]),2)*norm(abs([V]),2));
    end
    if(Gammat > Gamma)
        Gamma = Gammat;
    end
    
    % Set quantities in bounds
    e0 = 2*eps*(9*s+14)*Gamma;
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
        w = (Y*wcoeff);
        beta(m+1) = sqrt(wcoeff'* double(mp(G,34)*mp(wcoeff,34))); % apply G to vector in quadruple precision
        vcoeff(:,j+1) = wcoeff/beta(m+1);
        v(:,m+1) = (Y*vcoeff(:,j+1));
        ucoeff(:,j+1) = T*vcoeff(:,j+1) -beta(m+1)*vcoeff(:,j);
        u(:,m+1) = (Y*ucoeff(:,j+1));
        
%         
%         Yw(m) = norm(abs(Y)*abs(wcoeff)) /  norm(Y*wcoeff);
%         Yv(m) = sqrt(abs(vcoeff(:,j))'*Gabs*abs(vcoeff(:,j))) / sqrt(abs(vcoeff(:,j)'*G*vcoeff(:,j)));
%         %Yv1(m) = (abs(vcoeff(:,j+1))'*Gabs*abs(vcoeff(:,j+1))) / (abs(vcoeff(:,j+1)'*G*vcoeff(:,j+1)));
%         Yu(m) = sqrt(abs(ucoeff(:,j))'*Gabs*abs(ucoeff(:,j))) /  sqrt(abs(ucoeff(:,j)'*G*ucoeff(:,j)));
%        % Yuv(m) = abs(vcoeff(:,j))'*Gabs*abs(ucoeff(:,j)) /  (abs(ucoeff(:,j)'*G*vcoeff(:,j)));
%         YBv(m) = sqrt(abs(T)*abs(vcoeff(:,j)))'*Gabs*(abs(T)*abs(vcoeff(:,j))) /  sqrt(norm(abs(full(T)))^2*(vcoeff(:,j)'*G*vcoeff(:,j)));
%        
%         
%         
%         % Store measured quantities and their bounds
%         if(m>1)
%             
%             Gamma_norm(m) = Yw(m);
%             Gamma_orth(m) = max([Yw(m-1), Yw(m), Yv(m),Yu(m)]);
%             
%             Gamma_colbound(m) = max([Yw(m), Yu(m), YBv(m), Yv(m), Yv(m-1)]);
% 
%             Gamma_diffbound(m) = max(Gamma_diffbound(m-1), max([Yw(m), YBv(m), Yv(m), Yu(m), Yw(m-1), YBv(m-1), Yv(m-1), Yu(m-1)]));% max([Yw(m),Yv(m),Yu(m),Yuv(m),Yv(m-1)]); %
% 
%             
        Yw(m) = (abs(wcoeff)'*Gabs*abs(wcoeff)) /  (abs(wcoeff'*G*wcoeff));
        Yv(m) = (abs(vcoeff(:,j))'*Gabs*abs(vcoeff(:,j))) / (abs(vcoeff(:,j)'*G*vcoeff(:,j)));
        Yv1(m) = (abs(vcoeff(:,j+1))'*Gabs*abs(vcoeff(:,j+1))) / (abs(vcoeff(:,j+1)'*G*vcoeff(:,j+1)));
        Yu(m) = (abs(ucoeff(:,j))'*Gabs*abs(ucoeff(:,j))) /  (abs(ucoeff(:,j)'*G*ucoeff(:,j)));
        Yuv(m) = abs(vcoeff(:,j))'*Gabs*abs(ucoeff(:,j)) /  (abs(ucoeff(:,j)'*G*vcoeff(:,j)));
        YBv(m) = (abs(T)*abs(vcoeff(:,j)))'*Gabs*(abs(T)*abs(vcoeff(:,j))) /  (norm(abs(full(T)))^2*(vcoeff(:,j)'*G*vcoeff(:,j)));
        
        
        % Store measured quantities and their bounds
        if(m>1)
            
            Gamma_norm(m) = Yw(m);
            Gamma_orth(m) = max([Yw(m-1), Yw(m), Yv(m),Yu(m)]);
            Gamma_colbound(m) = max([Yw(m), Yu(m), YBv(m), Yv(m), Yv(m-1)]);
            Gamma_diffbound(m) = max(Gamma_diffbound(m-1), max([Yw(m), YBv(m), Yv(m), Yu(m), Yw(m-1), YBv(m-1), Yv(m-1), Yu(m-1)]));
            
            results.Gamma_norm(m)=Gamma_norm(m);
            results.Gamma_orth(m)=Gamma_orth(m);
            results.Gamma_colbound(m)=Gamma_colbound(m);
            results.Gamma_diffbound(m)=Gamma_diffbound(m);
            
            results.orthactual(m) = mp(beta(m+1),64)*abs(mp(v(:,m),64)'*mp(v(:,m+1),64));
            results.orthbound(m) =  e0*sigma;
            results.orthboundtight(m) = 2*eps*(9*s+14)*sqrt(Gamma_orth(m))*sigma;
            
            results.normactual(m) = abs(mp(mp(v(:,m+1),64)'*mp(v(:,m+1),64)-1,64));
            results.normbound(m) = e0/2;
            results.normboundtight(m) = eps*(9*s+14)*sqrt(Gamma_norm(m));
            
            results.colactual(m) = norm(mp( mp(A,64)*mp(v(:,m),64)-(mp(beta(m+1),64)*mp(v(:,m+1),64) + mp(beta(m),64)*mp(v(:,m-1),64) + mp(alpha(m),64)*mp(v(:,m),64)),64));
            results.colbound(m) = e1*sigma;
            results.colboundtight(m) = eps*( (N+2*s+5)*theta + (4*s+9)*tau + 10*s+16)*sqrt(Gamma_colbound(m))*sigma;
            
            results.diffactual(m) = abs(mp(beta(m+1),64)^2 + mp(alpha(m),64)^2 + mp(beta(m),64)^2 - norm(mp(A,64)*mp(v(:,m),64))^2);
            results.diffbound(m) = 2*m*(3*e0+2*e1)*sigma^2;
            results.diffboundtight(m) = 2*m*(3*2*eps*(9*s+14)*sqrt(Gamma_diffbound(m))+2*eps*( (N+2*s+5)*theta + (4*s+9)*tau + 10*s+16)*sqrt(Gamma_diffbound(m)))*sigma^2;
            
            
        end
        
    end
    
    % Increment outer iteration count
    k = k+1;
end


end
