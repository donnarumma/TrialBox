function [ U, V] = sim_sspca_noG( X, params, V) 
 
%==========================================================================
% *** Sparse Structured Principal Component Analysis (SSPCA) ***
% 
% [ U, V, params] = sspca( X, spG, params ) 
%
%
% solve: Min_{U,V} [ |X- U*V^t|_fro ^2 / ( 2*normalization ) + lambda * sum_k Omega(V^k) ]
%        under the constraint |U^k|_2 <= 1
%
% INPUT:
%
%   X:   data matrix of size n x p; contains n observations (the rows) described by p features (the columns) 
%
%   spG: matrix *** SPARSE *** of size p x ng; contains the description of the
%        group structure used to define the regularization term Omega.
%
%        spG(j,g) is zero if the variable j is not in the group indexed by g
%                 otherwise, accounts for the weight of the variable j in the group indexed by g
%
%        Example: if spG = speye(p), then Omega is the standard L_{params.normparam} norm (see more examples in the demo files)
%
%  params: structure that contains the different parameters of sspca (see the complete list of fields below, with their corresponding meaning).
%
%        Most important fields:
%
%                   params.lambda: regularization parameter (required)
%                   params.r     : number of dictionary elements (required)
%                   params.posV  : add positivity constraint on V
%                   params.posU  : add positivity constraint on U
%
%
% OUTPUT:
%
%   decomposition (U,V) such that U*V^t approximates well X 
%
%       U is n x params.r
%       V is p x params.r
%         
% For more details, see
%
%   (2009) R. Jenatton, G. Obozinski and F. Bach.  Structured sparse principal component analysis.
%
%
% Copyright (c) 2010 Rodolphe Jenatton. All rights reserved.
%==========================================================================

if any( isnan(X(:)) ) || any( isinf(X(:)) ) || ~isa(X,'double'),
    error('*** Bad data matrix X (contains some Inf or/and NaN values, or is not double) ***');
end

[n,p] = size(X);


%==========================================================================
% Process the params
%==========================================================================
if isfield(params,'lambda'), % Regularization parameter
    lambda = params.lambda;
else
    error('*** No lambda provided ***');
end
%--------------------------------------------------------------------------
if isfield(params,'r'), % Number of dictionary elements
    r = params.r;
else
    error('*** No dictionary size provided ***');
end
%--------------------------------------------------------------------------
% Number of different nonzero patterns among the r dictionary elements
% (in the paper, referred to as the "Shared structure" setting)
% For instance, if r=36 and m=12, there are 12 different nonzero patterns, each
% one represented by r/m=3 dictionary elements.
%
% The different nonzero patterns are organized as follows in V
% V(,1)...V(:,r/m)           : first nonzero pattern
% V(,r/m+1)...V(:,2*r/m)     : second nonzero pattern
% ...
% V(,r*(m-1)/m + 1)...V(:,r) : m-th nonzero pattern
if isfield(params,'m') && ( mod(r,params.m)==0 ), 
    m           = params.m;    
else
    m = r;
end
%--------------------------------------------------------------------------
if isfield(params,'max_it') % Maximum number of iterations (i.e., loops over U,V, Eta)
    max_it = params.max_it;
else
    max_it        = 300; 
    params.max_it = max_it;
end
%--------------------------------------------------------------------------
if isfield(params,'it0'),  % Display cost function every it0 iterations (for no verbose, set it0 to Inf)
    it0 = params.it0;
else
    it0        = 30; 
    params.it0 = it0;
end
%--------------------------------------------------------------------------
if isfield(params,'max_it_U'), % Maximum number of iterations for the update of U
    params_U.max_it = params.max_it_U;
else
    params_U.max_it = 5; 
    params.max_it_U = params_U.max_it;
end
%--------------------------------------------------------------------------
if isfield(params,'max_it_V'),  % Maximum number of iterations for the update of V
    params_V.max_it = params.max_it_V;
else
    params_V.max_it = 5;
    params.max_it_V = params_V.max_it;
end

%--------------------------------------------------------------------------
% Stopping criterion: the algorithm stops if the relative decrease in the
% cost function becomes smaller than min_delta_cost
if isfield(params,'min_delta_cost'),
    min_delta_cost = params.min_delta_cost;
else
    min_delta_cost        = 0.1;
    params.min_delta_cost = min_delta_cost;
end
%--------------------------------------------------------------------------
if isfield(params,'initialization_seed'),% Control the seed for random initializations
    initialization_seed = params.initialization_seed;
else
    initialization_seed        = 0;
    params.initialization_seed = initialization_seed;
end
%--------------------------------------------------------------------------
if isfield(params,'normalization'), % Normalization of the cost function
    normalization = params.normalization;
else
    normalization = norm( X, 'fro' )^2 / 2;
end
params_V.normalization = normalization;
%--------------------------------------------------------------------------
if isfield(params,'posU') && params.posU,% Add positivity constraint on U
    params_U.pos = true;
else
    params_U.pos = false;
end
%--------------------------------------------------------------------------
if isfield(params,'posV') && params.posV,% Add positivity constraint on V
    params_V.pos = true;
else
    params_V.pos = false;
end
%--------------------------------------------------------------------------
% Parameter that controls the convexity of Omega
% Omega is a mixed norm "L_normparam - L_2"
% If normparam is in (0,1), then Omega is concave (sparsity is yielded more agressively)
if isfield(params,'normparam'), 
    params_Eta.normparam = params.normparam;
else
    params_Eta.normparam = 0.5;
end
params.normparam = params_Eta.normparam;
%--------------------------------------------------------------------------
% The optimization is based on a variational reformulation of Omega 
% The parameter epsilon deals with the nonsmoothness that appears in this reformulation
% (for more detail, see the aforementioned paper)
%
% Note that, when displaying the cost function, "cost w/o smoothing" refers
% to the objective function *without* the variational reformulation of Omega
if isfield(params,'epsilon'),
    params_Eta.epsilon = params.epsilon;
else
    params_Eta.epsilon = 1e-7;
end
%==========================================================================
% Intialization
%==========================================================================
if isfield(params,'U') && ( size(params.U,1) == n ) && ( size(params.U,2) == r ), % We can restart from a given U
    U = params.U;
else
    past_seed = randn('state');
    randn('state',initialization_seed);
    
    if isfield(params,'posU') && params.posU,
        U = abs(randn(n,r));
    else
        U = randn(n,r);
    end
    
    randn('state',past_seed);
end

% if isfield(params,'V') && ( size(params.V,1) == p ) && ( size(params.V,2) == r ), % We can restart from a given V
%     V = params.V;
% else
%     past_seed = randn('state');
%     randn('state',initialization_seed);
%     
%     if isfield(params,'posV') && params.posV,
%         V = abs(randn(p,r));
%     else
%         V = randn(p,r);
%     end
%     
%     randn('state',past_seed);
% end
%--------------------------------------------------------------------------
% We make sure that (U,V) are feasible points
%--------------------------------------------------------------------------
if isfield(params,'posU') && params.posU,
    U( U < 0 ) = 0.0;
end
U = bsxfun( @times, U, 1./sqrt( sum( U.^2, 1 ) ) );
% if isfield(params,'posV') && params.posV,
%     V( V < 0 ) = 0.0;
% end
%==========================================================================

t          = 1;
delta_cost = Inf;
past_cost  = Inf;

%==========================================================================

while ( t<=max_it ) && ( delta_cost > min_delta_cost ),
    
    display_cost_function;
    
    U               = mexUpdateU( V, X, U, params_U);
    
    V               = mexUpdateV( V, X, U, reshape(inv_Zeta,[p,r]), lambda, params_V);
    
    t = t + 1;
    
end 
%==========================================================================
% Nested functions
%==========================================================================
% Note that, when displaying the cost function, "cost w/o smoothing" refers
% to the objective function *without* the variational reformulation of Omega
    function display_cost_function
    
        if ( mod(t,it0) == 0 ),

          loss = norm( X - U*V', 'fro' )^2 / (2*normalization);
            
          Omega_exact  = 0; %mexGetOmega(reshape(V,[r/m*p,m]),spG_squared,params_Eta);
          Omega_approx = 0; %get_Omega_from_eta;

          
%           cost_exact  = loss + lambda*Omega_exact;
%           cost_approx = loss + lambda*Omega_approx;
           cost_exact  = loss;
           cost_approx = loss;

          if ( t > it0 ),

            delta_cost = 100*(past_cost - cost_approx)/past_cost;

          end

          past_cost  = cost_approx;

          fprintf( '*** Iteration: %3.0d, cost: %6.4f [cost w/o smoothing: %6.4f], cost variation: %6.4f ***\n',...
                                     t,  cost_approx, cost_exact, delta_cost) 

        end
        
    end
%==========================================================================
    function output = get_Omega_from_eta
        
        beta           = params_Eta.normparam / ( 2 - params_Eta.normparam );
        
        norm_eta       = sum( sum( Eta .^ beta, 1 ) .^ (1/beta) );

        quadratic_term = sum( sum( reshape(V.^2,[r/m*p,m]) .* inv_Zeta, 1) );
        
        output         = 0.5 * quadratic_term + 0.5 * norm_eta;
        
    end
%==========================================================================
end
