function P = phi(L,n,options)
%PHI Initializes \phi_i(L) for i=0,...,n for matrix or vector L.
%
% NOTE: This code was developed to experiment with different methods for initializing
% phi functions. For a more robust PADE' based phi function initilizer see 
% EXPINT (http://www.math.ntnu.no/num/expint/matlab.php)
%
% Parameters
%   L           Matrix/Vector
%               Square matrix or array of scaler values. 
%  
%   n           Integer
%               highest phi function to initialize.
%
%   options     struct (optional)
%               By Default scaling and squaring is used when L is a Matrix. options 
%               can be used to switch to contour integral.
%
%               Fields:
%                   'algorithm' - 's' for scaling & squaring, "c" for contour integral
%
%                   'M'         - If options.algorithm='s' then M denotes number of 
%                                 discretization points for contour. 
%                               - If options.algorithm='c' then M denotes number of 
%                                 Taylor terms to use in sum. 
%
%                   'R'         - contour radius
%
%                   'z0'        - contour center               
% Output
%   P       Matrix or List
%           If L is a Matrix, phi{i} = phi_{i-1}(L) where size(P{i}) = size(L)
%           If L is a vector, phi(i,j) = \phi_{i-1}(L(j)) where size(P) = [n+1,size(L)]

if(nargin == 2)
    options = struct('algorithm','s');
end

% Scaler case
if(isvector(L))
    P = phi_scaler(L,n);
% Matrix Case    
elseif(ismatrix(L))
    if(strcmp(options.algorithm,'s'))
        P = phi_matrix_scalesquare(L,n,options);
    elseif(strcmp(options.algorithm,'c'))
        P = phi_matrix_contour(L,n,options);
    end
else
    fprintf('Invalid Input. L must be Vector or a Matrix\n');
end

end

% ==============================================================================

% Scaler Algorithms

% ==============================================================================


function P=phi_scaler(L,n)
%P Evaluates \phi_i(L) for i=0,...,n and scaler/vector L using 
% Recursion relation and Cauchy Integral Formula.
%
% Input Parameters
%   L - (array) of \Lambda values L = [\Lambda_1, \Lambda_2, ..., \Lambda_N]
%   n - highest phi function to initialize.
% Output Parameters
%   P - (array) array of dimensions (n+1) x length(L). P(i,j) = \phi_{i-1}(L(j))

P = zeros(n+1,length(L));

% Cauchy Integral Settings
tol = 1;            % if abs(L(i)) < tol then use contour integral, otherwise direct formula
R = 2*tol;          % Contour radius
M = 32;             % Number of descritization points
z   = R*exp(2*pi*1i*(0:(M-1))/M);  % circular contour 

% Split small and large L
direct_index  = abs(L) > tol;
contour_index = ~direct_index;

% direct formula
Ld = L(direct_index);
if(~isempty(Ld))
    P(1,direct_index) = exp(Ld);
    for i=2:n+1
        P(i,direct_index) = (P(i-1,direct_index) - 1/factorial(i-2))./Ld;
    end
end

% cauchy integral formula
Lc = L(contour_index);
if(~isempty(Lc))
    Lc = reshape(Lc,1,length(Lc));      % reshape into 1xNc vector
    zc = reshape(z,1,1,M);              % reshape into 1x1xM vector
    Lz = bsxfun(@plus,Lc,zc);           % form 1xNcxM vector where Lz(1,i,j) = Lc(i) + ac(j)
    PC = zeros(n+1,length(Lc),M);

    PC(1,:,:) = exp(Lz);
    for i=2:n+1
        PC(i,:,:) = (PC(i-1,:,:) - 1/factorial(i-2))./Lz;
    end
    P(:,contour_index) = (1/M)*sum(PC,3);
    % remove imaginary part from real L
    real_index = contour_index & imag(L) == 0; 
    P(:,real_index) = real(P(:,real_index));
end

end

% ==============================================================================

% Matrix Algorithms

% ==============================================================================

function P = phi_matrix_contour(A,n,params)
%PHI_matrix_contour Evaluates \phi_i(A) for matrix A and i=0,...,n using Cauchy Integral Formula.
%
% Input Parameters
%   A       - (matrix) matrix agrument
%   n       - highest phi function to initialize.
%   params  - (struct,optional) if specificed must contain fields: 
%             "M"  - number of descrete points
%             "R"  - contour radius
%             "z0" - contour center
% Output Parameters
%   P - (cell array) nx1 cell array of dimensions where  P{i} = \phi_{i-1}(A)

if(nargin == 3 && sum(isfield(params,{'M','z0','R'})) == 3)
    % Explicit Parameters
    M  = params.M;
    z0 = params.z0;
    R  = params.R;
else
    % Automatic parameters
    % For details see H.A Ashi (2007, p.477) "Comparison of methods for evaluating functions of a matrix exponential"
    e = 1e-13;                          % relative error we hope achieve
    ev = eig(A);                        % eigenvalues of A
    z0 = (max(ev) + min(ev))/2;         % contour center
    R  = (max(ev) - z0)+10;             % contour radius
    M  = max(60,ceil(exp(1)*R - log(e) - 3/2*log(R) - 3/2 - log(sqrt(2*pi)))); % Number of discrete points for trap rule:
end

P = cell(n+1,1); % struct to store phi functions

% Initialize scaler phi values
z = R*exp(2*pi*1i*(0:(M-1))/M);
p = phi_scaler(z + z0,n);

% Init matrix phi values
I = eye(size(A));
for i=1:M
    Ai = (I + (z0*I - A)/z(i))\I;
    for j=1:n+1
        if(i==1)
            P{j} = p(j,i)/M*Ai;
        else
            P{j} = P{j} + p(j,i)/M*Ai;
        end
    end    
end

% cast to real
if(isreal(A))
    for i=1:n+1
        P{i} = real(P{i});
    end
end
end


function P = phi_matrix_scalesquare(A,n,params)
%INITPHI Evaluates \phi_i(A) for matrix A and i=0,...,n using Taylor Based Scaling and Squaring Procedure
% Input Parameters
%   A       - (matrix) matrix agrument
%   n       - highest phi function to initialize.
%   params  - (struct), optional. contains fields: 
%             "M"  - number of taylor terms to sum
% Output Parameters
%   P - (cell array) nx1 cell array of dimensions where  P{i} = \phi_{i-1}(A)e

% Set number of Taylor Series Terms to sum
if(nargin == 3 && isfield(params,'M'))
    M = params.M;
else
    M = 20;
end

P = cell(n+1,1);  % Struct to store phi functions
s = max(0,ceil(log(max(norm(A,1),norm(A,Inf))/1)/log(2.d0)));
A = A/2^s;

% Initial Taylor Series Defenition
for m=0:n  
    % Horner's Scheme
    P{m+1} = A*(1/factorial(M+m));
    for k=(M-1):-1:1
        P{m+1} = A/factorial(k+m) + A*P{m+1};
    end
    P{m+1} = (1/factorial(m))*A^0 + P{m+1};
end

% === Alternative Algorithm for Taylor Series Defenition (Faster, may be less stable in certain cases)
% Ap = cell(M+1,1); % powers of A
% Ap{1} = A^0;
% for i=1:M+1
%     Ap{i+1} = A*Ap{i};
% end
% for m=0:n  
%     P{m+1} = Ap{1}/factorial(m);
%     for k=1:M
%         P{m+1} = P{m+1} + Ap{k+1}/factorial(m+k);
%     end
% end

% Scale up by Powers of Two
for i=1:s
    PS = P;
    for m=0:n
        if(m == 0)
            PS{1} =  PS{1}^2;
        else
            PS{m+1} = (1/2^m) * P{1} * P{m+1};
            for l=1:m
                 PS{m+1} = PS{m+1} + (1/2^m) * P{l+1} / factorial(m-l);
            end
        end
    end
    P = PS;
end
end