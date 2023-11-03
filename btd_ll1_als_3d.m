function [A, B, C, const, output] = btd_ll1_als_3d(T, R, L, options, varargin)
% BTD_LL1_ALS_3D Computes the (Lr,Lr,1) Block-Term Decomposition (BTD) of a
% 3rd-order tensor using the alternative least squares algorithm (ALS).
%
% INPUT:
%   T (3D array): original 3-D tensor with dimensions I_1 x I_2 x I_3
%   R (integer): Number of components (sources) for the BTD decomposition
%   L (integer/vector): Rank of the factor matrices in the first and second mode 
%              (in your syntax, you can choose L to be equal for all 
%               components (int) or set it per component (vec))
%   options (struct, optional) : optimization options containing:
%          - th_relerr : relative error threshold
%          - maxiter   : max number of iterations
%   variable optional input (3-cell array OR 3 arrays, optional): initialization for 
%           the factor matrices
%
% OUTPUT:
%   A (2D matrix): mode-1 factor matrix with normalized columns (I_1 x (R x Lr))
%   B (2D matrix): mode-2 factor matrix with normalized columns (I_2 x (R x Lr))
%   C (2D matrix): mode-3 factor matrix with normalized columns (I_3 x (R x Lr))
%   const (vector):  vector containing the respective weights 
%   output (struct) : optimization options containing:
%          - relerr : relative error achieved, defined as Frobenius norm of 
%           the residual of the decomposition OVER Frobenius norm of the 
%           original tensor.
%          - numiter : the number of iterations the algorithm ran for
% 
% Number of iterations
maxiter = options.maxiter;
% Threshold for relative error
th_relerr = options.th_relerr;
% Check if the initialization for the factor matrices was given
init = [];
if ~isempty(varargin)
    if length(varargin) == 1    % Given as cell
        init = varargin{:}; 
    else                        % Given as matrices 
        init = varargin;
    end
end
% Initialize the three factor matrices 
if isempty(init)    % Randomly if not initialization was given
    A = randn(size(T, 1), R*L);
    B = randn(size(T, 2), R*L);
    C = randn(size(T, 3), R*1);
else                % Otherwise use the given initialization

    A = init{1};
    B = init{2};
    C = init{3};
end
% Normalize the columns of the initial factor matrices
A = normc(A);
B = normc(B);
C = normc(C);
T1 = hidden_mode_n_matricization(T, 1);
T2 = hidden_mode_n_matricization(T, 2);
T3 = hidden_mode_n_matricization(T, 3);
I1 = size(T,1);
I2 = size(T,2);
I3 = size(T,3);

for idxiter = 1:maxiter
   BC = [];
   B_new = reshape(B, I2, L, R);
    for dim=1:R
        BC = [BC, kron(C(:,dim),B_new(:,:,dim))];
    end
    A = T1*pinv(BC');
    A = normc(A);

    AC = [];
    A_new = reshape(A, I1, L, R);
    for dim=1:R
        AC = [AC, kron(C(:,dim),A_new(:,:,dim))];
    end
    B = T2*pinv(AC');
    B = normc(B);
    AB = [];
    B_new = reshape(B, I2, L, R);
    A_new = reshape(A, I1, L, R);
    for dim=1:R
        %AB = [AB, hidden_khatri_rao(B_new(:,:,dim),A_new(:,:,dim))];
        kr = hidden_khatri_rao(B_new(:,:,dim),A_new(:,:,dim));
        AB  = [AB kr* ones(L,1)];
    end
    C = T3*pinv(AB');
    const= vecnorm(C);
    C = normc(C);
    lamuda = diag(const);
    T3_est = C*lamuda*(AB)';
    
    % Calculate the relative error between the estimate and the true
    % mode-3 unfolding of the tensor
    relerr(idxiter) = norm((T3 - T3_est),"fro") / norm(T3,"fro");
    idxiter
    if relerr(idxiter) < th_relerr
        break;
    end
end
output.relerr = relerr;
output.numiter = idxiter;

% Warning if maximum number of iterations was reached
if idxiter == options.maxiter
    warning(['The ALS algorithm reached the maximum number of ' num2str(options.maxiter) ' iterations.'])
end
end
