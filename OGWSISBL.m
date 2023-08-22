function [x_est,Azimuth] = OGWSISBL( paras )
% This procedure is used to perform off-grid DOA estimation of the array
% This program was last updated in 07-24-2023 , written by Haosheng Fu
% From "Off-Grid Error Calibration for DOA Estimation Based on Sparse Bayesian Learning with Weighted Sinc Interpolation"
% IEEE Transactions on Vehicular Technology
% doi:10.1109/TVT.2023.3298965
%--------------------------------------------------------------------------
% Input :           Y - Measuremented Signal;
%                   r - Grid resolution
%                   p - Array element positions
%                 tol - Iterative tolerance (User set)
%             maxiter - Maximum number of iterations (User set)
% Output:       x_est - Signal reconstruction results 
%             Azimuth - Grid updata results (Eliminate off-grid errors) 
%--------------------------------------------------------------------------
Y = paras.Y; r = paras.resolution; p = paras.position;
tol = paras.tol; maxiter = paras.maxiter;
%--------------------------------------------------------------------------
search_area=[ -90 : r : 90 ];  GridNum = length(search_area);
[~, du, W] = BHW_SincUpdata1D(GridNum,GridNum,'kaiser');
a_search = search_area * pi / 180.;
H = exp( -1i * pi * p' * sin( a_search ) );
grid_u = search_area;   temp_u = grid_u.';
I = W .* eye( size( W, 1 ) );  A = H * I;
[ M, T ]=size(Y);  N = size(A,2);
%% Parameter Initial
a = 1e-4;  b = 1e-4;  rho = 0.01;
sigma2 = 10^(-2)*(norm(Y))^2/(M*T);
alpha0 = 1/sigma2;
alpha = sum(abs(A'*Y),2)/(M*T);
data = zeros(N,1); 
converged = false;
iter_beta = 1;
iter = 0;
alpha0seq = zeros(maxiter,2);
B = H*du;
%% off-grid DOA estimation
while ~converged
    iter = iter + 1;
    alpha_last = alpha;
    
    Phi = A + B * diag(data);
   
    C = 1 / alpha0 * eye(M) + Phi * diag(alpha) * Phi';
    Cinv = inv(C);
    Sigma = diag(alpha) - diag(alpha) * Phi' * Cinv * Phi * diag(alpha);
    mu = alpha0 * Sigma * Phi' * Y;
    
    

    for ll = 1:size(mu,1)
        mu_norm(ll,1) = norm(mu(ll,:));
    end

    gamma1 = 1 - real(diag(Sigma)) ./ (alpha);

    % update alpha
    musq = sum( mu.*conj(mu), 2) + T*real(diag(Sigma));%mean(abs(mu).^2, 2);
    alpha = ( -T+ sqrt(  T^2 + 4*rho* real(musq) ) ) / (  2*rho   );%musq + real(diag(Sigma));
 
    % update alpha0
    resid = Y - Phi * mu;
    alpha0 = (T * M + a - 1) / (norm(resid, 'fro')^2 + T / alpha0 * sum(gamma1) + b);
    alpha0seq(iter) = alpha0;

    % stopping criteria
    err = norm(alpha - alpha_last)/norm(alpha_last);
    if err < tol || iter >= maxiter
        converged = true;
        iter_beta = 5;
    end
    % update grids
    data_last = data;
    data = update_grid(A,T,Y,B,mu,Sigma,data_last,1,iter_beta);
    temp_u = grid_u.' + r*data;

end
x_est = abs(mu_norm).^2 + abs(diag(Sigma));
Azimuth = temp_u;
end




