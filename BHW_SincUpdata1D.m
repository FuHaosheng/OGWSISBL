function [S0,Su_d,W] = BHW_SincUpdata1D(P,N,WindowsType)
% This function is used to generate the weighted window
% This program was last updated in 07-24-2023 , written by Haosheng Fu
% Input :           P - Length of window;
%                   N - Data length
%         WindowsType - Optional window type
% Output:          S0 - Zero-order approximation of the window function corresponds to I in (9) in the article 
%                Su_d - First-order approximation of the window function corresponds to S_dash in (10) in the article 
%                   W - Weighted window corresponds to W in (9) in the article
%-----------------------------------------------------------------------------------------------------------------------
% Blackman-Harris Window parameters
w0 = 0.42323; w1 = 0.49755; w2 = 0.07922; 
%% Generated weighted window
S0 = eye(P);  n = 0:1:N-1;  w = zeros(P,1);
switch WindowsType
    case 'RectWindow'
        w(1:N,1) = ones(N,1);
    case 'Three-termBHW'
        w(1:N,1) = w0*ones(1,N) + w1*cos(2*pi*(2*n)/N) + w2*cos(2*pi*(4*n)/N);
        w = w.';
    case 'blackmanharris'
        w(1:N,1) = blackmanharris(N);%'symmetric' or 'periodic'
    case 'barthannwin'
        w(1:N,1) = barthannwin(N);
    case 'blackman'
        w(1:N,1) = blackman(N);%'symmetric' or 'periodic'
    case 'bohmanwin'
        w(1:N,1) = bohmanwin(N);
    case 'chebwin'
        w(1:N,1) = chebwin(N,15);
    case 'flattopwin'
        w(1:N,1) = flattopwin(N);%'symmetric' or 'periodic'
    case 'gausswin'
        w(1:N,1) = gausswin(N,0.5);
    case 'hamming'
        w(1:N,1) = hamming(N);%'symmetric' or 'periodic'
    case 'hann'
        w(1:N,1) = hann(N);%'symmetric' or 'periodic'
    case 'kaiser'
        w(1:N,1) = kaiser(N,1);
    case 'nuttallwin'
        w(1:N,1) = nuttallwin(N);%'symmetric' or 'periodic'
    case 'parzenwin'
        w(1:N,1) = parzenwin(N);
    case 'taylorwin'
        w(1:N,1) = taylorwin(N,6,-20);
    case 'tukeywin'
        w(1:N,1) = tukeywin(N,0.5); 
end
n_ux = 0:1:P-1;
p1_ux = zeros(1,P);  p2_ux = zeros(1,P);
for i = 1:P
p1_ux(i) = w(i)*(-1)^(n_ux(i))*(1./n_ux(i));%w(i)*
p2_ux(i) = w(i)*(-1)^(-1*n_ux(i))*(1./(-1*n_ux(i)));
end
P1_ux = toeplitz(p1_ux);  P2_ux = toeplitz(p2_ux);
Su_d = triu(ones(P,P)).*P1_ux + tril(ones(P,P)).*P2_ux;
Su_d(logical(eye(size(Su_d)))) = zeros(P,1);%
W = toeplitz(w);

end



