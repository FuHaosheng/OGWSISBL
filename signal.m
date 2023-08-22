function [Y] = signal(p, theta, SNR, K, Mode)
% This procedure is used to generate the measurement signal Y for simulation
% This program was last updated in 07-24-2023 , written by Haosheng Fu
%--------------------------------------------------------------------------
% Input :    position - Array element positions;
%               theta - Angle grid
%                 SNR - Signal to Noise Ratio (SNR)
%                 K   - Number of snapshots
%                Mode - Source type Settings (Coherent or Incoherent)
% Output:       Y     - Measurement Signal
%--------------------------------------------------------------------------
M = length(p);  N_alpha=length(theta);
A=exp(-1i*pi*p'*sin(theta*pi/180));
switch Mode
%     case 'EqualitySource'
%         S=(randn(N_alpha,K)+1j*randn(N_alpha,K));  % complex signals
%         Vj=diag(sqrt(  10^(SNR/10)./diag(1/K*(S*S') ) ) );
%         s=Vj*S;
%         S=[s(1,:);s(1,:);s(1,:)];       
    case 'CoherentSource'
%         alpha1 = -15 + (10)*rand(1,1);
%         alpha2 = 10 + (10)*rand(1,1); 
%         TrueDOAs = round([alpha1 alpha2]*100)/100;% alpha3
        S = ( randn( N_alpha, K ) + 1i * randn( N_alpha, K ) );  % complex signals
        Vj = diag(sqrt(  10^(SNR/10)./diag(1/K*(S*S') ) ) );
        s = Vj * S;
        S = [ (0.5 + 1i*0.8) * s(1, :); (-0.7-1i*0.7) * s(1, :) ];               %пе╨е
    case 'UncoherentSource'
%         alpha1 = -35 + (10)*rand(1,1);
%         alpha2 = -5 + (10)*rand(1,1); alpha3 = 20 + (10)*rand(1,1);
%         TrueDOAs = round([alpha1 alpha2 alpha3]*100)/100;% alpha3
        S = randn( N_alpha, K ) + 1i * randn(N_alpha, K);  % complex signals
        Vj = diag( sqrt(  10^( SNR/10 ) ./ diag( 1/K * (S * S') ) ) );
        S = Vj * S;
end

noise = sqrt(1/2) * ( randn(M, K) + 1i * randn(M, K) );
Y = A * S + noise;

end

