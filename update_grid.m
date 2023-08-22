function data = update_grid(H,T,Y,Su,mu,Sigma,data,r,iter_beta)

BHB = (Su)' * (Su);  
G = (mu * mu') + T * Sigma;

P1 = real((conj(BHB) .* G) );%conjÇó¸´¹²éî

v1 = zeros(size(H,2),1);
for t = 1:T
    v1 = v1 + real(conj(mu(:,t)) .* ((Su)' * (Y(:,t) - H * mu(:,t))));%real(conj(mu(:,t)) .* (Su' * (Y(:,t) - H * mu(:,t) - Sv*diag(mu(:,t))*PHI2*g(:,2)))) ;
end
v1 = v1 - T * real(diag(Su' * H * Sigma));

temp1 =  P1 \ v1;
% for ii = 1:length(temp1)
%     if isnan(temp1(ii,1)) == 1
%         temp1(ii,1) = 0;
%     end
% end

if any(abs(temp1) > r/2) || any(diag(P1) == 0)
   for i = 1:iter_beta
       for n = 1:length(temp1)
           temp_beta = data;
           temp_beta(n) = 0;
           data(n) = (v1(n) -P1(n,:) * temp_beta) / P1(n,n);
           if data(n) > r/2
               data(n) = r/2;%*randn(1)
           end
           if data(n) < -r/2
               data(n) = -r/2;
           end
           if P1(n,n) == 0
               data(n) = 0;
           end
       end
   end
else
   data = temp1;
end
%                
% temp1 = data(idx);
end

