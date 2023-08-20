load('hw13.mat');
close all;
x = x2;
figure
plot(x)

L = 100;      
K = 5;                                       
T = length(x);

alpha = ones(K,1);
tau = ((0:K-1)*floor((T-1)/K))'+1;

idx_z = (0:(L-1)).'+(1:(T-L+1));
Z = x(idx_z);



itrMAX=10;
obj_fun = zeros(1,itrMAX);
%% 
for itr=1:10   
    
    idx_y = tau'+(0:(L-1)).';
    Y = x(idx_y);
    S_hat = Y*(alpha);
    S_hat = S_hat/norm(S_hat,'fro');

    rho = S_hat.'*Z;
    for k = 1:K
        [alpha(k),tau(k)] = max(abs(rho));
        rho(max(tau(k)-L+1,1):min(T-L+1,tau(k)+L-1)) = 0;
    end
    
    x_hat = zeros(1,T);
    idx_x = tau'+(0:(L-1)).';
    x_hat(idx_x) = S_hat*alpha';
    
    obj_fun(itr) = norm(x - x_hat , 'fro')/norm(x,'fro');
end

[tau,I] = sort(tau);
alpha = alpha(I);


figure
subplot(3,1,1)
plot(S_hat)

subplot(3,1,2)
stem(tau , alpha ,'^')
xlim([0 T]);

subplot(3,1,3)
plot(x);
hold on
plot(x_hat,'r')
