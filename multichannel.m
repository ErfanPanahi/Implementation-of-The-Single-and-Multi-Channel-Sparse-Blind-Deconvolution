load('hw13.mat');
close all;

figure
subplot(2,1,1)
plot(X(1,:))
subplot(2,1,2)
plot(X(2,:))



L = 100;      
K = 5;                                        
[N,T] = size(X);

A = randn(N,N);
A = A./repmat(vecnorm(A),N,1);
S_hat = zeros(L,N);
alpha = zeros(K,N);
tau = zeros(K,N);

itrMAX=10;
obj_fun=zeros(1,itrMAX);
%%
for itr=1:itrMAX
    Z = pinv(A)*X;
    
    for n = 1:N
        [S_hat(:,n),alpha(:,n),tau(:,n),er] = SBD(Z(n,:) , L , K);
        zz = zeros(1,T);
        idx = tau(:,n)'+(0:(L-1)).';
        zz(idx) = S_hat(:,n)*alpha(:,n)';
        Z(n,:)=zz;
    end

    A = X*pinv(Z);
    A = A./repmat(vecnorm(A),N,1);

    X_hat = A*Z;
    obj_fun(itr) = norm(X - X_hat , 'fro')/norm(X,'fro');

end

%%
X1=A(:,1)*Z(1,:);
X2=A(:,2)*Z(2,:);


close all
figure
subplot(4,1,1)
plot(S_hat(:,1))
xlim([0 T])
subplot(4,1,2)
stem(tau(:,1) , alpha(:,1) ,'^')
xlim([0 T])
subplot(4,1,3)
plot(1:T,X(1,:),'LineWidth',4)
hold on
plot(1:T,X1(1,:),'r','LineWidth',2)
subplot(4,1,4)
plot(1:T,X(2,:),'LineWidth',4)
hold on
plot(1:T,X1(2,:),'r','LineWidth',2)

figure
subplot(4,1,1)
plot(S_hat(:,2))
xlim([0 T])
subplot(4,1,2)
stem(tau(:,2) , alpha(:,2) ,'^')
xlim([0 T])
subplot(4,1,3)
plot(1:T,X(1,:),'LineWidth',4)
hold on
plot(1:T,X2(1,:),'r','LineWidth',2)
subplot(4,1,4)
plot(1:T,X(2,:),'LineWidth',4)
hold on
plot(1:T,X2(2,:),'r','LineWidth',2)


figure
subplot(2,1,1)
plot(X(1,:),'LineWidth',4);hold on;plot(X_hat(1,:),'r','LineWidth',2)
subplot(2,1,2)
plot(X(2,:),'LineWidth',4);hold on;plot(X_hat(2,:),'r','LineWidth',2)