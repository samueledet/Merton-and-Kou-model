%Kou model
S= csvread('japan_stockprice.csv');
dt = 1/252;
t = linspace(0,(length(S)-1)*dt,length(S));
R = diff(log(S),1);
epsilon = 0.02;
jumpupindex = find(R > epsilon);
jumpdownindex = find(R < -epsilon);
lambdahat=(length(jumpupindex)+length(jumpdownindex))/((length(S)-1)*dt);

%diffusion
diffusionindex = find(abs(R) <=epsilon);
Rdiffusion = R(diffusionindex);
sigmahat = std(Rdiffusion)/sqrt(dt);
muhat = (2*mean(Rdiffusion)+(sigmahat^2)*dt)/(2*dt);

%probabilities of upward and downward jumps
phat = length(jumpupindex)/(length(jumpupindex)+length(jumpdownindex));
qhat = length(jumpdownindex)/(length(jumpupindex)+length(jumpdownindex));

%average return of upward and downward jump
Rjumpup = R(jumpupindex);
Rjumpdown = R(jumpdownindex);
eta1hat = 1/mean(Rjumpup);
eta2hat  = 1/mean(Rjumpdown);

%Kou model
theta0 = [muhat sigmahat lambdahat phat qhat eta1hat eta2hat];
Logkou=@(mu, sigma, lambda, p, q, eta1, eta2)-sum(log(logkoupdf(R,dt,mu,sigma,lambda,p,q,eta1,eta2)));
options = optimset('MaxFunEvals',1000, 'MaxIter',1000);
[theta,fval] = fminsearch(@(theta)Logkou(theta(1), theta(2), theta(3), theta(4), theta(5), theta(6), theta(7)), theta0, options);

%Simulated log return
Rsim = logkournd(dt, theta(1), theta(2), theta(3), theta(4), theta(5), theta(6), theta(7), t,1);
Logkous = @(mu,sigma,lambda,p,q,eta1,eta2)-sum(log(logkoupdf(Rsim,dt,mu,sigma,lambda,p,q,eta1,eta2)));
thetas = fminsearch(@(theta)Logkous(theta(1), theta(2), theta(3), theta(4), theta(5), theta(6), theta(7)), theta, options);

%Display simulated estimated parameters
disp(['mu: ' num2str([theta0(1), theta(1), thetas(1)])])
disp(['sigma: ' num2str([theta0(2), theta(2), thetas(2)])])
disp(['lambda: ' num2str([theta0(3), theta(3), thetas(3)])])
disp(['p: ' num2str([theta0(4), theta(4), thetas(4)])])
disp(['q: ' num2str([theta0(5), theta(5), thetas(5)])])
disp(['eta1: ' num2str([theta0(6), theta(6), thetas(6)])])
disp(['eta2: ' num2str([theta0(7), theta(7), thetas(7)])])

%Kolmogorov test
[h,p] = kstest2(R,Rsim)

%plot density function of empirical and simulated log returns
[f,xi] = ksdensity(R);
plot(xi,f,'-r','Linewidth',2);
hold on
[g,yi] = ksdensity(Rsim);
plot(yi,g,'-b','Linewidth',2)
legend('Empirical','Kou','Location','northwest')
title('Density function of empirical and simulated log return')
hold off



