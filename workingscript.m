%Merton model
S= csvread('japan_stockprice.csv');
dt = 1/252;
t = linspace(0,(length(S)-1)*dt,length(S));
R = diff(log(S),1);
epsilon = 0.02;
jumpindex = find(abs(R)>epsilon);
lambdahat=length(jumpindex)/((length(S)-1)*dt);
Rjumps = R(jumpindex);
diffusionindex = find(abs(R)<=epsilon);
Rdiffusion = R(diffusionindex);
sigmahat = std(Rdiffusion)/sqrt(dt);
muhat = (2*mean(Rdiffusion)+(sigmahat^2)*dt)/(2*dt);
sigma_jhat = sqrt((var(Rjumps)-sigmahat^2*dt));
mu_jhat = mean(Rjumps)-(muhat-sigmahat^2/2)*dt;

%the merton model
theta0 = [muhat sigmahat lambdahat mu_jhat sigma_jhat];
Logmerton=@(mu, sigma, lambda, mu_j, sigma_j)-sum(log(logmertonpdf(R,dt,mu,sigma,lambda,mu_j,sigma_j)));
options = optimset('MaxFunEvals',10000);
[theta,fval] = fminsearch(@(theta)Logmerton(theta(1), theta(2), theta(3), theta(4), theta(5)), theta0, options);
Rsim = logmertonrnd(dt, theta(1), theta(2), theta(3), theta(4), theta(5), t,1);
Logmertons = @(mu,sigma,lambda,mu_j,sigma_j)-sum(log(logmertonpdf(Rsim,dt,mu,sigma,lambda,mu_j,sigma_j)));
thetas = fminsearch(@(theta)Logmertons(theta(1), theta(2), theta(3), theta(4), theta(5)), theta, options);
disp(['mu: ' num2str([theta0(1), theta(1), thetas(1)])])
disp(['sigma: ' num2str([theta0(2), theta(2), thetas(2)])])
disp(['lambda: ' num2str([theta0(3), theta(3), thetas(3)])])
disp(['muj: ' num2str([theta0(4), theta(4), thetas(4)])])
disp(['sigmaj: ' num2str([theta0(5), theta(5), thetas(5)])])

%Kolmogorov test
[h,p] = kstest2(R,Rsim)

%plot density function of empirical and simulated log returns
[f,xi] = ksdensity(R);
plot(xi,f,'-r','Linewidth',2);
hold on
[g,yi] = ksdensity(Rsim);
plot(yi,g,'-b','Linewidth',2)
legend('Empirical','Merton','Location','northwest')
title('Density function of empirical and simulated log return')
hold off


