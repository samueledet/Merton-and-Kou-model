function Rsim=logkournd(dt,mu,sigma,lambda,p,q,eta1,eta2,t,Ns)
dN = poissrnd(lambda*dt, length(t)-1, Ns);  %poisson
dW = sqrt(dt).*normrnd(0,1,length(t)-1, Ns);  %brownian
Y = (p/eta1 - q/eta2)*dN + 2*(p/eta1^2 + q/eta2^2)*dN.*randn(length(t)-1,Ns);
Rsim = (mu - sigma^2/2)*dt + sigma*dW + Y;
end
