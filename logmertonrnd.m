function Rsim=logmertonrnd(dt,mu,sigma,lambda,mu_j,sigma_j,t,Ns)
dN = poissrnd(lambda*dt, length(t)-1, Ns);  %poisson
Y = mu_j*dN + sigma_j*sqrt(dN).*randn(length(t)-1,Ns);  %jump
dW = sqrt(dt).*normrnd(0,1,length(t)-1, Ns);  %brownian
Rsim = (mu - sigma^2/2)*dt + sigma*dW + Y;
end
