function pdflog = logkoupdf(r,dt,mu,sigma,lambda,p,q,eta1,eta2)
    if lambda > 0
        f1 = ((1-lambda*dt)/(sigma*sqrt(dt)))*(1/sqrt(2*pi*sigma^2*dt))*exp(-(r-(mu-sigma^2/2)*dt).^2 /(2*sigma^2*dt));
        up = p*eta1*exp(0.5*sigma^2*eta1^2*dt)*exp(-(r-(mu-sigma^2/2)*dt)*eta1).*normalize((r-(mu-sigma^2/2)*dt - (eta1*sigma^2*dt))/(sigma*sqrt(dt)));
        down = q*eta2*exp(0.5*sigma^2*eta2^2*dt)*exp(-(r-(mu-sigma^2/2)*dt)*eta2).*normalize(-(r-(mu-sigma^2/2)*dt + (eta2*sigma^2*dt))/(sigma*sqrt(dt)));
        %up = p*eta1*exp(0.5*sigma^2*eta1^2*dt)*exp(-(r-(mu-sigma^2/2)*dt)*eta1)*int((r-(mu-sigma^2/2)*dt - (eta1*sigma^2*dt))/(sigma*sqrt(dt)),[0,100000]);
        %down = q*eta2*exp(0.5*sigma^2*eta2^2*dt)*exp(-(r-(mu-sigma^2/2)*dt)*eta2)*int(-(r-(mu-sigma^2/2)*dt + (eta2*sigma^2*dt))/(sigma*sqrt(dt)),[-100000,0]);
        
        pdflog = f1 + (lambda*dt*(up+down));
    else
        pdflog = 1/sqrt(2*pi*sigma^2)*exp(-(r-(mu-sigma^2/2)*dt).^2/(2*(sigma^2*dt)));
    end
        
end

