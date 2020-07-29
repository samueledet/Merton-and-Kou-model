function pdflog = logmertonpdf(r,dt,mu,sigma,lambda,mu_j,sigma_j)
    if lambda > 0
        nterm = 100;
        term = zeros(length(r), nterm);  %matrix 5333*100
        for k = 1: nterm
            poisson = (lambda*dt)^k/prod(1:k)*exp(-lambda*dt);
            normal = 1/sqrt(2*pi*(sigma^2*dt + sigma_j^2*k))*exp(-(r-((mu-sigma^2/2)*dt + mu_j*k)).^2/(2*(sigma^2*dt + sigma_j^2*k)));
            term(:,k) = poisson*normal;
        end
        %pdflog = sum(term,2);
        pdflog = sum([1/sqrt(2*pi*sigma^2*dt)*exp(-(r-(mu-sigma^2/2)*dt).^2/(2*sigma^2*dt))*exp(-lambda*dt) term], 2);
    else
        pdflog = 1/sqrt(2*pi*sigma^2)*exp(-(r-(mu-sigma^2/2)*dt).^2/(2*(sigma^2*dt)));
    end
end

