% expected amplitude for summing two sines with a random phase relationship
% whereby the difference in phases is normally distrobuted

% indefinite integral of this function wrt x gives amplitude for given sigma
Integrand = @(sigma, x) 1./(sigma.*sqrt(pi)).*exp(-0.5*(x./sigma).^2).*sqrt(1+cos(x));

% loop over sigma values and calculate indefinite integral for each
sigmas = 0.1:0.1:10;

Integrals = zeros(size(sigmas));

% range of integration
phi = -100:0.01:100;

% for plotting integrands
figure
hold on

for index = 0:range(size(sigmas))
    currentIntegrand = Integrand(sigmas(index+1), phi);
    plot(phi, currentIntegrand)
    
    currentIntegral = -trapz(currentIntegrand, phi);
    Integrals(index+1) = currentIntegral;
end

figure
plot(sigmas, Integrals)
title('Expected Amplitude for Sum of 2 Sines with Gaussian Phase Difference Probability')
xlabel('Std. Dev. of Phase Difference Gaussian')
ylabel('Amplitude')
grid on