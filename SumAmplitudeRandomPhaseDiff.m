% expected amplitude for summing two sines with a random phase relationship
% whereby the difference in phases is normally distrobuted

% indefinite integral of this function wrt x gives amplitude for given sigma
Integrand = @(sigma, x) 1./(sigma.*sqrt(pi)).*exp(-0.5*(x./sigma).^2).*sqrt(1+cos(x));

% we'll loop over sigma values and calculate indefinite integral for each
sigmas = 0.1:0.1:10;
Integrals = zeros(size(sigmas));

phi = -100:0.01:100;    % range of integration

% for plotting integrands
%figure
%hold on

for index = 0:range(size(sigmas))
    currentIntegrand = Integrand(sigmas(index+1), phi);
    %plot(phi, currentIntegrand)
    
    currentIntegral = -trapz(currentIntegrand, phi);
    Integrals(index+1) = currentIntegral;
end

% plot ampltidue vs. sigma
figure
plot(sigmas, Integrals)
title('Expected Amplitude for Sum of 2 Sines with Gaussian Phase Difference Probability')
xlabel('Std. Dev. of Phase Difference Gaussian')
ylabel('Amplitude')
grid on

% fit resulting curve to normal distrobution
initGuess = [1,1,1];    % initial guess for amplitude, std dev, and +c
[fit, residual] = fminsearch(@(p) norm( p(1).*exp(-0.5.*((sigmas(:))./p(2)).^2) + p(3) -Integrals(:)), initGuess)
fit
residual

% plot normal fit
hold on
plot(sigmas, fit(1).*exp(-0.5.*((sigmas(:))./fit(2)).^2) + fit(3));


% logistic fit
initGuess = [0.7,2,2,1.2722];    % initial guess for amplitude, exponent, x0, and +c
[fit, residual,exitflag] = fminsearch(@(p) norm( p(1)./( 1 + exp(p(2)*(sigmas-p(3))) ) + p(4) -Integrals(:)), initGuess)
fit
residual
exitflag

% plot logistic fit
hold on
plot(sigmas, fit(1)./(1 + exp(fit(2)*(sigmas-fit(3)))) + fit(4));
%     ------> doesn't work for some reason
