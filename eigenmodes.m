% Following Stone (1970)

clr

% more points in k
% numerical parameters
tolerance = 1e-4;
eps = 1e-3;
iterlimit = 100;

sigma0 = -1i;
k = 01;[0.2:0.05:2];
l = 15;[0:10];
Ri = 0.5;

z = linspace(0,1,100)';
dz = z(2)-z(1);
sigma_c = nan(length(k),length(l),length(Ri));
sigma_found = sigma0;

for rr = 1:length(Ri)
    for ll = 1:length(l)
        for kk = 1:length(k)
        
            % initialize
            phi = nan([length(z) 1]);
            phi(1) = 0;
            phi(2) = 1; % set amplitude
            phi(end) = 1;  
            
            iter = 1; %count = 1;
            sigma = sigma_found;
            while abs(phi(end)) >= tolerance
                h = sigma + k(kk)*z; 
                a = (k(kk)^2*(Ri(rr)*h.^4 + (3-Ri(rr))*h.^2 -2) + l(ll)^2 * (Ri(rr)*h.^4 + (1-Ri(rr))*h.^2))./(h.^2 .* (1-h.^2).^2);
  
                % solve the equation
                for n=2:length(z)-1
%                     h = sigma + k(kk)*z(n);
%                     np1 = (1-h.^2)/(2*dz^2) - 1/dz * (k(kk)./h - 1i*l);
%                     nn  = -1 * ( (1-h.^2)/dz^2 + Ri*(k(kk)^2 + l^2) + 2*1i*k(kk)*l./h);
%                     nm1 = (1-h.^2)./(2*dz^2) + 1/dz * (k(kk)./h - 1i*l);
%                     phi(n+1) = - (nn*phi(n) + nm1 * phi(n-1))./ (np1);
                    phi(n+1) = (2 - dz^2*a(n))*phi(n) - phi(n-1);
                end

                % check with lower BC and change sigma if needed
                % shooting method - (sigma1,phi1) - (sigma2,phi2) linear fit and extrapolate to find phi2=0
                if abs(phi(end) - 0) >= tolerance
                    if iter>= 2 && phi(end) ~= philast
                        slope = ((sigma-sigmalast)./(phi(end)-philast));
                        sigmalast = sigma;
                        sigma = sigma - slope * phi(end);
                        %count = 1;
                    else
                        sigmalast = sigma;
                        sigma = sigma - eps;
                        %count=count+1;
                    end
                end
                iter = iter+1;
                philast = phi(end);

                %[k l iter phi(end) sigma]

                figure(1);
                plot(phi);
                phisave(kk,ll,rr,:) = phi;
                sigma_c(kk,ll,rr) = sigma;
                if iter >=iterlimit %|| max(phi) > 1e10
                    phisave(kk,ll,rr,:) = nan(size(phi));
                    sigma_c(kk,ll,rr) = NaN;
                    sigma_found = sigma0;
                    warning('Reached max. no. of iterations. Terminating.'); 
                    break;
                else
                    sigma_found = sigma;
                end
            end           
            %pause;
           if isnan(sigma_c(kk,ll,rr))
               warning('bad value found');
               break;
           else
               fprintf('\n Found sigma_c = %.3f, iteration = %d', -1*imag(sigma_c(kk,ll,rr)), iter);
           end
        end
    end
end
sigma_c = -addnan(imag(sigma_c),1);

fprintf('\n');
%%

if max(sigma_c > 1), sigma_c = addnan(sigma_c,1); end

colors = distinguishable_colors(length(l));

for ii=1:length(l)
    plot(k,sigma_c(:,ii),'Color',colors(ii,:));
    hold on
    leg{ii} = num2str(l(ii));
end
legend(cellstr(leg));
% w = real(h./sqrt(1-h) .* ((1-h)./(1+h)).^(1i*l/2/k) .* phi);