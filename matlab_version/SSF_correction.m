function [S] = SSF_correction(f_v, theta, lamda, r)
    S=zeros(length(theta),1);
    for i1 = 1: length(theta)
        q = 4 * pi * sin(theta(i1)/2) / lamda;
        x = q * r;
        psi = sin(x) / x;
        sigma = 3*(sin(x)/x^3 - cos(x)/x^2);
        alpha = f_v/(1-f_v) * ((1 + 3*f_v/(1-f_v)) * sigma + 3*psi) + cos(x);
        beta = f_v/(1-f_v) * x * sigma + sin(x);
        S(i1)=1 / (alpha^2 + beta^2);
        % 
        % %different formula
        % % n = f_v *3/4 / pi /r^3;
        % 
        % p = 4*pi*sin(theta(i1)/2) / lamda;
        % 
        % alpha = (1 + 2*f_v)^2 / (1 - f_v)^4;
        % beta = -6*f_v*(1 + f_v/2)^2 / (1 - f_v)^4;
        % delta = alpha*f_v/2;
        % %delta = f_v*(1 + 2*f_v)^2 /(2*(1 - f_v)^2)
        % u=2*p*r;
        % 
        % if(abs(theta(i1)) <  10^-8)
        %     C = -24 * f_v * (alpha/3 + beta/4 + delta/6);
        % else
        %     C = 24 *(f_v) * ((alpha+beta+delta) * cos(u)/u^2 -(alpha+ 2*beta + 4*delta) * sin(u)/u^3 -  2*(beta+6*delta)*cos(u)/u^4 + 2*beta/u^4 + 24*delta*sin(u)/u^5 + 24*delta*(cos(u)-1)/u^6);
        % end
        % 
        % S(i1) = 1/(1-C);

    end
end