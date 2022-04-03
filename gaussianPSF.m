%2D gaussian function, also calulates radius @ FWHM:
function [PSF,max] = gaussianPSF(x0,y0,sigma) 
    syms x y
    AMP = (2*pi*sigma^2)^-1; % Normalized Scaling Factor
    PSF = AMP*exp(-((x-x0)^2+(y-y0)^2)*(2*sigma^2)^-1); 
    max = subs(PSF, [x, y], [x0,y0]);    
end
