function [fitresult, gof, fout] = lsq2DGaussian(image, guess, display)
%Least Squares Gaussian fit

[sx,sy]=size(image);

[X,Y]=meshgrid(1:sy,1:sx);
[xData, yData, zData] = prepareSurfaceData( X,Y, image );

model ='A+N*(2*pi*sigma^2)^-1*exp(-((x-xo)^2+(y-yo)^2)*((2*sigma^2)^-1))';
ft = fittype( model, 'independent', {'x', 'y'}, 'dependent', 'z','coefficients',{'A','N','xo','yo','sigma'} ); %creates a symmetrical 2d gaussian fit
opts = fitoptions(ft); %Performing NonLinear Least Squares fitting, normalized

opts.Lower = guess.Lower;  %Setting upper and lower bounds, initial guess for fitting       
opts.StartPoint=guess.StartPoint;
opts.Upper=guess.Upper;

[fitresult, gof,fout] = fit( [xData, yData], zData,ft, opts);

if display == 1
    figure(2)
    % % Plot fit with data.
    plot( fitresult, [xData, yData], zData);
    % Label axes
    xlabel( 'X' );
    ylabel( 'Y' );
    zlabel( 'im' );
    grid on
    colormap(jet);
end
