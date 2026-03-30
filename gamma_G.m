function y = gamma_G(x,h)
    % delay function for bursting of red blood cells
    a1 = 24*33; % shape parameter
    b1 = 1/3; % scale parameter, mean should be alpha_G
    % mean = a1*b1
    %x = linspace(0,20*24,1000);
    % From Martcheva p305
    pi1 = 1-gamcdf(x,a1,b1); 
    pi2 = 1-gamcdf(x+h,a1,b1);
    %plot(x/24,y)

    y = 1/h*log(pi1./pi2);
    id = find(y>.6,1,'first');
    y(id:end) = .6;
end
