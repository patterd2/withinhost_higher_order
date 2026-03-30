function y = phi(x,IT,s)

phi0 = 0.5*(1 + tanh((0-IT)/s));
y = (0.5*(1 + tanh((x-IT)/s)) - phi0)/(1 - phi0);

end
