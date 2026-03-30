function y = betaHV(x)
% infectivity level based on gametocycte load
   y = 0.03*(x.^0.6)./(1+0.035*x.^0.6); % functional form taken from REF? need to update
end
