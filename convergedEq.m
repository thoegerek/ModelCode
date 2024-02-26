function [conv,ideq] = convergedEq(eqs,st,x0)
%Finds the equilibrium that x0 converges to from a vector of equilibria eqs and their stability st
%eqs: 1-d vecotr of equilibria
%st: flags for stability of each equilibrium (-1 = stable) ,same dim as eqs
%x0: starting value of interest

unsorted = eqs;
[eqs,idx] = sort(eqs);
st = st(idx);

above = (eqs-x0)>0;
below = (eqs-x0)<0;

stabove = st(above);
if ~isempty(stabove) && stabove(1) == -1
    conv = eqs(above);
    conv = conv(1);
else
    conv = eqs(below);
    conv = conv(end);
end

ideq = find(unsorted == conv);
