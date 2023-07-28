function [v] = katz(A,t, rh_vector)
m = size(A,1);
n = size(A,2);

resolvent =speye(m,n) - t*A;
v = resolvent \ rh_vector;

end


