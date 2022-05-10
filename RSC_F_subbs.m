function oo=RSC_F_subbs(expr,n)
%% calculate subs(expr,1) for n times
if (fix(n)~=n || n<0)
  error('Input subbs time should be a natural number.')
elseif (n==0)
  oo=expr;
else
  for ii=1:1:n
    expr=subs(expr,1);
  end
  oo=expr;
end
end