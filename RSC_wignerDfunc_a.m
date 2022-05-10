function ex=RSC_wignerDfunc_a(pp)
%% a series for ^__ & b series for alpha, beta, gamma.
%% _a for all standard calc
a1=sym(pp(1)); a2=sym(pp(2)); a3=sym(pp(3));
b1=sym(pp(4)); b2=sym(pp(5)); b3=sym(pp(6));
%%
ep1=exp(-1i*a2*b1);
ep2=exp(-1i*a3*b3);
ss=factorial(a1+a2)*factorial(a1-a2)*factorial(a1+a3)*factorial(a1-a3);
ss=sqrt(ss);
mm=0;
for zz=0:1:(2*max([abs(a1),abs(a2),abs(a3)]))
  if (a1+a3-zz<0 || a2-a3+zz<0 || a1-a2-zz<0)
    continue
  end
  mm=mm+((-1)^(a2-a3+zz))*((cos(b2/2))^(2*a1+a3-a2-2*zz))*...
      ((sin(b2/2))^(a2-a3+2*zz))/...
      (factorial(a1+a3-zz)*factorial(zz)*factorial(a2-a3+zz)*factorial(a1-a2-zz));
end
ex=ep1*ep2*ss*mm;
ex=simplify(ex);
end