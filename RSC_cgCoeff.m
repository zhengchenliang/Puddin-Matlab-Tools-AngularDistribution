function ex=RSC_cgCoeff(pp,num2sym)
%% l1, l2, l; m1, m2, m; num2sym for sym()
l1=pp(1); l2=pp(2); l=pp(3);
m1=pp(4); m2=pp(5); m=pp(6);
%%
if (m==m1+m2)
  delta=1;
else
  delta=0;
end
sup=(2*l+1)*factori(l1+l2-l)*factori(l1-l2+l)*factori(-l1+l2+l);
sdn=factori(l1+l2+l+1);
ss=sqrt(sup/sdn);
dd=factori(l1+m1)*factori(l1-m1)*factori(l2+m2)*...
    factori(l2-m2)*factori(l+m)*factori(l-m);
dd=sqrt(dd);
mm=0;
for zz=0:1:(2*max([l1,l2,l,m1,m2,m]))
  if (l1+l2-l-zz<0 || l1-m1-zz<0 || l2+m2-zz<0 || l-l2+m1+zz<0 || l-l1-m2+zz<0)
    continue
  end
  mm=mm+((-1)^zz)/(factori(l-l2+m1+zz)*factori(l-l1-m2+zz)*...
      factori(zz)*factori(l1+l2-l-zz)*factori(l1-m1-zz)*factori(l2+m2-zz));
end
ex=delta*ss*dd*mm;
if (num2sym==0)
  ex=sym(ex);
end
end
function oo=factori(xx)
  if (xx<=0 || fix(xx)~=xx)
    oo=gamma(xx+1);
  else
    oo=factorial(xx);
  end
end