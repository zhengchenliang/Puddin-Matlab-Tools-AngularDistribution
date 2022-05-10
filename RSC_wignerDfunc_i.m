function ex=RSC_wignerDfunc_i(pp)
%% a series for ^__ & b series for alpha, beta, gamma.
%% _i for partial fast calc
a1=sym(pp(1)); a2=sym(pp(2)); a3=sym(pp(3));
b1=sym(pp(4)); b2=sym(pp(5)); b3=sym(pp(6));
%%
if (a1==0)
  ex=1;
end
if (a1==1)
  if (a2==1)
    if (a3==1)
      ex=(1+cos(b2))*(exp(-1i*(b1+b3)))/2;
    elseif (a3==0)
      ex=-(sin(b2)*exp(-1i*b1))/sqrt(2);
    elseif (a3==-1)
      ex=(1-cos(b2))*(exp(-1i*(b1-b3)))/2;
    else
    end
  elseif (a2==0)
    if (a3==1)
      ex=(sin(b2)*exp(-1i*b3))/sqrt(2);
    elseif (a3==0)
      ex=cos(b2);
    elseif (a3==-1)
      ex=-(sin(b2)*exp(1i*b3))/sqrt(2);
    else
    end
  elseif (a2==-1)
    if (a3==1)
      ex=(1-cos(b2))*(exp(1i*(b1-b3)))/2;
    elseif (a3==0)
      ex=(sin(b2)*exp(1i*b1))/sqrt(2);
    elseif (a3==-1)
      ex=(1+cos(b2))*(exp(1i*(b1+b3)))/2;
    else
    end
  else
  end
elseif (a1==2)
  if (a2==2)
    if (a3==2)
      ex=(exp(-2*1i*(b1+b3)))*((1+cos(b2))/2)^2;
    elseif (a3==1)
      ex=-(exp(-1i*(2*b1+b3)))*(1+cos(b2))*(sin(b2))/2;
    elseif (a3==0)
      ex=(exp(-2*1i*b1))*(sin(b2)^2)*sqrt(3/8);
    elseif (a3==-1)
      ex=-(exp(-1i*(-2*b1+b3)))*(1-cos(b2))*(sin(b2))/2;
    elseif (a3==-2)
      ex=(exp(2*1i*(-b1+b3)))*((1-cos(b2))/2)^2;
    else
    end
  elseif (a2==1)
    if (a3==2)
      ex=(exp(-1i*(b1+2*b3)))*(1+cos(b2))*(sin(b2))/2;
    elseif (a3==1)
      ex=(exp(-1i*(b1+b3)))*((cos(b2)^2)-(1-cos(b2))/2);
    elseif (a3==0)
      ex=-(exp(-1i*b1))*(sin(2*b2))*sqrt(3/8);
    elseif (a3==-1)
      ex=(exp(1i*(-b1+b3)))*(-(cos(b2)^2)+(1+cos(b2))/2);
    elseif (a3==-2)
      ex=-(exp(1i*(-b1+2*b3)))*(1-cos(b2))*(sin(b2))/2;
    else
    end
  elseif (a2==0)
    if (a3==2)
      ex=(exp(-2*1i*b3))*(sin(b2)^2)*sqrt(3/8);
    elseif (a3==1)
      ex=(exp(-1i*b3))*(sin(2*b2))*sqrt(3/8);
    elseif (a3==0)
      ex=(3*(cos(b2)^2)-1)/2;
    elseif (a3==-1)
      ex=-(exp(1i*b3))*(sin(2*b2))*sqrt(3/8);
    elseif (a3==-2)
      ex=(exp(2*1i*b3))*(sin(b2)^2)*sqrt(3/8);
    else
    end
  elseif (a2==-1)
    if (a3==2)
      ex=(exp(1i*(b1-2*b3)))*(1-cos(b2))*(sin(b2))/2;
    elseif (a3==1)
      ex=(exp(1i*(b1-b3)))*(-(cos(b2)^2)+(1+cos(b2))/2);
    elseif (a3==0)
      ex=(exp(1i*b1))*(sin(2*b2))*sqrt(3/8);
    elseif (a3==-1)
      ex=(exp(1i*(b1+b3)))*((cos(b2)^2)-(1-cos(b2))/2);
    elseif (a3==-2)
      ex=-(exp(1i*(b1+2*b3)))*(1+cos(b2))*(sin(b2))/2;
    else
    end
  elseif (a2==-2)
    if (a3==2)
      ex=(exp(2*1i*(b1-b3)))*((1-cos(b2))/2)^2;
    elseif (a3==1)
      ex=(exp(1i*(2*b1-b3)))*(1-cos(b2))*(sin(b2))/2;
    elseif (a3==0)
      ex=(exp(2*1i*b1))*(sin(b2)^2)*sqrt(3/8);
    elseif (a3==-1)
      ex=(exp(1i*(2*b1+b3)))*(1+cos(b2))*(sin(b2))/2;
    elseif (a3==-2)
      ex=(exp(2*1i*(b1+b3)))*((1+cos(b2))/2)^2;
    else
    end
  else
  end
else
end
end