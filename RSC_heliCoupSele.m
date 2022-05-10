function [qq,qcoeff]=RSC_heliCoupSele(pp,num2sym) %#ok<*AGROW>
%% j1, j2, j; p1, p2, p; c1, c2, c; num2sym for sym()
%% qq row for [j,lam1,lam2,LL,SS]; qcoeff row for F2G coefficients
j1=pp(1); j2=pp(2); j=pp(3);
p1=pp(4); p2=pp(5); p=pp(6);
c1=pp(7); c2=pp(8); c=pp(9);
%%
if (abs(c1)~=1 || abs(c2)~=1 || abs(c)~=1)
  error('charge configuration incorrect.')
end
if (abs(p1)~=1 || abs(p2)~=1 || abs(p)~=1)
  error('Parity configuration incorrect.')
end
if (j1<0 || j2<0 || j<0 || rem(2*rem(j1,1),1)~=0 || ...
    rem(2*rem(j2,1),1)~=0 || rem(2*rem(j,1),1)~=0)
  error('Total angular momentum unphysical.')
end
%%
qq=[]; qcoeff=[];
for lam1=-j1:1:j1
  for lam2=-j2:1:j2
    del=lam1-lam2;
    for LL=(rem(j,1)):1:j
      nfactor=sqrt((2*LL+1)/(2*j+1));
      for SS=(rem(j,1)):1:j
        if (abs(LL-SS)>j || LL+SS<j || abs(j1-j2)>SS || j1+j2<SS)
          continue
        end
        cgt1=RSC_cgCoeff([LL,SS,j,0,del,del],1);
        cgt2=RSC_cgCoeff([j1,j2,SS,lam1,-lam2,del],1);
        cgs1=RSC_cgCoeff([LL,SS,j,0,del,del],0);
        cgs2=RSC_cgCoeff([j1,j2,SS,lam1,-lam2,del],0);
        ncgt=nfactor*cgt1*cgt2;
        ncgs=nfactor*cgs1*cgs2;
        if (ncgt~=0)
          qq=[qq;[j,lam1,lam2,LL,SS]];
          if (num2sym==0)
            qcoeff=[qcoeff;ncgs];
          else
            qcoeff=[qcoeff;ncgt];
          end
        end
      end
    end
  end
end
for ii=length(qcoeff):-1:1
  if (isnan(qcoeff(ii)))
    qcoeff(ii)=[];
    qq(ii,:)=[];
  end
end
end