function [tt,ttsym]=RSC_heliCoupDFcj(pp,sX) %#ok<*AGROW,*ASGLU>
%% j1, j2, j; p1, p2, p; c1, c2, c; sX for alpha, beta, gamma; num2sym for sym()
%% tt row for [j1,j2,j,lam1,lam2,M]; ttsym row for DF term
%% only two-body decay once, parity sym considered.
%% Pend: iii) A2I DM; vi) SR;
j1=pp(1); j2=pp(2); j=pp(3);
p1=pp(4); p2=pp(5); p=pp(6);
c1=pp(7); c2=pp(8); c=pp(9); 
x1=sX(1); x2=sX(2); % x3=sX(3);
%%
[ss,sssym]=RSC_heliCoupD(pp,sX); [sz,~]=size(ss);
tt=ss; ttsym=[];
for ii=1:1:sz
  spick=ss(ii,:);
  eval(['syms FJ',num2str(spick(3)),'L',sname(spick(4)),...
      'N',sname(spick(5)),'c;'])
  eval(['ttsym=[ttsym; sssym(ii)*FJ',num2str(spick(3)),'L',sname(spick(4)),...
      'N',sname(spick(5)),'c];'])
end
end
function oo=sname(ii)
if (ii<0)
  oo=['m',num2str(abs(ii))];
else
  oo=['p',num2str(ii)];
end
end