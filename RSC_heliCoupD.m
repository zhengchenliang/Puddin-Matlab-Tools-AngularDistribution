function [ss,sssym]=RSC_heliCoupD(pp,sX) %#ok<*AGROW>
%% j1, j2, j; p1, p2, p; c1, c2, c; sX for alpha, beta, gamma; num2sym for sym()
%% ss row for [j1,j2,j,lam1,lam2,M]; sssym row for D-matrix element
%% only two-body decay once, parity sym considered.
%% Pend: iii) A2I DM; vi) SR;
j1=pp(1); j2=pp(2); j=pp(3);
p1=pp(4); p2=pp(5); p=pp(6);
c1=pp(7); c2=pp(8); c=pp(9); 
x1=sX(1); x2=sX(2); % x3=sX(3);
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
%% P conservation
rr=[]; rrsym=[];
[qq,~]=RSC_heliCoupSele(pp,0); [qs,~]=size(qq);
qmemo=[-1000,-1000,-1000];
qpx=p*p1*p2*(-1)^(j-j1-j2); % Parity conservation factor
for ii=1:1:qs
  qpick=qq(ii,:);
  qpicked=qpick(1:3);
  qpicked_sym=[qpicked(1),-qpicked(2),-qpicked(3)]; % Parity conservation
  veto=0; [qm,~]=size(qmemo);
  for jj=1:1:qm
    if (qpicked==qmemo(jj,:))
      veto=1;
    end
  end
  qmemo=[qmemo;qpicked];
  qmemo=[qmemo;qpicked_sym];
  if (veto==0)
    for MM=-j:1:j
      lam1=qpicked(2); lam2=qpicked(3); del=lam1-lam2;
      syms AA
      AA=RSC_wignerDfunc_a([j,MM,del,x1,x2,0])+...
         RSC_wignerDfunc_a([j,MM,-del,x1,x2,0])*qpx; % two-body only
      rr=[rr;[j1,j2,qpicked,MM]];
      rrsym=[rrsym;simplify(AA)];
    end
  end
end
%% CP conservation
if (j1==j2)
  ss=[]; sssym=[];
  [rs,~]=size(rr);
  rmemo=[-1000,-1000,-1000];
  for ii=1:1:rs
    rpick=rr(ii,:);
    rpicked=rpick(4:6);
    rpicked_sym=[rpicked(2),rpicked(1),rpicked(3)]; % CP conservation
    if (rpicked(1)~=rpicked(2))
      rmemo=[rmemo;rpicked_sym];
    end
  end
  [rm,~]=size(rmemo);
  veto=-1000*ones(rs,3);
  for ii=1:1:rs
    rpick=rr(ii,:);
    rpicked=rpick(4:6);
    for jj=1:1:rm
      if (rpicked==rmemo(jj,:))
        veto(ii,:)=[rpicked(1)+rpicked(2),abs(rpicked(1)-rpicked(2)),rpicked(3)];
      end
    end
  end
  rpx=c1*c2*(-1)^j; % CP conservation factor
  vetos=zeros(1,rs);
  for ii=1:1:rs-1
    rpickii=rr(ii,:);
    rsympickii=rrsym(ii,:);
    for jj=(ii+1):1:rs
      rsympickjj=rrsym(jj,:);
      if (veto(ii,1)==veto(jj,1) && veto(ii,2)==veto(jj,2) &&...
          veto(ii,3)==veto(jj,3) && sum(veto(ii,:))~=-3000)
        vetos(ii)=1;
        vetos(jj)=1;
        ss=[ss;rpickii];
        sssym=[sssym;rsympickii+rpx*rsympickjj];
      end
    end
    if (vetos(ii)==0)
      ss=[ss;rpickii];
      sssym=[sssym;rsympickii];
    end
  end
  rpickrs=rr(rs,:);
  rsympickrs=rrsym(rs,:);
  if (vetos(rs)==0)
    ss=[ss;rpickrs];
    sssym=[sssym;rsympickrs];
  end
  % sssym=conj(sssym); % the bra counts for the conj()
else
  ss=rr;
  sssym=rrsym;
  % sssym=conj(sssym); % the bra counts for the conj()
end
end