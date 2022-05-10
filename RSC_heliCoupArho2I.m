function [vv,vvsym] = RSC_heliCoupArho2I(ww,wwsym,rhomm) %#ok<*AGROW>
%% Input: a standard JM.set-DF.expr amplitude stat & density matrix
%% Output: a single JM.1-I.all intensity expr
if (isempty(ww))
  vv=ww; vvsym=wwsym;
else
[wm,wn]=size(wwsym);
rmin=min(ww(:,6));
rmax=max(ww(:,6));
vvsym=sym(zeros(1,wn));
for kk=1:1:wn % per sub-process (always 1)
  tempi=0;
  rnum=1;
  for rr=rmin:1:rmax % summing M
    snum=1;
    for ss=rmin:1:rmax % summing M'
      tempr=[];
      temps=[];
      for ii=1:1:wm
        if (ww(ii,6)==rr)
          tempr=[tempr,wwsym(ii,kk)];
        end
        if (ww(ii,6)==ss)
          temps=[temps,wwsym(ii,kk)];
        end % lam selected
      end % M, M' selected
      temprs=0;
      for jj=1:1:length(tempr)
        for hh=1:1:length(temps)
          temprs=temprs+rhomm(rnum,snum)*tempr(jj)*conj(temps(hh));
        end
      end % I term calculated
      tempi=tempi+temprs;
      snum=snum+1;
    end
    rnum=rnum+1;
  end
  vvsym(kk)=tempi;
end
vv=ww(1,:);
end
end