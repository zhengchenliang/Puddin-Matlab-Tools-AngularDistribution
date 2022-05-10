function [ww,wwsym] = RSC_heliCoupAmpSkim(vv,vvsym)
%% Input: a standard JM.set-DF.expr amplitude (incl. DclF/DFcl etc.)
%% Output: a skimmed JM.set-DF.expr amplitude without zero
ww=vv; wwsym=vvsym;
vm=length(vvsym);
for ii=vm:-1:1
  if (vvsym(ii)==0)
    wwsym(ii)=[];
    ww(ii,:)=[];
  end
end
end

