function [ww,wwsym]=RSC_heliCoupChain2A(vv,vvsym) %#ok<*AGROW>
%% Input: a standard JM.set-DF.expr chain (incl. DclF/DFcl etc.)
%% Output: a single ww sequence of chain decay wwsym amplitudes
[vm,vn]=size(vvsym);
for ii=1:1:vm
  for jj=1:1:vn
    if (vvsym(ii,jj)==-1000)
      vvsym(ii,jj)=0;
    end
  end
end
%% Standard JM.set-DF.expr chain multiply-by-step
for kk=1:1:vn
  vpick=vv(:,(7*kk-6):(7*kk-1));
  vspick=vvsym(:,kk);
  vsmat=[]; vpmat=[];
  for ii=1:1:vm
    if (vspick(ii)~=0)
      vsmat=[vsmat;vspick(ii)];
      vpmat=[vpmat;vpick(ii,:)];
    end
  end % v_mat series: step n
  vskm=length(vsmat);
  if (kk==1)
    vsmul=vsmat;
    vpmul=vpmat;
    continue
  end
  vsmull=[];
  vpmull=[];
  for ii=1:1:vskm
    for jj=1:1:length(vsmul)
      vsmull=[vsmull;vsmat(ii)*vsmul(jj)];
      vpmull=[vpmull;[vpmul(jj,:),-1000,vpmat(ii,:)]];
    end
  end % v_mull series: step n+1
  vsmul=vsmull;
  vpmul=vpmull; % v_mul series: step (n+1)-1
end
ww=vpmul;
wwsym=vsmul;
end