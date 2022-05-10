function [uu,uusym] = RSC_heliCoupDcjChain(varargin) %#ok<*AGROW>
%% Chain any D: pp1, sX1, pp2, sX2, ...
%% uu & uusym per 6+1 col for ss & sssym (D)
lv=length(varargin);
if (fix(lv/2)~=lv/2)
  error('Input chain incorrect.')
end
uu=[]; uusym=[];
temp=cell(1,lv/2); tempsym=temp;
ts=zeros(1,lv/2);
for ii=1:1:(lv/2)
  [temp{ii},tempsym{ii}]=RSC_heliCoupDcj(varargin{2*ii-1},varargin{2*ii});
  [ts(ii),~]=size(temp{ii});
end
tsm=max(ts);
for ii=1:1:(lv/2)
  temp{ii}=[temp{ii};-1000*ones(tsm-ts(ii),6)];
  tempsym{ii}=[tempsym{ii};-1000*ones(tsm-ts(ii),1)];
  uu=[uu,temp{ii}]; uu=[uu,-1000*ones(tsm,1)];
  uusym=[uusym,tempsym{ii}];
end
end

