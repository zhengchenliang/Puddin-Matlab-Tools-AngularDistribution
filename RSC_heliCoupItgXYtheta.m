function [vv,vvsym] = RSC_heliCoupItgXYtheta(XYZ,XYZsym,sxXY,sxYZ,draw,times)
%% 1-D XYZ integral #4 to #1
[mx,my]=size(XYZsym);
[~,nn]=size(XYZ);
nn=ceil(nn/7);
vv=XYZ;
vvsym=sym(zeros(mx,my));
for jj=1:1:my
  XYZpick=XYZsym(:,jj);
  for ii=1:1:mx
    if (XYZpick(ii)==-1000)
      XYZpick(ii)=0;
    end
    XYZtmp=int(XYZpick(ii),sxXY(1),-pi,pi);
    XYZtmp=int(XYZtmp,sxYZ(1),-pi,pi);
    XYZtmp=int(XYZtmp,sxYZ(2),0,pi);
    vvsym(ii,jj)=XYZtmp;
    if (draw==0)
      draw_x=0:pi/1000:pi;
      if (XYZtmp~=0)
        figure()
        plot(draw_x,RSC_F_subbs(subs(XYZtmp,sxXY(2),draw_x),times*nn),'k'); grid on
        xlabel(string(sxXY(2))); ylabel('A/I (a.u.)');
        title(['chr=',num2str(times)])
      end
    end
  end
end
end