function [vv,vvsym] = RSC_heliCoupItgYZvarphi(XYZ,XYZsym,sxXY,sxYZ,draw,times)
%% 1-D XYZ integral #4 to #1 (draw=0 to draw, times for subbs times)
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
    XYZtmp=int(XYZpick(ii),sxXY(2),0,pi);
    XYZtmp=int(XYZtmp,sxYZ(2),0,pi);
    XYZtmp=int(XYZtmp,sxXY(1),-pi,pi);
    vvsym(ii,jj)=XYZtmp;
    if (draw==0)
      draw_x=0:pi/1000:pi;
      if (XYZtmp~=0)
        figure()
        plot(draw_x,RSC_F_subbs(subs(XYZtmp,sxYZ(1),draw_x),times*nn),'k'); grid on
        xlabel(string(sxYZ(1))); ylabel('A/I (a.u.)');
        title(['chr=',num2str(times)])
      end
    end
  end
end
end