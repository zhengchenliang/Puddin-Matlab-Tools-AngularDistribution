clc; clear
%% Set Name
syms varphi_Xj varphi_Xy ...
     theta_Xj theta_Xy ...
     varphi_Jp varphi_Jm varphi_Yp varphi_Ym ...
     theta_Jp theta_Jm theta_Yp theta_Ym
%% Configure
ppXJY=[1,1,0,-1,-1,1,-1,-1,1];
ppJpm=[1,1,1,-1,-1,-1,-1,1,-1]; ppYpm=ppJpm;
sxXJ=[varphi_Xj,theta_Xj]; sxXY=[varphi_Xy,theta_Xy];
sxJp=[varphi_Jp,theta_Jp]; sxJm=[varphi_Jm,theta_Jm];
sxYp=[varphi_Yp,theta_Yp]; sxYm=[varphi_Ym,theta_Ym];
%% DF Chain
[XJp,XJpsym]=RSC_heliCoupDcjFChain(ppXJY,sxXJ,ppJpm,sxJp);
[XJm,XJmsym]=RSC_heliCoupDcjFChain(ppXJY,sxXJ,ppJpm,sxJm);
[XYp,XYpsym]=RSC_heliCoupDcjFChain(ppXJY,sxXY,ppYpm,sxYp);
[XYm,XYmsym]=RSC_heliCoupDcjFChain(ppXJY,sxXY,ppYpm,sxYm);
%% Amplitude
[XJpA,XJpAsym]=RSC_heliCoupChain2A(XJp,XJpsym);
[XJmA,XJmAsym]=RSC_heliCoupChain2A(XJm,XJmsym);
[XYpA,XYpAsym]=RSC_heliCoupChain2A(XYp,XYpsym);
[XYmA,XYmAsym]=RSC_heliCoupChain2A(XYm,XYmsym);
%% Amplitude Skim
[XJpA,XJpAsym]=RSC_heliCoupAmpSkim(XJpA,XJpAsym);
[XJmA,XJmAsym]=RSC_heliCoupAmpSkim(XJmA,XJmAsym);
[XYpA,XYpAsym]=RSC_heliCoupAmpSkim(XYpA,XYpAsym);
[XYmA,XYmAsym]=RSC_heliCoupAmpSkim(XYmA,XYmAsym);
%% Integral: Amplitude
[XJv,XJvsym]=RSC_heliCoupItgXYvarphi(XJpA,XJpAsym,sxXJ,sxJp,0,1);
[XJt,XJtsym]=RSC_heliCoupItgXYtheta(XJpA,XJpAsym,sxXJ,sxJp,0,1);
[Jpv,Jpvsym]=RSC_heliCoupItgYZvarphi(XJpA,XJpAsym,sxXJ,sxJp,0,1);
[Jpt,Jptsym]=RSC_heliCoupItgYZtheta(XJpA,XJpAsym,sxXJ,sxJp,0,1);
%% Density Matrix
if (ppXJY(3)==0)
  rho=1;
elseif (ppXJY(3)==1)
  rho=[0 0 0;
       0 1 0;
       0 0 0];
elseif (ppXJY(3)==2)
  rho=[0 0 0 0 0;
       0 0 0 0 0;
       0 0 1 0 0;
       0 0 0 0 0;
       0 0 0 0 0];
end
%% Intensity
[XJpI,XJpIsym]=RSC_heliCoupArho2I(XJpA,XJpAsym,rho);
[XJmI,XJmIsym]=RSC_heliCoupArho2I(XJmA,XJmAsym,rho);
[XYpI,XYpIsym]=RSC_heliCoupArho2I(XYpA,XYpAsym,rho);
[XYmI,XYmIsym]=RSC_heliCoupArho2I(XYmA,XYmAsym,rho);
%% Integral: Intensity
[XJvi,XJvisym]=RSC_heliCoupItgXYvarphi(XJpI,XJpIsym,sxXJ,sxJp,0,10);
[XJti,XJtisym]=RSC_heliCoupItgXYtheta(XJpI,XJpIsym,sxXJ,sxJp,0,10);
[Jpvi,Jpvisym]=RSC_heliCoupItgYZvarphi(XJpI,XJpIsym,sxXJ,sxJp,0,10);
[Jpti,Jptisym]=RSC_heliCoupItgYZtheta(XJpI,XJpIsym,sxXJ,sxJp,0,10);
