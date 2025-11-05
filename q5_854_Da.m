%Luftdruck Resonanzen ******sehr gut***************
%[y, Ts]=hol_DWD(3); %1=2000..2009; 2=2010..2019
load yDWD, Ts=3600;
L=length(y); figure (1), zeig_sp2(y,-Ts,'log');
fGW0=5.8546e-6; %HW95=fGW0*360*3600; optimiert

tD=23.93447192; fD=1/(tD*3600); 
tJ=365.25636042; fJ=1/tJ/24/3600; 
[sp,f]=zeig_sp2(y,Ts,-200000.3);
k=1;while f(k)<fGW0,k=k+1;end,k=k-15:k+15; plot(1e6*f(k),sp(k))
t=(1:Ts:21/fJ)'; t=2*pi*t(1:length(y));
if 1==100000000000000000000000, erg=zeros(65,2); 
for k3=1:65 %opt. Phasenverschiebung für AM suchen
AM=sin(t*fD+k3/10); %normieren
y1=y.*AM; [sp,~]=zeig_sp2(y1,Ts,-200000.3); %plot(f(k),sp(k))
erg(k3,1)=k3/10; erg(k3,2)=max(sp(k));
end, plot(erg(1:3,2)) %Optimum bei k3=9.7
end

y0=y;%.*sin(t*fD+0.97); %++++++++++++++++++++++++++++

fGW=fGW0-5.3e-6; %Offset
y= filtTrans_n_Hz(y0,Ts,3e4,fGW0,fGW,0.45e-6); %bis 2*200nHz wegen Mod-Index
y=decimate(decimate(y,10),10); Ts=100*Ts; zeig_sp2(y,-Ts,'log'); 
[sp,f]=zeig_sp2(y,Ts,-300000.3);k=1;while f(k)<fGW,k=k+1;end,q1=k-13:k+14;
plot(1e9*f(q1),sp(q1)), title('verschobene fGW')
t=(1:Ts:21/fJ)'; tt=1:(length(y)); fZF=1/Ts/10; %~0.3e-6;
    dftt=1e-20*t(tt); t=2*pi*t(tt); %41 Jahre FM darstellen
    tt=(tt)'*2*pi*fZF*Ts; %fixe ZF wegen variablem Oszillator
    stt=sin(tt); ctt=cos(tt);
j=0.08e-9; h1=window_sinc_filter(5e5,fZF-j,fZF+j,1/Ts,'bandpass','Blackman');
h4=window_sinc_filter(2e5,0.2e-9,1,1/Ts,'low','blackman');
h5=window_sinc_filter(1e4,10e-9,1,1/Ts,'low','blackman');
%j=5*0.08e-9;
%h6=window_sinc_filter(5e5,fZF+fJ-j,fZF+fJ+j,1/Ts,'bandpass','Blackman');
%h7=window_sinc_filter(5e5,fZF-fJ-j,fZF-fJ+j,1/Ts,'bandpass','Blackman');
ax=0.917; px=6.266; fx=40.9065e-9; D1=0.88; ww=-3.3e-04; phas=-5.92;
ax2=1.591; px2=5.08; fx2=16.777e-9; % C
ax3=0.389; px3=0.62; fx3=9.63e-10; %Frequenzzähler M
ax4=0; px4=4.8; fx4=12.2e-10; % Frequenzzähler L 
ax5=0.93; px5=2.84; fx5=4.451e-9; % D
ax6=2.316; px6=4.329; fx6=9.2185e-9; % F
ax7=0; px7=0.6712; fx7=12.3575e-10; % Frequenzzähler  Jahre
ax8=2.073; px8=0.766; fx8=57.127e-9; % 
ax9=0; px9=3; fx9=102e-9; % E
ax10=0; px10=3; fx10=6.12776e-9; % Jahre
ax11=0; px11=3; fx11=0.98e-9; % H
ax12=0; px12=3; fx12=4.8e-9; % J
ax13=0; px13=3; fx13=1.24e-9; % K
ax14=0; px14=3; fx14=1.39276e-10; % Frequenzzähler
ax15=0; px15=3; fx15=1.20464e-10; % Frequenzzähler
axJ=1.795; pxJ=0.797; fxJ=31.496e-9; %31.69 nHz

dp=0.1; dp2=0.1; dp3=1e-3; dp4=1e-3; dp5=0.1; dp6=0.1; dp7=1e-3; dp8=0.01; dpJ=0.1; dp9=0.1; dp10=0.1; dp11=0.01; dp12=0.1; dp13=0.1; dp14=1e-3; dp15=1e-3; %Schrittweite
df=3e-12; df2=3e-12; df3=1e-13; df4=1e-13; df5=3e-12; df6=3e-12; df7=1e-13; df8=3e-12; df9=3e-12; df10=3e-12; df11=3e-12; df12=3e-12; df13=3e-12; df14=1e-13; df15=1e-13; dfJ=1e-13;
da=1e-3; da2=1e-2; da3=1e-5; da4=1e-5; da5=1e-3; da6=1e-3; da7=1e-5; da8=1e-3; daJ=1e-4; da9=1e-3; da10=1e-2; da11=1e-2; da12=1e-3; da13=1e-3; da14=1e-5; da15=1e-5; dD=0.01;

R=zeros(595,52); erg=zeros(4,2); Limi=-9; k3alt=1; z0=zeros(400,3); 
for k3=1:595
z2=ax*sin(px+t*fx)+ax2*sin(px2+t*fx2)+ax3*sin(px3+t*fx3)+ax4*sin(px4+t*fx4)+ax5*sin(px5+t*fx5)+ax6*sin(px6+t*fx6)+ax7*sin(px7+t*fx7)+ax8*sin(px8+t*fx8);
z2=z2+axJ*sin(pxJ+t*fxJ)+ax9*sin(px9+t*fx9)+ax10*sin(px10+t*fx10)+ax11*sin(px11+t*fx11)+ax12*sin(px12+t*fx12)+ax13*sin(px13+t*fx13)+ax14*sin(px14+t*fx14);
z2=z2+ax15*sin(px15+t*fx15);

%PSD berechnen
%h6=window_sinc_filter(1e3,80e-9,1,1/Ts,'low','blackman');
%y1=conv(ys,h6,'same').*stt+conv(yc,h6,'same').*ctt;
%j=16*1024;[sp,f]= pwelch(cat(1,y1,zeros(3e5,1)),rectwin(j),j/8,16*j,1/Ts);plot(1e9*f,sp)
%k3=1/Ts/length(sp)/2*sum(sp(26100:26320)-0.03); k=f(26100)-f(26320); %PSD

if ax~=0, z2=z2-ax*sin(px+t*fx); %nach 3 Iterationen wieder einsetzen
for j=1:123
for k=1:3, pxx=px+dp*(k-2); %sehr tiefe Frequenz
    z=ax*sin(pxx+t*fx)+z2; 
    yc=t.*(fGW+D1*dftt)+z; ys=y.*sin(yc); yc=y.*cos(yc); 
    yc=conv(ys,h4,'same').*stt+conv(yc,h4,'same').*ctt;
    erg(k,:)=[k-2,max(yc)];
end, p=polyfit(erg(1:3,1),erg(1:3,2),2); %max suchen
if p(1)<0, k=-p(2)/2/p(1); px=px+dp*max(-19,min(19,k));
else
    if erg(1,2)<erg(3,2), px=px+4*dp;
    else, px=px-4*dp;
    end
end
if abs(k)<1
    if dp<1e-7 || max(erg(1:3,2))<1.4*min(erg(1:3,2)), break
    else, dp=dp/2;
    end
end
end

for j=1:123
for k=1:3, fxx=fx+(k-2)*df; 
    z=ax*sin(px+t*fxx)+z2; 
    yc=t.*(fGW+D1*dftt)+z; ys=y.*sin(yc); yc=y.*cos(yc); 
    yc=conv(ys,h4,'same').*stt+conv(yc,h4,'same').*ctt;
    erg(k,:)=[k-2,max(yc)];
end, p=polyfit(erg(1:3,1),erg(1:3,2),2); %max suchen
if p(1)<0, k=-p(2)/2/p(1); fx=fx+df*max(-19,min(19,k));
else
    if erg(1,2)<erg(3,2), fx=fx+4*df;
    else, fx=fx-4*df;
    end
end
if abs(k)<1
    if df<2e-18 || max(erg(1:3,2))<1.4*min(erg(1:3,2)), break
    else, df=df/2;
    end
end
end

for j=1:123
for k=1:3, axx=ax+da*(k-2);
    z=axx*sin(px+t*fx)+z2; 
    yc=t.*(fGW+D1*dftt)+z; ys=y.*sin(yc); yc=y.*cos(yc); 
    yc=conv(ys,h4,'same').*stt+conv(yc,h4,'same').*ctt;
    erg(k,:)=[k-2,max(yc)];
end, p=polyfit(erg(1:3,1),erg(1:3,2),2); %max suchen
if p(1)<0, k=-p(2)/2/p(1); ax=ax+da*max(-19,min(19,k));
else
    if erg(1,2)<erg(3,2), ax=ax+4*da;
    else, ax=ax-4*da;
    end
end
if abs(k)<1
    if da<1e-6 || max(erg(1:3,2))<1.4*min(erg(1:3,2)), break
    else, da=da/2;
    end
end
end%,jjj=log10(erg(2,2)); %fx.. ist abgeschlossen
z2=z2+ax*sin(px+t*fx); %nach 3 Iterationen wieder einsetzen
end

if axJ~=0, z2=z2-axJ*sin(pxJ+t*fxJ); %nach 3 Iterationen wieder einsetzen
for j=1:123
for k=1:3, pxx=pxJ+dpJ*(k-2);
    z=axJ*sin(pxx+t*fxJ)+z2; 
    yc=t.*(fGW+D1*dftt)+z; ys=y.*sin(yc); yc=y.*cos(yc); 
    yc=conv(ys,h4,'same').*stt+conv(yc,h4,'same').*ctt;
    erg(k,:)=[k-2,max(yc)];
end, p=polyfit(erg(1:3,1),erg(1:3,2),2); %max suchen
if p(1)<0, k=-p(2)/2/p(1); pxJ=pxJ+dpJ*max(-19,min(19,k));
else
    if erg(1,2)<erg(3,2), pxJ=pxJ+4*dpJ;
    else,  pxJ=pxJ-4*dpJ;
    end
end
if abs(k)<1
    if dpJ<1e-5 || max(erg(1:3,2))<1.4*min(erg(1:3,2)), break
    else, dpJ=dpJ/2;
    end
end
end

for j=1:123
for k=1:3, fxx=fxJ+(k-2)*dfJ; 
    z=axJ*sin(pxJ+t*fxx)+z2; 
    yc=t.*(fGW+D1*dftt)+z; ys=y.*sin(yc); yc=y.*cos(yc); 
    yc=conv(ys,h4,'same').*stt+conv(yc,h4,'same').*ctt;
    erg(k,:)=[k-2,max(yc)];
end, p=polyfit(erg(1:3,1),erg(1:3,2),2); %max suchen
if p(1)<0, k=-p(2)/2/p(1); fxJ=fxJ+dfJ*max(-19,min(19,k));
else
    if erg(1,2)<erg(3,2), fxJ=fxJ+4*dfJ;
    else,  fxJ=fxJ-4*dfJ;
    end
end
if abs(k)<1
    if dfJ<2e-15 || max(erg(1:3,2))<1.4*min(erg(1:3,2)), break
    else, dfJ=dfJ/2;
    end
end
end

for j=1:123
for k=1:3, axx=axJ+daJ*(k-2);
    z=axx*sin(pxJ+t*fxJ)+z2; 
    yc=t.*(fGW+D1*dftt)+z; ys=y.*sin(yc); yc=y.*cos(yc); 
    yc=conv(ys,h4,'same').*stt+conv(yc,h4,'same').*ctt;
    erg(k,:)=[k-2,max(yc)];
end, p=polyfit(erg(1:3,1),erg(1:3,2),2); %max suchen
if p(1)<0, k=-p(2)/2/p(1); axJ=axJ+daJ*max(-19,min(19,k));
else
    if erg(1,2)<erg(3,2), axJ=axJ+4*daJ;
    else, axJ=axJ-4*daJ;
    end
end
if abs(k)<1
    if daJ<1e-6 || max(erg(1:3,2))<1.4*min(erg(1:3,2)), break
    else, daJ=daJ/2;
    end
end
end %fxJ.. ist abgeschlossen
z2=z2+axJ*sin(pxJ+t*fxJ); %nach 3 Iterationen wieder einsetzen
end


if ax2~=0, z2=z2-ax2*sin(px2+t*fx2); %nach 3 Iterationen wieder einsetzen
for j=1:123
for k=1:3, pxx=px2+dp2*(k-2);
    z=ax2*sin(pxx+t*fx2)+z2; 
    yc=t.*(fGW+D1*dftt)+z; ys=y.*sin(yc); yc=y.*cos(yc); 
    yc=conv(ys,h4,'same').*stt+conv(yc,h4,'same').*ctt;
    erg(k,:)=[k-2,max(yc)];
end, p=polyfit(erg(1:3,1),erg(1:3,2),2); %max suchen
if p(1)<0, k=-p(2)/2/p(1); px2=px2+dp2*max(-19,min(19,k));
else
    if erg(1,2)<erg(3,2), px2=px2+4*dp2;
    else,  px2=px2-4*dp2;
    end
end
if abs(k)<1
    if dp2<1e-5 || max(erg(1:3,2))<1.4*min(erg(1:3,2)), break
    else, dp2=dp2/2;
    end
end
end

for j=1:123 %
for k=1:3, fxx=fx2+(k-2)*df2; 
    z=ax2*sin(px2+t*fxx)+z2; 
    yc=t.*(fGW+D1*dftt)+z; ys=y.*sin(yc); yc=y.*cos(yc); 
    yc=conv(ys,h4,'same').*stt+conv(yc,h4,'same').*ctt;
    erg(k,:)=[k-2,max(yc)];
end, p=polyfit(erg(1:3,1),erg(1:3,2),2); %max suchen
if p(1)<0, k=-p(2)/2/p(1); fx2=fx2+df2*max(-19,min(19,k));
else
    if erg(1,2)<erg(3,2), fx2=fx2+4*df2;
    else,  fx2=fx2-4*df2;
    end
end
if abs(k)<1
    if df2<2e-15 || max(erg(1:3,2))<1.4*min(erg(1:3,2)), break
    else, df2=df2/2;
    end
end
end

for j=1:123
for k=1:3, axx=ax2+da2*(k-2);
    z=axx*sin(px2+t*fx2)+z2; 
    yc=t.*(fGW+D1*dftt)+z; ys=y.*sin(yc); yc=y.*cos(yc); 
    yc=conv(ys,h4,'same').*stt+conv(yc,h4,'same').*ctt;
    erg(k,:)=[k-2,max(yc)];
end, p=polyfit(erg(1:3,1),erg(1:3,2),2); %max suchen
if p(1)<0, k=-p(2)/2/p(1); ax2=ax2+da2*max(-19,min(19,k));
else
    if erg(1,2)<erg(3,2), ax2=ax2+4*da2;
    else,  ax2=ax2-4*da2;
    end
end
if abs(k)<1
    if da2<1e-5 || max(erg(1:3,2))<1.4*min(erg(1:3,2)), break
    else, da2=da2/2;
    end
end
end%,jjj=log10(erg(2,2)); %fx2.. ist abgeschlossen
z2=z2+ax2*sin(px2+t*fx2); %nach 3 Iterationen wieder einsetzen
end

if ax3~=0, z2=z2-ax3*sin(px3+t*fx3); %nach 3 Iterationen wieder einsetzen
for j=1:123
for k=1:3, pxx=px3+dp3*(k-2);
    z=ax3*sin(pxx+t*fx3)+z2; 
    yc=t.*(fGW+D1*dftt)+z; ys=y.*sin(yc); yc=y.*cos(yc); 
    yc=conv(ys,h5,'same').*stt+conv(yc,h5,'same').*ctt;
    z0=T_zaehl_T(conv(yc,h1,'same'),Ts); %plot(z0(:,2))
    z0(:,3)=z0(:,2)-mean(z0(:,2));
    [p,~,mu]=polyfit(z0(:,1),z0(:,3),2);
    z0(:,3)=polyval(p,z0(:,1),[],mu); %plot(z0(:,1),z0(:,3),z0(:,1),y1)
    erg(k,:)=[k-2,sum(z0(:,3).^2)];
end, p=polyfit(erg(1:3,1),erg(1:3,2),2); %min suchen
if p(1)>0, k=-p(2)/2/p(1); px3=px3+dp3*max(-19,min(19,k));
else
    if erg(1,2)>erg(3,2), px3=px3+4*dp3;
    else,  px3=px3-4*dp3;
    end
end
if abs(k)<1
    if dp3<1e-7 || max(erg(1:3,2))<1.4*min(erg(1:3,2)), break
    else, dp3=dp3/2;
    end
end
end

for j=1:123
for k=1:3, fxx=fx3+(k-2)*df3; 
    z=ax3*sin(px3+t*fxx)+z2; 
    yc=t.*(fGW+D1*dftt)+z; ys=y.*sin(yc); yc=y.*cos(yc); 
    yc=conv(ys,h5,'same').*stt+conv(yc,h5,'same').*ctt;
    z0=T_zaehl_T(conv(yc,h1,'same'),Ts); %plot(z0(:,2))
    z0(:,3)=z0(:,2)-mean(z0(:,2));
    [p,~,mu]=polyfit(z0(:,1),z0(:,3),2);
    z0(:,3)=polyval(p,z0(:,1),[],mu); %plot(z0(:,1),z0(:,3),z0(:,1),y1)
    erg(k,:)=[k-2,sum(z0(:,3).^2)];
end, p=polyfit(erg(1:3,1),erg(1:3,2),2); %min suchen
if p(1)>0, k=-p(2)/2/p(1); fx3=fx3+df3*max(-19,min(19,k));
else
    if erg(1,2)>erg(3,2), fx3=fx3+4*df3;
    else,  fx3=fx3-4*df3;
    end
end
if abs(k)<1
    if df3<1e-17 || max(erg(1:3,2))<1.4*min(erg(1:3,2)), break
    else, df3=df3/2;
    end
end
end

for j=1:123
for k=1:3, axx=ax3+da3*(k-2);
    z=axx*sin(px3+t*fx3)+z2; 
    yc=t.*(fGW+D1*dftt)+z; ys=y.*sin(yc); yc=y.*cos(yc); 
yc=conv(ys,h5,'same').*stt+conv(yc,h5,'same').*ctt;
    z0=T_zaehl_T(conv(yc,h1,'same'),Ts); %plot(z0(:,2))
    z0(:,3)=z0(:,2)-mean(z0(:,2));
    [p,~,mu]=polyfit(z0(:,1),z0(:,3),2);
    z0(:,3)=polyval(p,z0(:,1),[],mu); %plot(z0(:,1),z0(:,3),z0(:,1),z0(:,3))
    erg(k,:)=[k-2,sum(z0(:,3).^2)];
end, p=polyfit(erg(1:3,1),erg(1:3,2),2); %min suchen
if p(1)>0, k=-p(2)/2/p(1); ax3=ax3+da3*max(-19,min(19,k));
else
    if erg(1,2)>erg(3,2), ax3=ax3+4*da3;
    else,  ax3=ax3-4*da3;
    end
end
if abs(k)<1
    if da3<1e-7 || max(erg(1:3,2))<1.4*min(erg(1:3,2)), break
    else, da3=da3/2;
    end
end
end, jjj=log10(erg(2,2)); %fx3.. ist abgeschlossen
z2=z2+ax3*sin(px3+t*fx3); %nach 3 Iterationen wieder einsetzen
end

if ax4~=0, z2=z2-ax4*sin(px4+t*fx4); %nach 3 Iterationen wieder einsetzen
for j=1:123
for k=1:3, pxx=px4+dp4*(k-2);
    z=ax4*sin(pxx+t*fx4)+z2; 
    yc=t.*(fGW+D1*dftt)+z; ys=y.*sin(yc); yc=y.*cos(yc); 
yc=conv(ys,h5,'same').*stt+conv(yc,h5,'same').*ctt;
    z0=T_zaehl_T(conv(yc,h1,'same'),Ts); %plot(z0(:,2))
    z0(:,3)=z0(:,2)-mean(z0(:,2));
    [p,~,mu]=polyfit(z0(:,1),z0(:,3),2);
    z0(:,3)=polyval(p,z0(:,1),[],mu); %plot(z0(:,1),z0(:,3),z0(:,1),z0(:,3))
    erg(k,:)=[k-2,sum(z0(:,3).^2)];
end, p=polyfit(erg(1:3,1),erg(1:3,2),2); %min suchen
if p(1)>0, k=-p(2)/2/p(1); px4=px4+dp4*max(-19,min(19,k));
else
    if erg(1,2)>erg(3,2), px4=px4+4*dp4;
    else,  px4=px4-4*dp4;
    end
end
if abs(k)<1
    if dp4<1e-7 || max(erg(1:3,2))<1.4*min(erg(1:3,2)), break
    else, dp4=dp4/2;
    end
end
end

for j=1:123
for k=1:3, fxx=fx4+(k-2)*df4; 
    z=ax4*sin(px4+t*fxx)+z2; 
    yc=t.*(fGW+D1*dftt)+z; ys=y.*sin(yc); yc=y.*cos(yc); 
yc=conv(ys,h5,'same').*stt+conv(yc,h5,'same').*ctt;
    z0=T_zaehl_T(conv(yc,h1,'same'),Ts); %plot(z0(:,2))
    z0(:,3)=z0(:,2)-mean(z0(:,2));
    [p,~,mu]=polyfit(z0(:,1),z0(:,3),2);
    z0(:,3)=polyval(p,z0(:,1),[],mu); %plot(z0(:,1),z0(:,3),z0(:,1),z0(:,3))
    erg(k,:)=[k-2,sum(z0(:,3).^2)];
end, p=polyfit(erg(1:3,1),erg(1:3,2),2); %min suchen
if p(1)>0, k=-p(2)/2/p(1); fx4=fx4+df4*max(-19,min(19,k));
else
    if erg(1,2)>erg(3,2), fx4=fx4+4*df4;
    else,  fx4=fx4-4*df4;
    end
end
if abs(k)<1
    if df4<1e-18 || max(erg(1:3,2))<1.4*min(erg(1:3,2)), break
    else, df4=df4/2;
    end
end
end

for j=1:123
for k=1:3, axx=ax4+da4*(k-2);
    z=axx*sin(px4+t*fx4)+z2; 
    yc=t.*(fGW+D1*dftt)+z; ys=y.*sin(yc); yc=y.*cos(yc); 
yc=conv(ys,h5,'same').*stt+conv(yc,h5,'same').*ctt;
    z0=T_zaehl_T(conv(yc,h1,'same'),Ts); %plot(z0(:,2))
    z0(:,3)=z0(:,2)-mean(z0(:,2));
    [p,~,mu]=polyfit(z0(:,1),z0(:,3),2);
    z0(:,3)=polyval(p,z0(:,1),[],mu); %plot(z0(:,1),z0(:,3),z0(:,1),z0(:,3))
    erg(k,:)=[k-2,sum(z0(:,3).^2)];
end, p=polyfit(erg(1:3,1),erg(1:3,2),2); %min suchen
if p(1)>0, k=-p(2)/2/p(1); ax4=ax4+da4*max(-19,min(19,k));
else
    if erg(1,2)>erg(3,2), ax4=ax4+4*da4;
    else,  ax4=ax4-4*da4;
    end
end
if abs(k)<1
    if da4<1e-7 || max(erg(1:3,2))<1.4*min(erg(1:3,2)), break
    else, da4=da4/2;
    end
end
end,jjj=log10(erg(2,2)); %fx4.. ist abgeschlossen
z2=z2+ax4*sin(px4+t*fx4); %nach 3 Iterationen wieder einsetzen
end


if ax5~=0, z2=z2-ax5*sin(px5+t*fx5); %nach 3 Iterationen wieder einsetzen
for j=1:123
for k=1:3, pxx=px5+dp5*(k-2);
    z=ax5*sin(pxx+t*fx5)+z2; 
    yc=t.*(fGW+D1*dftt)+z; ys=y.*sin(yc); yc=y.*cos(yc); 
    yc=conv(ys,h4,'same').*stt+conv(yc,h4,'same').*ctt;
    erg(k,:)=[k-2,max(yc)];
end, p=polyfit(erg(1:3,1),erg(1:3,2),2); %max suchen
if p(1)<0, k=-p(2)/2/p(1); px5=px5+dp5*max(-19,min(19,k));
else
    if erg(1,2)<erg(3,2), px5=px5+4*dp5;
    else,  px5=px5-4*dp5;
    end
end
if abs(k)<1
    if dp5<1e-5 || max(erg(1:3,2))<1.4*min(erg(1:3,2)), break
    else, dp5=dp5/2;
    end
end
end

for j=1:123
for k=1:3, fxx=fx5+(k-2)*df5; 
    z=ax5*sin(px5+t*fxx)+z2; 
    yc=t.*(fGW+D1*dftt)+z; ys=y.*sin(yc); yc=y.*cos(yc); 
    yc=conv(ys,h4,'same').*stt+conv(yc,h4,'same').*ctt;
    erg(k,:)=[k-2,max(yc)];
end, p=polyfit(erg(1:3,1),erg(1:3,2),2); %max suchen
if p(1)<0, k=-p(2)/2/p(1); fx5=fx5+df5*max(-19,min(19,k));
else
    if erg(1,2)<erg(3,2), fx5=fx5+4*df5;
    else,  fx5=fx5-4*df5;
    end
end
if abs(k)<1
    if df5<2e-15 || max(erg(1:3,2))<1.4*min(erg(1:3,2)), break
    else, df5=df5/2;
    end
end
end

for j=1:123
for k=1:3, axx=ax5+da5*(k-2);
    z=axx*sin(px5+t*fx5)+z2; 
    yc=t.*(fGW+D1*dftt)+z; ys=y.*sin(yc); yc=y.*cos(yc); 
    yc=conv(ys,h4,'same').*stt+conv(yc,h4,'same').*ctt;
    erg(k,:)=[k-2,max(yc)];
end, p=polyfit(erg(1:3,1),erg(1:3,2),2); %max suchen
if p(1)<0, k=-p(2)/2/p(1); ax5=ax5+da5*max(-19,min(19,k));
else
    if erg(1,2)<erg(3,2), ax5=ax5+4*da5;
    else,  ax5=ax5-4*da5;
    end
end
if abs(k)<1
    if da5<1e-5 || max(erg(1:3,2))<1.4*min(erg(1:3,2)), break
    else, da5=da5/2;
    end
end
end %fx5.. ist abgeschlossen
z2=z2+ax5*sin(px5+t*fx5); %nach 3 Iterationen wieder einsetzen
end

if ax6~=0, z2=z2-ax6*sin(px6+t*fx6); %nach 3 Iterationen wieder einsetzen
for j=1:123
for k=1:3, pxx=px6+dp6*(k-2);
    z=ax6*sin(pxx+t*fx6)+z2; 
    yc=t.*(fGW+D1*dftt)+z; ys=y.*sin(yc); yc=y.*cos(yc); 
    yc=conv(ys,h4,'same').*stt+conv(yc,h4,'same').*ctt;
    erg(k,:)=[k-2,max(yc)];
end, p=polyfit(erg(1:3,1),erg(1:3,2),2); %max suchen
if p(1)<0, k=-p(2)/2/p(1); px6=px6+dp6*max(-19,min(19,k));
else
    if erg(1,2)<erg(3,2), px6=px6+4*dp6;
    else,  px6=px6-4*dp6;
    end
end
if abs(k)<1
    if dp6<1e-5 || max(erg(1:3,2))<1.4*min(erg(1:3,2)), break
    else, dp6=dp6/2;
    end
end
end

for j=1:123
for k=1:3, fxx=fx6+(k-2)*df6; 
    z=ax6*sin(px6+t*fxx)+z2; 
    yc=t.*(fGW+D1*dftt)+z; ys=y.*sin(yc); yc=y.*cos(yc); 
    yc=conv(ys,h4,'same').*stt+conv(yc,h4,'same').*ctt;
    erg(k,:)=[k-2,max(yc)];
end, p=polyfit(erg(1:3,1),erg(1:3,2),2); %max suchen
if p(1)<0, k=-p(2)/2/p(1); fx6=fx6+df6*max(-19,min(19,k));
else
    if erg(1,2)<erg(3,2), fx6=fx6+4*df6;
    else,  fx6=fx6-4*df6;
    end
end
if abs(k)<1
    if df6<2e-15 || max(erg(1:3,2))<1.4*min(erg(1:3,2)), break
    else, df6=df6/2;
    end
end
end

for j=1:123
for k=1:3, axx=ax6+da6*(k-2);
    z=axx*sin(px6+t*fx6)+z2; 
    yc=t.*(fGW+D1*dftt)+z; ys=y.*sin(yc); yc=y.*cos(yc); 
    yc=conv(ys,h4,'same').*stt+conv(yc,h4,'same').*ctt;
    erg(k,:)=[k-2,max(yc)];
end, p=polyfit(erg(1:3,1),erg(1:3,2),2); %max suchen
if p(1)<0, k=-p(2)/2/p(1); ax6=ax6+da6*max(-19,min(19,k));
else
    if erg(1,2)<erg(3,2), ax6=ax6+4*da6;
    else,  ax6=ax6-4*da6;
    end
end
if abs(k)<1
    if da6<1e-5 || max(erg(1:3,2))<1.4*min(erg(1:3,2)), break
    else, da6=da6/2;
    end
end
end %fx6.. ist abgeschlossen
z2=z2+ax6*sin(px6+t*fx6); %nach 3 Iterationen wieder einsetzen
end


if ax7~=0, z2=z2-ax7*sin(px7+t*fx7); %nach 3 Iterationen wieder einsetzen
for j=1:123
for k=1:3, pxx=px7+dp7*(k-2);
    z=ax7*sin(pxx+t*fx7)+z2; 
    yc=t.*(fGW+D1*dftt)+z; ys=y.*sin(yc); yc=y.*cos(yc); 
yc=conv(ys,h5,'same').*stt+conv(yc,h5,'same').*ctt;
    z0=T_zaehl_T(conv(yc,h1,'same'),Ts); %plot(z0(:,2))
    z0(:,3)=z0(:,2)-mean(z0(:,2));
    [p,~,mu]=polyfit(z0(:,1),z0(:,3),2);
    z0(:,3)=polyval(p,z0(:,1),[],mu); %plot(z0(:,1),z0(:,3),z0(:,1),z0(:,3))
    erg(k,:)=[k-2,sum(z0(:,3).^2)];
end, p=polyfit(erg(1:3,1),erg(1:3,2),2); %min suchen
if p(1)>0, k=-p(2)/2/p(1); px7=px7+dp7*max(-19,min(19,k));
else
    if erg(1,2)>erg(3,2), px7=px7+4*dp7;
    else,  px7=px7-4*dp7;
    end
end
if abs(k)<1
    if dp7<1e-7 || max(erg(1:3,2))<1.4*min(erg(1:3,2)), break
    else, dp7=dp7/2;
    end
end
end

for j=1:123
for k=1:3, fxx=fx7+(k-2)*df7; 
    z=ax7*sin(px7+t*fxx)+z2; 
    yc=t.*(fGW+D1*dftt)+z; ys=y.*sin(yc); yc=y.*cos(yc); 
yc=conv(ys,h5,'same').*stt+conv(yc,h5,'same').*ctt;
    z0=T_zaehl_T(conv(yc,h1,'same'),Ts); %plot(z0(:,2))
    z0(:,3)=z0(:,2)-mean(z0(:,2));
    [p,~,mu]=polyfit(z0(:,1),z0(:,3),2);
    z0(:,3)=polyval(p,z0(:,1),[],mu); %plot(z0(:,1),z0(:,3),z0(:,1),z0(:,3))
    erg(k,:)=[k-2,sum(z0(:,3).^2)];
end, p=polyfit(erg(1:3,1),erg(1:3,2),2); %min suchen
if p(1)>0, k=-p(2)/2/p(1); fx7=fx7+df7*max(-19,min(19,k));
else
    if erg(1,2)>erg(3,2), fx7=fx7+4*df7;
    else,  fx7=fx7-4*df7;
    end
end
if abs(k)<1
    if df7<2e-17 || max(erg(1:3,2))<1.4*min(erg(1:3,2)), break
    else, df7=df7/2;
    end
end
end

for j=1:123
for k=1:3, axx=ax7+da7*(k-2);
    z=axx*sin(px7+t*fx7)+z2; 
    yc=t.*(fGW+D1*dftt)+z; ys=y.*sin(yc); yc=y.*cos(yc); 
yc=conv(ys,h5,'same').*stt+conv(yc,h5,'same').*ctt;
    z0=T_zaehl_T(conv(yc,h1,'same'),Ts); %plot(z0(:,2))
    z0(:,3)=z0(:,2)-mean(z0(:,2));
    [p,~,mu]=polyfit(z0(:,1),z0(:,3),2);
    z0(:,3)=polyval(p,z0(:,1),[],mu); %plot(z0(:,1),z0(:,3),z0(:,1),z0(:,3))
    erg(k,:)=[k-2,sum(z0(:,3).^2)];
end, p=polyfit(erg(1:3,1),erg(1:3,2),2); %min suchen
if p(1)>0, k=-p(2)/2/p(1); ax7=ax7+da7*max(-19,min(19,k));
else
    if erg(1,2)>erg(3,2), ax7=ax7+4*da7;
    else,  ax7=ax7-4*da7;
    end
end
if abs(k)<1
    if da7<1e-7 || max(erg(1:3,2))<1.4*min(erg(1:3,2)), break
    else, da7=da7/2;
    end
end
end, jjj=log10(erg(2,2)); %fx7.. ist abgeschlossen
z2=z2+ax7*sin(px7+t*fx7); %nach 3 Iterationen wieder einsetzen
end


if ax8~=0, z2=z2-ax8*sin(px8+t*fx8); %nach 3 Iterationen wieder einsetzen
for j=1:123
for k=1:3, pxx=px8+dp8*(k-2);
    z=ax8*sin(pxx+t*fx8)+z2; 
    yc=t.*(fGW+D1*dftt)+z; ys=y.*sin(yc); yc=y.*cos(yc); 
    yc=conv(ys,h4,'same').*stt+conv(yc,h4,'same').*ctt;
    erg(k,:)=[k-2,max(yc)];
end, p=polyfit(erg(1:3,1),erg(1:3,2),2); %max suchen
if p(1)<0, k=-p(2)/2/p(1); px8=px8+dp8*max(-19,min(19,k));
else
    if erg(1,2)<erg(3,2), px8=px8+4*dp8;
    else,  px8=px8-4*dp8;
    end
end
if abs(k)<1
    if dp8<1e-5 || max(erg(1:3,2))<1.4*min(erg(1:3,2)), break
    else, dp8=dp8/2;
    end
end
end

for j=1:123
for k=1:3, fxx=fx8+(k-2)*df8; 
    z=ax8*sin(px8+t*fxx)+z2; 
    yc=t.*(fGW+D1*dftt)+z; ys=y.*sin(yc); yc=y.*cos(yc); 
    yc=conv(ys,h4,'same').*stt+conv(yc,h4,'same').*ctt;
    erg(k,:)=[k-2,max(yc)];
end, p=polyfit(erg(1:3,1),erg(1:3,2),2); %max suchen
if p(1)<0, k=-p(2)/2/p(1); fx8=fx8+df8*max(-19,min(19,k));
else
    if erg(1,2)<erg(3,2), fx8=fx8+4*df8;
    else,  fx8=fx8-4*df8;
    end
end
if abs(k)<1
    if df8<2e-15 || max(erg(1:3,2))<1.4*min(erg(1:3,2)), break
    else, df8=df8/2;
    end
end
end

for j=1:123
for k=1:3, axx=ax8+da8*(k-2);
    z=axx*sin(px8+t*fx8)+z2; 
    yc=t.*(fGW+D1*dftt)+z; ys=y.*sin(yc); yc=y.*cos(yc); 
    yc=conv(ys,h4,'same').*stt+conv(yc,h4,'same').*ctt;
    erg(k,:)=[k-2,max(yc)];
end, p=polyfit(erg(1:3,1),erg(1:3,2),2); %max suchen
if p(1)<0, k=-p(2)/2/p(1); ax8=ax8+da8*max(-19,min(19,k));
else
    if erg(1,2)<erg(3,2), ax8=ax8+4*da8;
    else,  ax8=ax8-4*da8;
    end
end
if abs(k)<1
    if da8<1e-5 || max(erg(1:3,2))<1.4*min(erg(1:3,2)), break
    else, da8=da8/2;
    end
end
end%,jjj=log10(erg(2,2)); %fx8.. ist abgeschlossen
z2=z2+ax8*sin(px8+t*fx8); %nach 3 Iterationen wieder einsetzen
end


if ax9~=0, z2=z2-ax9*sin(px9+t*fx9); %nach 3 Iterationen wieder einsetzen
for j=1:123
for k=1:3, pxx=px9+dp9*(k-2);
    z=ax9*sin(pxx+t*fx9)+z2; 
    yc=t.*(fGW+D1*dftt)+z; ys=y.*sin(yc); yc=y.*cos(yc); 
    yc=conv(ys,h4,'same').*stt+conv(yc,h4,'same').*ctt;
    erg(k,:)=[k-2,max(yc)];
end, p=polyfit(erg(1:3,1),erg(1:3,2),2); %max suchen
if p(1)<0, k=-p(2)/2/p(1); px9=px9+dp9*max(-19,min(19,k));
else
    if erg(1,2)<erg(3,2), px9=px9+4*dp9;
    else,  px9=px9-4*dp9;
    end
end
if abs(k)<1
    if dp9<1e-5 || max(erg(1:3,2))<1.4*min(erg(1:3,2)), break
    else, dp9=dp9/2;
    end
end
end

for j=1:123
for k=1:3, fxx=fx9+(k-2)*df9; 
    z=ax9*sin(px9+t*fxx)+z2; 
    yc=t.*(fGW+D1*dftt)+z; ys=y.*sin(yc); yc=y.*cos(yc); 
    yc=conv(ys,h4,'same').*stt+conv(yc,h4,'same').*ctt;
    erg(k,:)=[k-2,max(yc)];
end, p=polyfit(erg(1:3,1),erg(1:3,2),2); %max suchen
if p(1)<0, k=-p(2)/2/p(1); fx9=fx9+df9*max(-19,min(19,k));
else
    if erg(1,2)<erg(3,2), fx9=fx9+4*df9;
    else,  fx9=fx9-4*df9;
    end
end
if abs(k)<1
    if df9<2e-15 || max(erg(1:3,2))<1.4*min(erg(1:3,2)), break
    else, df9=df9/2;
    end
end
end

for j=1:123
for k=1:3, axx=ax9+da9*(k-2);
    z=axx*sin(px9+t*fx9)+z2; 
    yc=t.*(fGW+D1*dftt)+z; ys=y.*sin(yc); yc=y.*cos(yc); 
    yc=conv(ys,h4,'same').*stt+conv(yc,h4,'same').*ctt;
    erg(k,:)=[k-2,max(yc)];
end, p=polyfit(erg(1:3,1),erg(1:3,2),2); %max suchen
if p(1)<0, k=-p(2)/2/p(1); ax9=ax9+da9*max(-19,min(19,k));
else
    if erg(1,2)<erg(3,2), ax9=ax9+4*da9;
    else,  ax9=ax9-4*da9;
    end
end
if abs(k)<1
    if da9<1e-5 || max(erg(1:3,2))<1.4*min(erg(1:3,2)), break
    else, da9=da9/2;
    end
end
end%,jjj=log10(erg(2,2)); %fx9.. ist abgeschlossen
z2=z2+ax9*sin(px9+t*fx9); %nach 3 Iterationen wieder einsetzen
end


if ax10~=0, z2=z2-ax10*sin(px10+t*fx10); %nach 3 Iterationen wieder einsetzen
for j=1:123
for k=1:3, pxx=px10+dp10*(k-2);
    z=ax10*sin(pxx+t*fx10)+z2; 
    yc=t.*(fGW+D1*dftt)+z; ys=y.*sin(yc); yc=y.*cos(yc); 
    yc=conv(ys,h4,'same').*stt+conv(yc,h4,'same').*ctt;
    erg(k,:)=[k-2,max(yc)];
end, p=polyfit(erg(1:3,1),erg(1:3,2),2); %max suchen
if p(1)<0, k=-p(2)/2/p(1); px10=px10+dp10*max(-19,min(19,k));
else
    if erg(1,2)<erg(3,2), px10=px10+4*dp10;
    else,  px10=px10-4*dp10;
    end
end
if abs(k)<1
    if dp10<1e-5 || max(erg(1:3,2))<1.4*min(erg(1:3,2)), break
    else, dp10=dp10/2;
    end
end
end

for j=1:123
for k=1:3, fxx=fx10+(k-2)*df10; 
    z=ax10*sin(px10+t*fxx)+z2; 
    yc=t.*(fGW+D1*dftt)+z; ys=y.*sin(yc); yc=y.*cos(yc); 
    yc=conv(ys,h4,'same').*stt+conv(yc,h4,'same').*ctt;
    erg(k,:)=[k-2,max(yc)];
end, p=polyfit(erg(1:3,1),erg(1:3,2),2); %max suchen
if p(1)<0, k=-p(2)/2/p(1); fx10=fx10+df10*max(-19,min(19,k));
else
    if erg(1,2)<erg(3,2), fx10=fx10+4*df10;
    else,  fx10=fx10-4*df10;
    end
end
if abs(k)<1
    if df10<2e-15 || max(erg(1:3,2))<1.4*min(erg(1:3,2)), break
    else, df10=df10/2;
    end
end
end

for j=1:123
for k=1:3, axx=ax10+da10*(k-2);
    z=axx*sin(px10+t*fx10)+z2; 
    yc=t.*(fGW+D1*dftt)+z; ys=y.*sin(yc); yc=y.*cos(yc); 
    yc=conv(ys,h4,'same').*stt+conv(yc,h4,'same').*ctt;
    erg(k,:)=[k-2,max(yc)];
end, p=polyfit(erg(1:3,1),erg(1:3,2),2); %max suchen
if p(1)<0, k=-p(2)/2/p(1); ax10=ax10+da10*max(-19,min(19,k));
else
    if erg(1,2)<erg(3,2), ax10=ax10+4*da10;
    else,  ax10=ax10-4*da10;
    end
end
if abs(k)<1
    if da10<1e-5 || max(erg(1:3,2))<1.4*min(erg(1:3,2)), break
    else, da10=da10/2;
    end
end
end%,jjj=log10(erg(2,2)); %fx10.. ist abgeschlossen
z2=z2+ax10*sin(px10+t*fx10); %nach 3 Iterationen wieder einsetzen
end


if ax11~=0, z2=z2-ax11*sin(px11+t*fx11); %nach 3 Iterationen wieder einsetzen
for j=1:123
for k=1:3, pxx=px11+dp11*(k-2);
    z=ax11*sin(pxx+t*fx11)+z2; 
    yc=t.*(fGW+D1*dftt)+z; ys=y.*sin(yc); yc=y.*cos(yc); 
    yc=conv(ys,h4,'same').*stt+conv(yc,h4,'same').*ctt;
    erg(k,:)=[k-2,max(yc)];
end, p=polyfit(erg(1:3,1),erg(1:3,2),2); %max suchen
if p(1)<0, k=-p(2)/2/p(1); px11=px11+dp11*max(-19,min(19,k));
else
    if erg(1,2)<erg(3,2), px11=px11+4*dp11;
    else,  px11=px11-4*dp11;
    end
end
if abs(k)<1
    if dp11<1e-5 || max(erg(1:3,2))<1.4*min(erg(1:3,2)), break
    else, dp11=dp11/2;
    end
end
end

for j=1:123
for k=1:3, fxx=fx11+(k-2)*df11; 
    z=ax11*sin(px11+t*fxx)+z2; 
    yc=t.*(fGW+D1*dftt)+z; ys=y.*sin(yc); yc=y.*cos(yc); 
    yc=conv(ys,h4,'same').*stt+conv(yc,h4,'same').*ctt;
    erg(k,:)=[k-2,max(yc)];
end, p=polyfit(erg(1:3,1),erg(1:3,2),2); %max suchen
if p(1)<0, k=-p(2)/2/p(1); fx11=fx11+df11*max(-19,min(19,k));
else
    if erg(1,2)<erg(3,2), fx11=fx11+4*df11;
    else,  fx11=fx11-4*df11;
    end
end
if abs(k)<1
    if df11<2e-15 || max(erg(1:3,2))<1.4*min(erg(1:3,2)), break
    else, df11=df11/2;
    end
end
end

for j=1:123
for k=1:3, axx=ax11+da11*(k-2);
    z=axx*sin(px11+t*fx11)+z2; 
    yc=t.*(fGW+D1*dftt)+z; ys=y.*sin(yc); yc=y.*cos(yc); 
    yc=conv(ys,h4,'same').*stt+conv(yc,h4,'same').*ctt;
    erg(k,:)=[k-2,max(yc)];
end, p=polyfit(erg(1:3,1),erg(1:3,2),2); %max suchen
if p(1)<0, k=-p(2)/2/p(1); ax11=ax11+da11*max(-19,min(19,k));
else
    if erg(1,2)<erg(3,2), ax11=ax11+4*da11;
    else,  ax11=ax11-4*da11;
    end
end
if abs(k)<1
    if da11<1e-5 || max(erg(1:3,2))<1.4*min(erg(1:3,2)), break
    else, da11=da11/2;
    end
end
end%,jjj=log10(erg(2,2)); %fx11.. ist abgeschlossen
z2=z2+ax11*sin(px11+t*fx11); %nach 3 Iterationen wieder einsetzen
end


if ax12~=0, z2=z2-ax12*sin(px12+t*fx12); %nach 3 Iterationen wieder einsetzen
for j=1:123
for k=1:3, pxx=px12+dp12*(k-2);
    z=ax12*sin(pxx+t*fx12)+z2; 
    yc=t.*(fGW+D1*dftt)+z; ys=y.*sin(yc); yc=y.*cos(yc); 
    yc=conv(ys,h4,'same').*stt+conv(yc,h4,'same').*ctt;
    erg(k,:)=[k-2,max(yc)];
end, p=polyfit(erg(1:3,1),erg(1:3,2),2); %max suchen
if p(1)<0, k=-p(2)/2/p(1); px12=px12+dp12*max(-19,min(19,k));
else
    if erg(1,2)<erg(3,2), px12=px12+4*dp12;
    else,  px12=px12-4*dp12;
    end
end
if abs(k)<1
    if dp12<1e-5 || max(erg(1:3,2))<1.4*min(erg(1:3,2)), break
    else, dp12=dp12/2;
    end
end
end

for j=1:123
for k=1:3, fxx=fx12+(k-2)*df12; 
    z=ax12*sin(px12+t*fxx)+z2; 
    yc=t.*(fGW+D1*dftt)+z; ys=y.*sin(yc); yc=y.*cos(yc); 
    yc=conv(ys,h4,'same').*stt+conv(yc,h4,'same').*ctt;
    erg(k,:)=[k-2,max(yc)];
end, p=polyfit(erg(1:3,1),erg(1:3,2),2); %max suchen
if p(1)<0, k=-p(2)/2/p(1); fx12=fx12+df12*max(-19,min(19,k));
else
    if erg(1,2)<erg(3,2), fx12=fx12+4*df12;
    else,  fx12=fx12-4*df12;
    end
end
if abs(k)<1
    if df12<2e-15 || max(erg(1:3,2))<1.4*min(erg(1:3,2)), break
    else, df12=df12/2;
    end
end
end

for j=1:123
for k=1:3, axx=ax12+da12*(k-2);
    z=axx*sin(px12+t*fx12)+z2; 
    yc=t.*(fGW+D1*dftt)+z; ys=y.*sin(yc); yc=y.*cos(yc); 
    yc=conv(ys,h4,'same').*stt+conv(yc,h4,'same').*ctt;
    erg(k,:)=[k-2,max(yc)];
end, p=polyfit(erg(1:3,1),erg(1:3,2),2); %max suchen
if p(1)<0, k=-p(2)/2/p(1); ax12=ax12+da12*max(-19,min(19,k));
else
    if erg(1,2)<erg(3,2), ax12=ax12+4*da12;
    else,  ax12=ax12-4*da12;
    end
end
if abs(k)<1
    if da12<1e-5 || max(erg(1:3,2))<1.4*min(erg(1:3,2)), break
    else, da12=da12/2;
    end
end
end %fx12.. ist abgeschlossen
z2=z2+ax12*sin(px12+t*fx12); %nach 3 Iterationen wieder einsetzen
end


if ax13~=0, z2=z2-ax13*sin(px13+t*fx13); %nach 3 Iterationen wieder einsetzen
for j=1:123
for k=1:3, pxx=px13+dp13*(k-2);
    z=ax13*sin(pxx+t*fx13)+z2; 
    yc=t.*(fGW+D1*dftt)+z; ys=y.*sin(yc); yc=y.*cos(yc); 
    yc=conv(ys,h4,'same').*stt+conv(yc,h4,'same').*ctt;
    erg(k,:)=[k-2,max(yc)];
end, p=polyfit(erg(1:3,1),erg(1:3,2),2); %max suchen
if p(1)<0, k=-p(2)/2/p(1); px13=px13+dp13*max(-19,min(19,k));
else
    if erg(1,2)<erg(3,2), px13=px13+4*dp13;
    else,  px13=px13-4*dp13;
    end
end
if abs(k)<1
    if dp13<1e-5 || max(erg(1:3,2))<1.4*min(erg(1:3,2)), break
    else, dp13=dp13/2;
    end
end
end

for j=1:123
for k=1:3, fxx=fx13+(k-2)*df13; 
    z=ax13*sin(px13+t*fxx)+z2; 
    yc=t.*(fGW+D1*dftt)+z; ys=y.*sin(yc); yc=y.*cos(yc); 
    yc=conv(ys,h4,'same').*stt+conv(yc,h4,'same').*ctt;
    erg(k,:)=[k-2,max(yc)];
end, p=polyfit(erg(1:3,1),erg(1:3,2),2); %max suchen
if p(1)<0, k=-p(2)/2/p(1); fx13=fx13+df13*max(-19,min(19,k));
else
    if erg(1,2)<erg(3,2), fx13=fx13+4*df13;
    else,  fx13=fx13-4*df13;
    end
end
if abs(k)<1
    if df13<2e-15 || max(erg(1:3,2))<1.4*min(erg(1:3,2)), break
    else, df13=df13/2;
    end
end
end

for j=1:123
for k=1:3, axx=ax13+da13*(k-2);
    z=axx*sin(px13+t*fx13)+z2; 
    yc=t.*(fGW+D1*dftt)+z; ys=y.*sin(yc); yc=y.*cos(yc); 
    yc=conv(ys,h4,'same').*stt+conv(yc,h4,'same').*ctt;
    erg(k,:)=[k-2,max(yc)];
end, p=polyfit(erg(1:3,1),erg(1:3,2),2); %max suchen
if p(1)<0, k=-p(2)/2/p(1); ax13=ax13+da13*max(-19,min(19,k));
else
    if erg(1,2)<erg(3,2), ax13=ax13+4*da13;
    else,  ax13=ax13-4*da13;
    end
end
if abs(k)<1
    if da13<1e-5 || max(erg(1:3,2))<1.4*min(erg(1:3,2)), break
    else, da13=da13/2;
    end
end
end %fx13.. ist abgeschlossen
z2=z2+ax13*sin(px13+t*fx13); %nach 3 Iterationen wieder einsetzen
end


if ax14~=0, z2=z2-ax14*sin(px14+t*fx14); %nach 3 Iterationen wieder einsetzen
for j=1:123
for k=1:3, pxx=px14+dp14*(k-2);
    z=ax14*sin(pxx+t*fx14)+z2; 
    yc=t.*(fGW+D1*dftt)+z; ys=y.*sin(yc); yc=y.*cos(yc); 
yc=conv(ys,h5,'same').*stt+conv(yc,h5,'same').*ctt;
    z0=T_zaehl_T(conv(yc,h1,'same'),Ts); %plot(z0(:,2))
    z0(:,3)=z0(:,2)-mean(z0(:,2));
    [p,~,mu]=polyfit(z0(:,1),z0(:,3),2);
    y1=polyval(p,z0(:,1),[],mu); %plot(z0(:,1),z0(:,3),z0(:,1),y1)
    erg(k,:)=[k-2,sum(y1.^2)];
end, p=polyfit(erg(1:3,1),erg(1:3,2),2); %min suchen
if p(1)>0, k=-p(2)/2/p(1); px14=px14+dp14*max(-19,min(19,k));
else
    if erg(1,2)>erg(3,2), px14=px14+4*dp14;
    else,  px14=px14-4*dp14;
    end
end
if abs(k)<1
    if dp14<1e-7 || max(erg(1:3,2))<1.4*min(erg(1:3,2)), break
    else, dp14=dp14/2;
    end
end
end

for j=1:123
for k=1:3, fxx=fx14+(k-2)*df14; 
    z=ax14*sin(px14+t*fxx)+z2; 
    yc=t.*(fGW+D1*dftt)+z; ys=y.*sin(yc); yc=y.*cos(yc); 
yc=conv(ys,h5,'same').*stt+conv(yc,h5,'same').*ctt;
    z0=T_zaehl_T(conv(yc,h1,'same'),Ts); %plot(z0(:,2))
    z0(:,3)=z0(:,2)-mean(z0(:,2));
    [p,~,mu]=polyfit(z0(:,1),z0(:,3),2);
    y1=polyval(p,z0(:,1),[],mu); %plot(z0(:,1),z0(:,3),z0(:,1),y1)
    erg(k,:)=[k-2,sum(y1.^2)];
end, p=polyfit(erg(1:3,1),erg(1:3,2),2); %min suchen
if p(1)>0, k=-p(2)/2/p(1); fx14=fx14+df14*max(-19,min(19,k));
else
    if erg(1,2)>erg(3,2), fx14=fx14+4*df14;
    else,  fx14=fx14-4*df14;
    end
end
if abs(k)<1
    if df14<2e-17 || max(erg(1:3,2))<1.4*min(erg(1:3,2)), break
    else, df14=df14/2;
    end
end
end

for j=1:123
for k=1:3, axx=ax14+da14*(k-2);
    z=axx*sin(px14+t*fx14)+z2; 
    yc=t.*(fGW+D1*dftt)+z; ys=y.*sin(yc); yc=y.*cos(yc); 
yc=conv(ys,h5,'same').*stt+conv(yc,h5,'same').*ctt;
    z0=T_zaehl_T(conv(yc,h1,'same'),Ts); %plot(z0(:,2))
    [p,~,mu]=polyfit(z0(:,1),z0(:,2),2);
    y1=polyval(p,z0(:,1),[],mu); %plot(z0(:,1),z0(:,3),z0(:,1),y1)
    erg(k,:)=[k-2,sum(y1.^2)];
end, p=polyfit(erg(1:3,1),erg(1:3,2),2); %min suchen
if p(1)>0, k=-p(2)/2/p(1); ax14=ax14+da14*max(-19,min(19,k));
else
    if erg(1,2)>erg(3,2), ax14=ax14+4*da14;
    else,  ax14=ax14-4*da14;
    end
end
if abs(k)<1
    if da14<1e-7 || max(erg(1:3,2))<1.4*min(erg(1:3,2)), break
    else, da14=da14/2;
    end
end
end %fx14.. ist abgeschlossen
z2=z2+ax14*sin(px14+t*fx14); %nach 3 Iterationen wieder einsetzen
end

if ax15~=0, z2=z2-ax15*sin(px15+t*fx15); %nach 3 Iterationen wieder einsetzen
for j=1:123
for k=1:3, pxx=px15+dp15*(k-2);
    z=ax15*sin(pxx+t*fx15)+z2; 
    yc=t.*(fGW+D1*dftt)+z; ys=y.*sin(yc); yc=y.*cos(yc); 
yc=conv(ys,h5,'same').*stt+conv(yc,h5,'same').*ctt;
    z0=T_zaehl_T(conv(yc,h1,'same'),Ts); %plot(z0(:,2))
    z0(:,3)=z0(:,2)-mean(z0(:,2));
    [p,~,mu]=polyfit(z0(:,1),z0(:,3),2);
    y1=polyval(p,z0(:,1),[],mu); %plot(z0(:,1),z0(:,3),z0(:,1),y1)
    erg(k,:)=[k-2,sum(y1.^2)];
end, p=polyfit(erg(1:3,1),erg(1:3,2),2); %min suchen
if p(1)>0, k=-p(2)/2/p(1); px15=px15+dp15*max(-19,min(19,k));
else
    if erg(1,2)>erg(3,2), px15=px15+4*dp15;
    else,  px15=px15-4*dp15;
    end
end
if abs(k)<1
    if dp15<1e-7 || max(erg(1:3,2))<1.4*min(erg(1:3,2)), break
    else, dp15=dp15/2;
    end
end
end

for j=1:123
for k=1:3, fxx=fx15+(k-2)*df15; 
    z=ax15*sin(px15+t*fxx)+z2; 
    yc=t.*(fGW+D1*dftt)+z; ys=y.*sin(yc); yc=y.*cos(yc); 
yc=conv(ys,h5,'same').*stt+conv(yc,h5,'same').*ctt;
    z0=T_zaehl_T(conv(yc,h1,'same'),Ts); %plot(z0(:,2))
    z0(:,3)=z0(:,2)-mean(z0(:,2));
    [p,~,mu]=polyfit(z0(:,1),z0(:,3),2);
    y1=polyval(p,z0(:,1),[],mu); %plot(z0(:,1),z0(:,3),z0(:,1),y1)
    erg(k,:)=[k-2,sum(y1.^2)];
end, p=polyfit(erg(1:3,1),erg(1:3,2),2); %min suchen
if p(1)>0, k=-p(2)/2/p(1); fx15=fx15+df15*max(-19,min(19,k));
else
    if erg(1,2)>erg(3,2), fx15=fx15+4*df15;
    else,  fx15=fx15-4*df15;
    end
end
if abs(k)<1
    if df15<2e-17 || max(erg(1:3,2))<1.4*min(erg(1:3,2)), break
    else, df15=df15/2;
    end
end
end

for j=1:123
for k=1:3, axx=ax15+da15*(k-2);
    z=axx*sin(px15+t*fx15)+z2; 
    yc=t.*(fGW+D1*dftt)+z; ys=y.*sin(yc); yc=y.*cos(yc); 
yc=conv(ys,h5,'same').*stt+conv(yc,h5,'same').*ctt;
    z0=T_zaehl_T(conv(yc,h1,'same'),Ts); %plot(z0(:,2))
    z0(:,3)=z0(:,2)-mean(z0(:,2));
    [p,~,mu]=polyfit(z0(:,1),z0(:,3),2);
    y1=polyval(p,z0(:,1),[],mu); %plot(z0(:,1),z0(:,3),z0(:,1),y1)
    erg(k,:)=[k-2,sum(y1.^2)];
end, p=polyfit(erg(1:3,1),erg(1:3,2),2); %min suchen
if p(1)>0, k=-p(2)/2/p(1); ax15=ax15+da15*max(-19,min(19,k));
else
    if erg(1,2)>erg(3,2), ax15=ax15+4*da15;
    else,  ax15=ax15-4*da15;
    end
end
if abs(k)<1
    if da15<1e-7 || max(erg(1:3,2))<1.4*min(erg(1:3,2)), break
    else, da15=da15/2;
    end
end
end %fx15.. ist abgeschlossen
z2=z2+ax15*sin(px15+t*fx15); %nach 3 Iterationen wieder einsetzen
end

%k=mean(z0(:,2))-Ts*5-5.5e-5; %Abweichung von 1/fZF
k=mean(z0(:,2))-Ts*5+ww; %Abweichung von 1/fZF
R(k3,1:31)=[fGW,D1,ax,px,fx,ax2,px2,fx2,ax3,px3,fx3,axJ,pxJ,fxJ,jjj,ax4,px4,fx4,ax5,px5,fx5,ax6,px6,fx6,ax7,px7,fx7,ax8,px8,fx8,k];
R(k3,32:52)=[ax9,px9,fx9,ax10,px10,fx10,ax11,px11,fx11,ax12,px12,fx12,ax13,px13,fx13,ax14,px14,fx14,ax15,px15,fx15];
if k3>k3alt+10 && mean(R(k3-5:k3,15))<Limi %ruhige Iteration? 
    plot(R(max(1,k3-240):k3,15)), drawnow%, disp(k)
    k3alt=k3; %Position speichern
    
for j=3:4 %D1 bestimmen
    for k=1:2 %yc soll bei 0 beginnen
        yc=t.*(fGW+D1*dftt)+z2+phas; ys=y.*sin(yc); yc=y.*cos(yc); 
        yc=conv(ys,h4,'same').*stt+conv(yc,h4,'same').*ctt; 
        erg(k,1:2)=[phas,yc(1)]; phas=phas+0.04;
    end
    phas=erg(1,1)-erg(1,2)*(erg(2,1)-erg(1,1))/(erg(2,2)-erg(1,2)); %genauer 0-Durchgang 
    yc=t.*(fGW+D1*dftt)+z2+phas; ys=y.*sin(yc); yc=y.*cos(yc); 
    yc=conv(ys,h4,'same').*stt+conv(yc,h4,'same').*ctt; %plot(yc(1:16))
    k1=1;
    for k=1:10:length(yc)
        z0(k1,1)=yc(k); k1=k1+1;
    end%, plot(z0(1:k1-1,1))
    erg(j,:)=[D1,z0(k1-1,1)-z0(1,1)]; D1=D1*1.002;
end
if erg(3,2)*erg(4,2)<0
D1=erg(3,1)-erg(3,2)*(erg(4,1)-erg(3,1))/(erg(4,2)-erg(3,2)); %genauer 0-Durchgang
else, D1=erg(3,1)+erg(3,2)*200;
end
    
    if k3>20, ww=ww-mean(R(max(10,k3-30):k3,31))/50; disp(ww),end
    %if abs(ww)<1e-16, beep, break, end
    %fGW=fGW+k*ww; %jjj=0; %Frequenz korrigieren
elseif k3>k3alt+20
    k3alt=k3; %Limi=Limi+0.2; 
end
end, R=R(1:k3-1,:); return


z0(:,3)=z0(:,2)-mean(z0(:,2));
[p,~,mu]=polyfit(z0(:,1),z0(:,3),4);
z0(:,4)=polyval(p,z0(:,1),[],mu); 
plot(z0(:,1),z0(:,3),z0(:,1),z0(:,4)), return


plot(z0(:,1),z0(:,3))
[b,a] = cheby1(6,0.01,0.2); 
plot(filtfilt(b,a,z0(:,2)-mean(z0(:,2))))

















if 1==100000000000000000, z2=ax2*sin(px2+t*fx2)+z;
while 1==1 %D1
for k=1:3, j=D1+dD*(k-2);
    yc=t.*(fGW+j*dftt)+z2; ys=y.*sin(yc); yc=y.*cos(yc); 
    yc=filtfilt(b,a,ys).*stt+filtfilt(b,a,yc).*ctt;
    ys=sync(conv(yc,h1,'same'),fZF,Ts); %plot(ys)
    y1=abs(hilbert(ys)); y1=y1(100:end-100); yc=1:length(y1);
    [p,~,mu]=polyfit(yc,y1,1); y2=y1-polyval(p,yc',[],mu); 
    [p,~,mu]=polyfit(yc,y2,1); y2=polyval(p,yc',[],mu); 
    erg(k,:)=[k-2,sum(y2.^2)];
end, p=polyfit(erg(1:3,1),erg(1:3,2),2);
if p(1)>0, k=-p(2)/2/p(1); D1=D1+dD*max(-19,min(19,k));
else, break
end
if abs(k)<0.4
    if dD<2e-6 || max(erg(1:3,2))<1.1*min(erg(1:3,2)), break
    else, dD=dD/2;
    end
end
end
end


%----------------------------------------
%------Frequenz darstellen-------------
y1=sin(t.*(fZF+D1*dftt));
z0=T_zaehl_T(y1,Ts); k=mean(z0(:,2)); 
h1=window_sinc_filter(2e3,0.2e-6,1,1/Ts,'low','blackman');
z0(:,2)=(k+conv(z0(:,2)-k,h1,'same'))/1e6; %Einheit = µHz
z0(:,1)=z0(:,1)/3600/24/tJ; 
plot(z0(:,1),z0(:,2)), title('rel. Schwingungsdauer des Oszi')
xlabel('Time (Years)'), ylabel('T (µs)')
for k=1:size(z0,1) %daraus Frequenz berechnen
    z0(k,3)=0.5/z0(k,2); %┬ÁHz
end, plot(z0(:,1),z0(:,3))
xlabel('Time (Years)'), ylabel('Frequency (┬ÁHz)'),title('Oszillatorfrequenz')
[~,erg]=findpeaks(z0(:,3),z0(:,1)); erg=erg-floor(erg);
k=mod(356.256*(0.5+mean(erg(2:length(erg)-2))),356.256); disp([k,k/30.41])

%-------------------------------------------------------

function erg=T_zaehl_T(y,Ts) %mit Zeitzuordnung
yc=y.*cat(1,0,y(1:end-1)); erg=zeros(300,2); %1=Zeit, 2=Zeitdiff
ys=(yc<0); %Vorzeichenwechsel
nn=1; j=1; while ys(j)==0, j=j+1; end %1. VZ-Wechsel
b=j-1+y(j-1)/(y(j-1)-y(j)); %genauer 0-Durchgang 
while j<length(ys)-30
    j=j+1; while ys(j)==0, j=j+1; end %nächster Wechsel
    a=j-1+y(j-1)/(y(j-1)-y(j)); %genauer 0-Durchgang 
    erg(nn,1)=Ts*(a+b)/2;
    erg(nn,2)=Ts*(a-b); nn=nn+1; b=a;
end 
erg=erg(1:nn-1,:);
end

function y=sync(yGW, f, Ts)
a=2*pi*Ts*f; b=sin(a); a=cos(a); L=length(yGW);
y=zeros(L,1); y1=y;
if yGW(1)>0, y(1)=50; %willkürlicher Startwert
else, y(1)=-50;
end
for j=2:L %vorwärts
y(j)=y(j-1)*a+y1(j-1)*b+yGW(j); %Add-Theorem
y1(j)=y1(j-1)*a-y(j-1)*b;
end%, figure(1),plot(y), title('Coherent Demodulation')
end

%frequenztransformierendes Bandfilter um f_xx
%mit wählbarer Flankensteilheit (n)
function ssb = filtTrans_n_Hz(y,Ts,n,f_ein,f_aus,B)
%f_xx und B(andbreite) in Hz
%n ist die Länge von h1
%wenn B<0 wird cheby-Filter verwendet
L=length(y); j=(1:L)'*2*pi*f_ein*Ts;
ys=y.*sin(j); yc=y.*cos(j);
if B>0
    h1=window_sinc_filter(n,B,1,1/Ts,'low','blackman');
    ys=conv(ys,h1,'same'); yc=conv(yc,h1,'same');
else
    [b,a] = cheby1(6,0.01,-B); %0.006
    yc=filtfilt(b,a,yc); ys=filtfilt(b,a,ys);
end
if f_ein~=f_aus
    j=(1:L)'*2*pi*f_aus*Ts; %verschobene MittenFrequenz
end
ys=ys.*sin(j); yc=yc.*cos(j);
ssb=2*(yc+ys).*wei(100,L); %Ts ist unverändert, B wie gewählt
%[erg f]=zeig_sp2(ssb,-Ts,'log');
end
%
function [y, Ts]=hol_DWD(a) %files hergestellt mit prepD2000d.m
%a=1 -> qxxx.mat für 2000..2009; 
%a=2 -> pxxx.mat für 2010..2019
oldpath=path;
path(oldpath,'C:\Users\herbe\Documents\Matlab\ScoX1\')
if a==1 || a==3
load 'q00183.mat' y, yc=y;
load 'q00282.mat' y, yc=yc+y;
load 'q00298.mat' y, yc=yc+y;
load 'q00303.mat' y, yc=yc+y;
load 'q00399.mat' y, yc=yc+y;
load 'q00427.mat' y, yc=yc+y;
load 'q00430.mat' y, yc=yc+y;
load 'q00656.mat' y, yc=yc+y;
load 'q00691.mat' y, yc=yc+y;
load 'q00722.mat' y, yc=yc+y;
load 'q00867.mat' y, yc=yc+y;
load 'q00880.mat' y, yc=yc+y;
load 'q01262.mat' y, yc=yc+y;
load 'q01270.mat' y, yc=yc+y;
load 'q01346.mat' y, yc=yc+y;
load 'q01420.mat' y, yc=yc+y;
load 'q01443.mat' y, yc=yc+y;
load 'q01468.mat' y, yc=yc+y;
load 'q01580.mat' y, yc=yc+y;
load 'q01639.mat' y, yc=yc+y;
load 'q01684.mat' y, yc=yc+y;
load 'q01691.mat' y, yc=yc+y;
load 'q01766.mat' y, yc=yc+y;
load 'q01957.mat' y, yc=yc+y;
load 'q01975.mat' y, yc=yc+y;
load 'q02115.mat' y, yc=yc+y;
load 'q02014.mat' y, yc=yc+y;
load 'q02171.mat' y, yc=yc+y;
load 'q02532.mat' y, yc=yc+y;
load 'q02559.mat' y, yc=yc+y;
load 'q02638.mat' y, yc=yc+y;
load 'q02667.mat' y, yc=yc+y;
load 'q02812.mat' y, yc=yc+y;
load 'q02925.mat' y, yc=yc+y;
load 'q03015.mat' y, yc=yc+y;
load 'q03023.mat' y, yc=yc+y;
load 'q03028.mat' y, yc=yc+y;
load 'q03032.mat' y, yc=yc+y;
load 'q03086.mat' y, yc=yc+y;
load 'q03093.mat' y, yc=yc+y;
load 'q03196.mat' y, yc=yc+y;
load 'q03287.mat' y, yc=yc+y;
load 'q03379.mat' y, yc=yc+y;
load 'q03552.mat' y, yc=yc+y;
load 'q03631.mat' y, yc=yc+y;
load 'q03660.mat' y, yc=yc+y;
load 'q03730.mat' y, yc=yc+y;
load 'q03761.mat' y, yc=yc+y;
load 'q03791.mat' y, yc=yc+y;
load 'q03811.mat' y, yc=yc+y;
load 'q03987.mat' y, yc=yc+y;
load 'q04336.mat' y, yc=yc+y;
load 'q04625.mat' y, yc=yc+y;
load 'q04642.mat' y, yc=yc+y;
load 'q04745.mat' y, yc=yc+y;
load 'q04887.mat' y, yc=yc+y;
load 'q04928.mat' y, yc=yc+y;
load 'q04931.mat' y, yc=yc+y;
load 'q05029.mat' y, yc=yc+y;
load 'q05142.mat' y, yc=yc+y;
load 'q05155.mat' y, yc=yc+y;
load 'q05426.mat' y, yc=yc+y;
load 'q05629.mat' y, yc=yc+y;
load 'q05705.mat' y, yc=yc+y;
load 'q05906.mat' y, yc=yc+y;
y=yc; 
end
if a==2 || a==3
load 'p00183.mat' y, ys=y;
load 'p00198.mat' y, ys=ys+y;
load 'p00232.mat' y, ys=ys+y;
load 'p00282.mat' y, ys=ys+y;
load 'p00303.mat' y, ys=ys+y;
load 'p00427.mat' y, ys=ys+y;
load 'p00430.mat' y, ys=ys+y;
load 'p00433.mat' y, ys=ys+y;
load 'p00591.mat' y, ys=ys+y;
load 'p00656.mat' y, ys=ys+y;
load 'p00691.mat' y, ys=ys+y;
load 'p00853.mat' y, ys=ys+y;
load 'p00880.mat' y, ys=ys+y;
load 'p00953.mat' y, ys=ys+y;
load 'p01001.mat' y, ys=ys+y;
load 'p01048.mat' y, ys=ys+y;
load 'p01078.mat' y, ys=ys+y;
load 'p01200.mat' y, ys=ys+y;
load 'p01262.mat' y, ys=ys+y;
load 'p01270.mat' y, ys=ys+y;
load 'p01303.mat' y, ys=ys+y;
load 'p01420.mat' y, ys=ys+y;
load 'p01443.mat' y, ys=ys+y;
load 'p01468.mat' y, ys=ys+y;
load 'p01550.mat' y, ys=ys+y;
load 'p01580.mat' y, ys=ys+y;
load 'p01587.mat' y, ys=ys+y;
load 'p01605.mat' y, ys=ys+y;
load 'p01612.mat' y, ys=ys+y;
load 'p01639.mat' y, ys=ys+y;
load 'p01684.mat' y, ys=ys+y;
load 'p01694.mat' y, ys=ys+y;
load 'p01757.mat' y, ys=ys+y;
load 'p01869.mat' y, ys=ys+y;
load 'p01975.mat' y, ys=ys+y;
load 'p02014.mat' y, ys=ys+y;
load 'p02023.mat' y, ys=ys+y;
load 'p02044.mat' y, ys=ys+y;
load 'p02115.mat' y, ys=ys+y;
load 'p02171.mat' y, ys=ys+y;
load 'p02261.mat' y, ys=ys+y;
load 'p02290.mat' y, ys=ys+y;
load 'p02429.mat' y, ys=ys+y;
load 'p02483.mat' y, ys=ys+y;
load 'p02597.mat' y, ys=ys+y;
load 'p02601.mat' y, ys=ys+y;
load 'p02638.mat' y, ys=ys+y;
load 'p02667.mat' y, ys=ys+y;
load 'p02712.mat' y, ys=ys+y;
load 'p02794.mat' y, ys=ys+y;
load 'p03287.mat' y, ys=ys+y;
load 'p04466.mat' y, ys=ys+y;
y=ys;
end
if a==3, y=cat(1,yc,ys); %Jahre: q=2000..2009; p=2010..2019;
end, y=y-mean(y); Ts=3600; %Sekunden
end

