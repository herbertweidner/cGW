%Luftdruck Resonanzen 
%[y, Ts]=hol_DWD(3); %1=2000..2009; 2=2010..2019
load yDWD, Ts=3600;
L=length(y); figure (1), zeig_sp2(y,-Ts,'log');
fGW0=4.47821445E-006; %HW95=fGW0*360*3600; optimiert

tD=23.93447192; fD=1/(tD*3600); 
tJ=365.25636042; fJ=1/tJ/24/3600; 
[sp,f]=zeig_sp2(y,Ts,-200000.3);
k=1;while f(k)<fGW0,k=k+1;end,k=k-3:k+3; plot(1e6*f(k),sp(k))
t=(1:Ts:21/fJ)'; t=2*pi*t(1:length(y));

y0=y.*(1.5+sin(t*fD+1.7)); %++++++++++++++++++++++++++++

Offset=3.9e-6; fGW=fGW0-Offset; korr=16e-11;
y= filtTrans_n_Hz(y0,Ts,3e4,fGW0,fGW,0.45e-6); %bis 2*200nHz wegen Mod-Index
y=decimate(decimate(y,10),10); Ts=100*Ts; zeig_sp2(y,-Ts,'log'); 
[sp,f]=zeig_sp2(y,Ts,-300000.3);k=1;while f(k)<fGW,k=k+1;end,q1=k-33:k+44;
plot(1e9*f(q1),sp(q1)), title('verschobene fGW')
t=(1:Ts:21/fJ)'; tt=1:(length(y)); fZF=1/Ts/10; %~1e-6;
    dftt=1e-20*t(tt); t=2*pi*t(tt); %41 Jahre FM darstellen
    tt=(tt)'*2*pi*fZF*Ts; %fixe ZF wegen variablem Oszillator
    stt=sin(tt); ctt=cos(tt);
j=0.4e-9; h1=window_sinc_filter(1e5,fZF-j,fZF+j,1/Ts,'bandpass','Blackman');
h4=window_sinc_filter(1e5,0.5e-9,1,1/Ts,'low','blackman');
%h4=window_sinc_filter(3e4,1.0e-9,1,1/Ts,'low','blackman');
h5=window_sinc_filter(1e4,10e-9,1,1/Ts,'low','blackman');
%j=5*0.08e-9;
%h6=window_sinc_filter(1e5,fZF+fJ-j,fZF+fJ+j,1/Ts,'bandpass','Blackman');
%h7=window_sinc_filter(5e5,fZF-fJ-j,fZF-fJ+j,1/Ts,'bandpass','Blackman');
%{
ax=0.2; px=3; fx=38.9735e-9; drift=100.85; % 
ax2=0; px2=3; fx2=59.375e-9; % C
ax3=0.2; px3=3; fx3=4e-10; %Frequenzzähler M
ax4=0; px4=3; fx4=4.29075e-10; % Frequenzzähler L 
ax5=0; px5=3; fx5=92.67406e-9; % D
ax6=0; px6=3; fx6=23.44711e-9; % F
ax7=0; px7=3; fx7=3.5202e-10; % Frequenzzähler  Jahre
ax8=0; px8=3; fx8=8.68867e-9; % 
ax9=0; px9=3; fx9=2.62448e-9; % E
ax10=0; px10=3; fx10=6.12776e-9; % Jahre
ax11=0; px11=3; fx11=0.98e-9; % H
ax12=0; px12=3; fx12=4.8e-9; % J
ax13=0; px13=3; fx13=1.24e-9; % K
ax14=0; px14=3; fx14=1.39276e-10; % Frequenzzähler
ax15=0; px15=3; fx15=1.20464e-10; % Frequenzzähler
axJ=1; pxJ=3; fxJ=31.68754e-9; %31.69 nHz
%}
ax=5.3613; px=2.30296; fx=38.50585e-9; drift=1100.875; 
ax2=4.966; px2=0.6514; fx2=14.1995e-9; % C
ax3=2.11032; px3=3.13838; fx3=10.6706e-10; %Frequenzzähler M
ax4=0.0202; px4=6.812; fx4=-5.51e-10; % Frequenzzähler L 
ax5=1.4995; px5=4.199; fx5=95.033e-9; % D
ax6=3.798; px6=5.704; fx6=3.3864e-9; % F
ax7=0.8652; px7=3.985; fx7=4.5032e-10; % Frequenzzähler  Jahre
ax8=2.2425; px8=0.732; fx8=74.2797e-9; % 
ax9=0; px9=0; fx9=63e-9; % 
ax10=0; px10=3; fx10=6.12776e-9; % Jahre
ax11=0; px11=3; fx11=0.98e-9; % H
ax12=0; px12=3; fx12=4.8e-9; % J
ax13=0; px13=3; fx13=1.24e-9; % K
ax14=0; px14=3; fx14=1.39276e-10; % Frequenzzähler
ax15=0; px15=3; fx15=1.20464e-10; % Frequenzzähler
axJ=5.2023; pxJ=4.9166; fxJ=31.67443e-9; %31.69 nHz

dD=1;
R=zeros(8695,52); erg=zeros(3,2); Limi=-11; k3alt=1; j=0;
for k3=1:size(R,1)
z2=ax*sin(px+t*fx)+ax2*sin(px2+t*fx2)+ax3*sin(px3+t*fx3)+ax4*sin(px4+t*fx4)+ax5*sin(px5+t*fx5)+ax6*sin(px6+t*fx6)+ax7*sin(px7+t*fx7)+ax8*sin(px8+t*fx8);
z2=z2+axJ*sin(pxJ+t*fxJ)+ax9*sin(px9+t*fx9)+ax10*sin(px10+t*fx10)+ax11*sin(px11+t*fx11)+ax12*sin(px12+t*fx12)+ax13*sin(px13+t*fx13)+ax14*sin(px14+t*fx14);
z2=z2+ax15*sin(px15+t*fx15);

for k=1:2, k1=drift+(k-2)*dD; %drift bestimmen
    yc=t.*(fGW+k1*dftt)+z2; ys=y.*sin(yc); yc=y.*cos(yc); 
    yi=conv(ys,h5,'same').*stt+conv(yc,h5,'same').*ctt;
    z0=T_zaehl_T(conv(yi,h1,'same'),Ts); %plot(z0(:,2))
    [p]=polyfit(z0(:,1),z0(:,2),1);
    erg(k,:)=[k1,p(1)];
end
drift=drift-erg(2,2)*dD/(erg(2,2)-erg(1,2)); %genauer 0-Durchgang 
%drift=min(350,max(drift,-150));
%y1=conv(ys,h1+h6+h7,'same').*stt+conv(yc,h1+h6+h7,'same').*ctt;
%zeig_sp2(y1-mean(y1),Ts,20000.3);

%[A,F,P,z2]=ITER(A,F,P,t,fGW,drift,y,z2,h4,dftt,stt,ctt);
[axJ,fxJ,pxJ,z2]=ITER(axJ,fxJ,pxJ,t,fGW,drift,y,z2,h4,dftt,stt,ctt);
[ax,fx,px,z2]=ITER(ax,fx,px,t,fGW,drift,y,z2,h4,dftt,stt,ctt);
[ax2,fx2,px2,z2]=ITER(ax2,fx2,px2,t,fGW,drift,y,z2,h4,dftt,stt,ctt);
[ax5,fx5,px5,z2]=ITER(ax5,fx5,px5,t,fGW,drift,y,z2,h4,dftt,stt,ctt);
[ax6,fx6,px6,z2]=ITER(ax6,fx6,px6,t,fGW,drift,y,z2,h4,dftt,stt,ctt);
[ax8,fx8,px8,z2]=ITER(ax8,fx8,px8,t,fGW,drift,y,z2,h4,dftt,stt,ctt);
[ax9,fx9,px9,z2]=ITER(ax9,fx9,px9,t,fGW,drift,y,z2,h4,dftt,stt,ctt);

%[A,F,P,z2]=ITER_z(A,F,P,t,fGW,drift,y,z2,h4,dftt,stt,ctt);
[ax3,fx3,px3,z2,j]=ITER_z(ax3,fx3,px3,t,fGW,drift,y,z2,h1,h5,Ts,dftt,stt,ctt,j);
[ax4,fx4,px4,z2,j]=ITER_z(ax4,fx4,px4,t,fGW,drift,y,z2,h1,h5,Ts,dftt,stt,ctt,j);
[ax7,fx7,px7,z2,j]=ITER_z(ax7,fx7,px7,t,fGW,drift,y,z2,h1,h5,Ts,dftt,stt,ctt,j);

k=mean(z0(:,2))-Ts*5; %disp(k) %Abweichung von 1/fZF
R(k3,1:31)=[fGW+Offset,drift,ax,px,fx,ax2,px2,fx2,ax3,px3,fx3,axJ,pxJ,fxJ,j,ax4,px4,fx4,ax5,px5,fx5,ax6,px6,fx6,ax7,px7,fx7,ax8,px8,fx8,k];
R(k3,32:52)=[ax9,px9,fx9,ax10,px10,fx10,ax11,px11,fx11,ax12,px12,fx12,ax13,px13,fx13,ax14,px14,fx14,ax15,px15,fx15];
if k3>k3alt+15 && mean(R(k3-5:k3,15))<Limi %ruhige Iteration? 
    plot(R(max(1,k3-240):k3,15)), drawnow, k3alt=k3; %Position speichern
    fGW=fGW+k*korr; %j=0; %Frequenz korrigieren
    disp(k)
elseif k3>k3alt+20
    Limi=Limi+0.2; k3alt=k3;
end
end, R=R(1:k3-1,:); %return
fGW0=fGW+Offset; disp(fGW0);

%PSD berechnen
h6=window_sinc_filter(2e3,240e-9,1,1/Ts,'low','blackman');
y1=conv(ys,h6,'same').*stt+conv(yc,h6,'same').*ctt;
yi=conv(ys,h6,'same').*stt-conv(yc,h6,'same').*ctt; %invers
[sp,~]=zeig_sp2(y1,Ts,2000.4); [spi,f]=zeig_sp2(yi,Ts,2000.4);
plot(f,sp,f,spi)
k=16*1024;[sp,f]= pwelch(cat(1,y1,zeros(3e5,1)),rectwin(k),k/8,16*k,1/Ts);
plot(1e9*f,sp/1e5), xlabel('Frequency (nHz)')
k3=1/Ts/length(sp)/2*sum(sp(26100:26320)-0.03); k=f(26100)-f(26320); %PSD
return

plot(z0(:,1),z0(:,3))
[b,a] = cheby1(6,0.01,0.2); 
plot(filtfilt(b,a,z0(:,2)-mean(z0(:,2))))

















if 1==100000000000000000, z2=ax2*sin(px2+t*fx2)+z;
while 1==1 %drift
for k=1:3, j=drift+dD*(k-2);
    yc=t.*(fGW+j*dftt)+z2; ys=y.*sin(yc); yc=y.*cos(yc); 
    yc=filtfilt(b,a,ys).*stt+filtfilt(b,a,yc).*ctt;
    ys=sync(conv(yc,h1,'same'),fZF,Ts); %plot(ys)
    y1=abs(hilbert(ys)); y1=y1(100:end-100); yc=1:length(y1);
    [p,~,mu]=polyfit(yc,y1,1); y2=y1-polyval(p,yc',[],mu); 
    [p,~,mu]=polyfit(yc,y2,1); y2=polyval(p,yc',[],mu); 
    erg(k,:)=[k-2,sum(y2.^2)];
end, p=polyfit(erg(:,1),erg(:,2),2);
if p(1)>0, k=-p(2)/2/p(1); drift=drift+dD*max(-19,min(19,k));
else, break
end
if abs(k)<0.4
    if dD<2e-6 || max(erg(:,2))<1.1*min(erg(:,2)), break
    else, dD=dD/2;
    end
end
end
end


%----------------------------------------
%------Frequenz darstellen-------------
y1=sin(t.*(fZF+drift*dftt));
z0=T_zaehl_T(y1,Ts); k=mean(z0(:,2)); 
h1=window_sinc_filter(5e3,240e-9,1,1/Ts,'low','blackman');
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

%Modulationen mit f(mod)>10 nHz sollte man die Funktion ITER suchen
%Langsamere Modulationen mit ITER_z
S
function [A,F,P,z2]=ITER(A,F,P,t,fGW,drift,y,z2,h4,dftt,stt,ctt)
if A==0, return, end 
dp=0.1; df=1e-12; da=1e-3; erg=zeros(3,2); 

z2=z2-A*sin(P+t*F); %nach 3 Iterationen wieder einsetzen
for j=1:123
for k=1:3, pxx=P+dp*(k-2);
    z=A*sin(pxx+t*F)+z2; 
    yc=t.*(fGW+drift*dftt)+z; ys=y.*sin(yc); yc=y.*cos(yc); 
    yi=conv(ys,h4,'same').*stt+conv(yc,h4,'same').*ctt;
    erg(k,:)=[k-2,max(yi)];
end, p=polyfit(erg(:,1),erg(:,2),2); %max suchen
if p(1)<0, k=-p(2)/2/p(1); P=P+dp*max(-19,min(19,k));
else
    if erg(1,2)<erg(3,2), P=P+4*dp;
    else,  P=P-4*dp;
    end
end
if abs(k)<1
    if dp<1e-5 || max(erg(:,2))<1.4*min(erg(:,2)), break
    else, dp=dp/2;
    end
end
end


for j=1:123
for k=1:3, fxx=F+(k-2)*df; 
    z=A*sin(P+t*fxx)+z2; 
    yc=t.*(fGW+drift*dftt)+z; ys=y.*sin(yc); yc=y.*cos(yc); 
    yi=conv(ys,h4,'same').*stt+conv(yc,h4,'same').*ctt;
    erg(k,:)=[k-2,max(yi)];
end, p=polyfit(erg(:,1),erg(:,2),2); %max suchen
if p(1)<0, k=-p(2)/2/p(1); F=F+df*max(-19,min(19,k));
else
    if erg(1,2)<erg(3,2), F=F+4*df;
    else, F=F-4*df;
    end
end
if abs(k)<1
    if df<2e-18 || max(erg(:,2))<1.4*min(erg(:,2)), break
    else, df=df/2;
    end
end
end

for j=1:123
for k=1:3, axx=A+da*(k-2);
    z=axx*sin(P+t*F)+z2; 
    yc=t.*(fGW+drift*dftt)+z; ys=y.*sin(yc); yc=y.*cos(yc); 
    yi=conv(ys,h4,'same').*stt+conv(yc,h4,'same').*ctt;
    erg(k,:)=[k-2,max(yi)];
end, p=polyfit(erg(:,1),erg(:,2),2); %max suchen
if p(1)<0, k=-p(2)/2/p(1); A=A+da*max(-19,min(19,k));
else
    if erg(1,2)<erg(3,2), A=A+4*da;
    else, A=A-4*da;
    end
end
if abs(k)<1
    if da<1e-6 || max(erg(:,2))<1.4*min(erg(:,2)), break
    else, da=da/2;
    end
end
end %F.. ist abgeschlossen
z2=z2+A*sin(P+t*F); %nach 3 Iterationen wieder einsetzen
end


function [A,F,P,z2,j]=ITER_z(A,F,P,t,fGW,drift,y,z2,h1,h5,Ts,dftt,stt,ctt,j)
if A==0, return, end 
dp=1e-4; df=1e-14; da=1e-5; erg=zeros(3,2); 

z2=z2-A*sin(P+t*F); %nach 3 Iterationen wieder einsetzen
for j=1:123
for k=1:3, pxx=P+dp*(k-2);
    z=A*sin(pxx+t*F)+z2; 
    yc=t.*(fGW+drift*dftt)+z; ys=y.*sin(yc); yc=y.*cos(yc); 
    yi=conv(ys,h5,'same').*stt+conv(yc,h5,'same').*ctt;
    z0=T_zaehl_T(conv(yi,h1,'same'),Ts); %plot(z0(:,2))
    z0(:,3)=z0(:,2)-mean(z0(:,2));
    [p,~,mu]=polyfit(z0(:,1),z0(:,3),2);
    z0(:,3)=polyval(p,z0(:,1),[],mu); %plot(z0(:,1),z0(:,3),z0(:,1),y1)
    erg(k,:)=[k-2,sum(z0(:,3).^2)];
end, p=polyfit(erg(:,1),erg(:,2),2); %min suchen
if p(1)>0, k=-p(2)/2/p(1); P=P+dp*max(-19,min(19,k));
else
    if erg(1,2)>erg(3,2), P=P+4*dp;
    else,  P=P-4*dp;
    end
end
if abs(k)<1
    if dp<1e-7 || max(erg(:,2))<1.4*min(erg(:,2)), break
    else, dp=dp/2;
    end
end
end

for j=1:123
for k=1:3, fxx=F+(k-2)*df; 
    z=A*sin(P+t*fxx)+z2; 
    yc=t.*(fGW+drift*dftt)+z; ys=y.*sin(yc); yc=y.*cos(yc); 
    yi=conv(ys,h5,'same').*stt+conv(yc,h5,'same').*ctt;
    z0=T_zaehl_T(conv(yi,h1,'same'),Ts); %plot(z0(:,2))
    z0(:,3)=z0(:,2)-mean(z0(:,2));
    [p,~,mu]=polyfit(z0(:,1),z0(:,3),2);
    z0(:,3)=polyval(p,z0(:,1),[],mu); %plot(z0(:,1),z0(:,3),z0(:,1),y1)
    erg(k,:)=[k-2,sum(z0(:,3).^2)];
end, p=polyfit(erg(:,1),erg(:,2),2); %min suchen
if p(1)>0, k=-p(2)/2/p(1); F=F+df*max(-19,min(19,k));
else
    if erg(1,2)>erg(3,2), F=F+4*df;
    else,  F=F-4*df;
    end
end
if abs(k)<1
    if df<1e-17 || max(erg(:,2))<1.4*min(erg(:,2)), break
    else, df=df/2;
    end
end
end

for j=1:123
for k=1:3, axx=A+da*(k-2);
    z=axx*sin(P+t*F)+z2; 
    yc=t.*(fGW+drift*dftt)+z; ys=y.*sin(yc); yc=y.*cos(yc); 
    yi=conv(ys,h5,'same').*stt+conv(yc,h5,'same').*ctt;
    z0=T_zaehl_T(conv(yi,h1,'same'),Ts); %plot(z0(:,2))
    z0(:,3)=z0(:,2)-mean(z0(:,2));
    [p,~,mu]=polyfit(z0(:,1),z0(:,3),2);
    z0(:,3)=polyval(p,z0(:,1),[],mu); %plot(z0(:,1),z0(:,3),z0(:,1),z0(:,3))
    erg(k,:)=[k-2,sum(z0(:,3).^2)];
end, p=polyfit(erg(:,1),erg(:,2),2); %min suchen
if p(1)>0, k=-p(2)/2/p(1); A=A+da*max(-19,min(19,k));
else
    if erg(1,2)>erg(3,2), A=A+4*da;
    else,  A=A-4*da;
    end
end
if abs(k)<1
    if da<1e-7 || max(erg(:,2))<1.4*min(erg(:,2)), break
    else, da=da/2;
    end
end
end, j=log10(erg(2,2)); %F.. ist abgeschlossen
z2=z2+A*sin(P+t*F); %nach 3 Iterationen wieder einsetzen
end

