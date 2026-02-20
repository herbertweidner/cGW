%ZetaPhoenicis, Carson-BW
%Autor = Herbert Weidner, 20-Feb-2026
%k=1.6697739*24*3600; fGW0=2/k; %[1.38630434624401e-05] astro
%------ Drift prüfen ----------------
%[y, Ts]=hol_DWD(3); %1=2000..2009; 2=2010..2019
clear, load yDWD, Ts=3600;
L=length(y); figure (1), zeig_sp2(y,-Ts,'log');
fGW0=13.862345e-6; %[1.38623565771659e-05] optimiert

%tD=23.93447192; fD=1/(tD*3600); 
tJ=365.25636042; fJ=1/tJ/24/3600; 
[sp,f]=zeig_sp2(y,-Ts,-200000.3); 
k=1;while f(k)<fGW0,k=k+1;end,k=k-37:k+37; plot(1e6*f(k),sp(k))

fGW=1.4e-6; Offset=fGW0-fGW; [b,a] = cheby1(6,0.01,0.02); 
y= filtTrans_n_Hz(y,Ts,5e4,fGW0,fGW,fGW-30e-9); %bis 2*1100nHz wegen Mod-Index
y=decimate(decimate(y,6),7); Ts=6*7*Ts; zeig_sp2(y,-Ts,'log'); 
[sp,f]=zeig_sp2(y,Ts,-300000.3);k=1;while f(k)<fGW,k=k+1;end,q1=k-39:k+39;
plot(1e9*f(q1),sp(q1)), title('verschobene fGW'), xlabel('Frequenz (nHz)')

t=(1:Ts:21/fJ)'; tt=1:(length(y)); fZF=1/Ts/40; %~160 nHz;
    dftt=1e-24*t(tt); dftt2=1e-30*t(tt).^2; %dftt3=1e-40*t(tt).^3;
    t=2*pi*t(tt); %21 Jahre FM darstellen
    tt=(tt)'*2*pi*fZF*Ts; %fixe ZF wegen variablem Oszillator
    stt=sin(tt); ctt=cos(tt);
j=0.2e-9; BF=window_sinc_filter(5e5,fZF-j,fZF+j,1/Ts,'bandpass','Blackman');
%LP_5=window_sinc_filter(3e5,0.5e-9,1,1/Ts,'low','blackman');
LP=window_sinc_filter(7e4,0.3e-9,1,1/(12*Ts),'low','blackman');

ax=2.08; px=1.25; fx=28.76e-9; drift=-2.953e6; drift2=-350.5;
ax2=1.51; px2=1.25; fx2=14.37e-9; %
ax3=1.40; px3=0.40; fx3=118.21e-9; %
axJ=5.25; pxJ=3.76; fxJ=31.27e-9; %31.69 nHz start with pxJ=3.7!!!!
ax4=2.85; px4=6.10; fx4=1.96e-9; % 

if 1==100000000000000000000000 %Richtung der Quelle bestimmen
    z2=axJ*sin(pxJ+t*fxJ); %nur Jahresrhythmus
    fGWdrift=fGW+drift*dftt+drift2*dftt2; %gemeinsam in allen Funktionen
    yi=sin(t.*fGWdrift+z2); z0=T_zaehl_T(yi,Ts);
    plot(z0(:,1),z0(:,2)), [b,a] = cheby1(6,0.01,0.1); 
    z0(:,1)=z0(:,1)/3600/24/tJ;
    z0(:,3)=filtfilt(b,a,z0(:,2)); plot(z0(:,1),z0(:,2),z0(:,1),z0(:,3))
    xlabel('Time (Jahre)'), ylabel('Halbe Schwingungsdauer (s)')
    z0(:,4:5)=NaN;
    z0(1,4)=15.237;z0(6,4)=8.71;z0(14,4)=0.653;z0(21,4)=-8.183;
    z0(29,4)=-14.405;z0(37,4)=-16.866;z0(45,4)=-14.947;z0(53,4)=-9.078;
    z0(60,4)=-0.717;z0(68,4)=7.751;z0(75,4)=14.554;z0(82,4)=17.332;

    z0(88,4)=15.237;z0(95,4)=8.71;z0(102,4)=0.653;z0(109,4)=-8.183;
    z0(117,4)=-14.405;z0(126,4)=-16.866;z0(134,4)=-14.947;z0(141,4)=-9.078;
    z0(149,4)=-0.717;z0(156,4)=7.751;z0(163,4)=14.554;z0(170,4)=17.332;

    z0(177,4)=15.237;z0(183,4)=8.71;z0(190,4)=0.653;z0(198,4)=-8.183;
    z0(205,4)=-14.405;z0(214,4)=-16.866;z0(222,4)=-14.947;z0(230,4)=-9.078;
    z0(237,4)=-0.717;z0(245,4)=7.751;z0(252,4)=14.554;z0(258,4)=17.332;
    z0(265,4)=15.237;
    plot(z0(:,1),z0(:,3),z0(:,1),1760*z0(:,4)+3.6e5)
z00=zeros(500,2); j=1;
    for k=1:size(z0,1)
        if ~isnan(z0(k,4))
            z00(j,:)=[z0(k,1),z0(k,4)]; j=j+1;
        end
    end, z00=z00(1:j-1,:);
    plot(2000+z0(:,1),1e9*Offset+5e8./z0(:,3)-1e9*fGW0), hold on
    plot(2000+z00(:,1),-9.4*z00(:,2),'r*'), hold off
    xlabel('Time (Year)'), ylabel('Frequency deviation (nHz)') 
    return
end

R=zeros(6950,20); %extend to (695000,20)?
erg=zeros(3,2); fakF=-2e-12; fakD1=-6e5; fakD2=-300;
for k3=1:size(R,1)
z2=ax*sin(px+t*fx)+ax2*sin(px2+t*fx2)+ax3*sin(px3+t*fx3)+ax4*sin(px4+t*fx4);
z2=z2+axJ*sin(pxJ+t*fxJ);

fGWdrift=fGW+drift*dftt+drift2*dftt2; %gemeinsam in allen Funktionen

[axJ,fxJ,pxJ,z2]=ITER(axJ,fxJ,pxJ,t,fGWdrift,y,z2,Ts,fZF,stt,ctt,LP);
[ax,fx,px,z2]=ITER(ax,fx,px,t,fGWdrift,y,z2,Ts,fZF,stt,ctt,LP);
[ax2,fx2,px2,z2]=ITER(ax2,fx2,px2,t,fGWdrift,y,z2,Ts,fZF,stt,ctt,LP);
[ax3,fx3,px3,z2]=ITER(ax3,fx3,px3,t,fGWdrift,y,z2,Ts,fZF,stt,ctt,LP);
[ax4,fx4,px4,z2]=ITER(ax4,fx4,px4,t,fGWdrift,y,z2,Ts,fZF,stt,ctt,LP);

yi=t.*fGWdrift+z2; ys=y.*sin(yi); yc=y.*cos(yi); 
yi=filtfilt(b,a,ys).*stt+filtfilt(b,a,yc).*ctt;
z0=T_zero(conv(yi,BF,'same'),Ts);
z00=diff(z0(:,1)); %plot(z00), drawnow
k=mean(z00)-0.5/fZF; %disp(k) %Abweichung von 1/fZF
fGW=fGW+k*fakF; %j=0; %Frequenz korrigieren
p=polyfit(1:size(z00),z00,2); drift=drift+p(2)*fakD1; drift2=drift2+p(1)*fakD2;
R(k3,:)=[fGW+Offset,drift,ax,px,fx,ax2,px2,fx2,ax3,px3,fx3,axJ,pxJ,fxJ,ax4,px4,fx4,k,p(1),p(2)];
end, R=R(1:k3-1,:); %return
fGW0=fGW+Offset; disp(fGW0);

%PSD berechnen
h6=window_sinc_filter(2e3,fZF-20e-9,1,1/Ts,'low','blackman');
y1=conv(ys,h6,'same').*stt+conv(yc,h6,'same').*ctt;
yi=conv(ys,h6,'same').*stt-conv(yc,h6,'same').*ctt; %invers
[sp,~]=zeig_sp2(y1,Ts,3000.4); [spi,f]=zeig_sp2(yi,Ts,3000.4);
plot(1e9*f,sp,1e9*f,spi)
k=16*1024;[sp,f]= pwelch(cat(1,y1,zeros(3e5,1)),rectwin(k),k/8,16*k,1/Ts);
plot(1e9*f,sp/1e5), xlabel('Frequency (nHz)')
return

%----------------------------------------

function erg=T_zero(y,Ts) %ohne Zeitzuordnung
yc=y.*cat(1,0,y(1:end-1)); erg=zeros(300,1); %1=Zeitpunkt der Nulldurchgänge
ys=(yc<0); %Vorzeichenwechsel
nn=1; j=1; while ys(j)==0, j=j+1; end %1. VZ-Wechsel
b=j-1+y(j-1)/(y(j-1)-y(j)); %genauer 0-Durchgang 
while j<length(ys)-30
    j=j+1; while ys(j)==0, j=j+1; end %nächster Wechsel
    a=j-1+y(j-1)/(y(j-1)-y(j)); %genauer 0-Durchgang 
    erg(nn,1)=Ts*(a+b)/2; nn=nn+1; b=a;
end 
erg=erg(1:nn-1,:);
end

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
if yGW(1)>0, y(1)=500; %willkürlicher Startwert
else, y(1)=-500;
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
%n ist die Länge von BF
%wenn B<0 wird cheby-Filter verwendet
L=length(y); j=(1:L)'*2*pi*f_ein*Ts;
ys=y.*sin(j); yc=y.*cos(j);
if B>0
    BF=window_sinc_filter(n,B,1,1/Ts,'low','blackman');
    ys=conv(ys,BF,'same'); yc=conv(yc,BF,'same');
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
%Modulationen mit f(mod)>10 nHz sollte man die Funktion ITER suchen
%Langsamere Modulationen mit ITER_low

function [A,F,P,Ph]=ITER(A,F,P,t,fGWdrift,y,Ph,Ts,fZF,stt,ctt,LP1)
if A==0, return, end 
dp=0.1; df=1e-12; da=0.01; erg=zeros(3,2); 
if P~=3.7
Ph=Ph-A*sin(P+t*F); %nach 3 Iterationen wieder einsetzen
for j=1:17
for k=1:3, xx=P+dp*(k-2);
    z=A*sin(xx+t*F)+Ph; 
    yi=t.*fGWdrift+z; ys=y.*sin(yi); yc=y.*cos(yi); 
    %Carson Bandbreite von sin(yi) prüfen!
    yi=conv(ys,LP1,'same').*stt+conv(yc,LP1,'same').*ctt;
    %yi=filtfilt(b,a,ys).*stt+filtfilt(b,a,yc).*ctt;
    ys=sync(yi,fZF,Ts); %plot(ys)
    y1=abs(hilbert(ys)); y1=y1(50:end-50); yc=1:length(y1);
    p=polyfit(yc,y1,1); erg(k,:)=[k-2,p(1)];
end, p=polyfit(erg(:,1),erg(:,2),2); %max suchen
if p(1)<0, k=-p(2)/2/p(1); P=P+dp*max(-9,min(9,k)); break
else
    if erg(1,2)<erg(3,2), P=P+2*dp;
    else,  P=P-2*dp;
    end
end
end
end


for j=1:17
for k=1:3, xx=F+(k-2)*df; 
    z=A*sin(P+t*xx)+Ph; 
    yc=t.*fGWdrift+z; ys=y.*sin(yc); yc=y.*cos(yc); 
    yi=conv(ys,LP1,'same').*stt+conv(yc,LP1,'same').*ctt;
    %yi=filtfilt(b,a,ys).*stt+filtfilt(b,a,yc).*ctt;
    ys=sync(yi,fZF,Ts); %plot(ys)
    y1=abs(hilbert(ys)); y1=y1(50:end-50); yc=1:length(y1);
    p=polyfit(yc,y1,1); erg(k,:)=[k-2,p(1)];
end, p=polyfit(erg(:,1),erg(:,2),2); %max suchen
if p(1)<0, k=-p(2)/2/p(1); F=F+df*max(-19,min(19,k)); break
else
    if erg(1,2)<erg(3,2), F=F+4*df;
    else, F=F-4*df;
    end
end
end

for j=1:17
for k=1:3, xx=A+da*(k-2);
    z=xx*sin(P+t*F)+Ph; 
    yc=t.*fGWdrift+z; ys=y.*sin(yc); yc=y.*cos(yc); 
    yi=conv(ys,LP1,'same').*stt+conv(yc,LP1,'same').*ctt;
    %yi=filtfilt(b,a,ys).*stt+filtfilt(b,a,yc).*ctt;
    ys=sync(yi,fZF,Ts); %plot(ys)
    y1=abs(hilbert(ys)); y1=y1(50:end-50); yc=1:length(y1);
    p=polyfit(yc,y1,1); erg(k,:)=[k-2,p(1)];
end, p=polyfit(erg(:,1),erg(:,2),2); %max suchen
if p(1)<0, k=-p(2)/2/p(1); A=A+da*max(-19,min(19,k)); break
else
    if erg(1,2)<erg(3,2), A=A+4*da;
    else, A=A-4*da;
    end
end
end %F.. ist abgeschlossen
Ph=Ph+A*sin(P+t*F); %nach 3 Iterationen wieder einsetzen
end


function [A,F,P,Ph]=ITER_low(A,F,P,t,fGWdrift,y,Ph,Ts,fZF,stt,ctt,LP1)
if A==0, return, end 
dp=0.1; df=3e-12; da=0.01; erg=zeros(65,2);
Ph=Ph-A*sin(P+t*F); %nach 3 Iterationen wieder einsetzen
for j=1:3 %Gesamtblock wiederholen
for k=1:65
    z=A*sin(t*F+k/10)+Ph; %alle Phasen k
    yc=t.*fGWdrift+z; ys=y.*sin(yc); yc=y.*cos(yc); 
    yi=conv(ys,LP1,'same').*stt+conv(yc,LP1,'same').*ctt;
    ys=sync(yi,fZF,Ts); %figure(1),plot(ys)
    y1=abs(hilbert(ys)); y1=y1(50:end-50); yc=(1:length(y1))';
    p=polyfit(yc,y1,1); ys=y1-polyval(p,yc); %plot(ys)
    erg(k,1)=sum(abs(ys));
end
[P,k]=min(erg(:,1)); P=k/10;
for k=1:3, xx=P+dp*(k-2);
    z=A*sin(xx+t*F)+Ph; 
    yc=t.*fGWdrift+z; ys=y.*sin(yc); yc=y.*cos(yc); 
    yi=conv(ys,LP1,'same').*stt+conv(yc,LP1,'same').*ctt;
    ys=sync(yi,fZF,Ts); %figure(1),plot(ys)
    y1=abs(hilbert(ys)); y1=y1(50:end-50); yc=(1:length(y1))';
    p=polyfit(yc,y1,1); ys=y1-polyval(p,yc); %plot(ys)
    erg(k,:)=[k-2,sum(abs(ys))];
end, p=polyfit(erg(1:3,1),erg(1:3,2),2); %min suchen
if p(1)>0, k=-p(2)/2/p(1); P=P+dp*max(-9,min(9,k));
else
    %P=P+pi/2;

end

while 1
for k=1:3, xx=F+(k-2)*df; 
    z=A*sin(P+t*xx)+Ph; 
    yc=t.*fGWdrift+z; ys=y.*sin(yc); yc=y.*cos(yc); 
    %yi=filtfilt(b,a,ys).*stt+filtfilt(b,a,yc).*ctt;
    yi=conv(ys,LP1,'same').*stt+conv(yc,LP1,'same').*ctt;
    ys=sync(yi,fZF,Ts); %figure(1),plot(ys)
    y1=abs(hilbert(ys)); y1=y1(50:end-50); yc=(1:length(y1))';
    p=polyfit(yc,y1,1); ys=y1-polyval(p,yc); %plot(ys)
    erg(k,:)=[k-2,sum(abs(ys))];
end, p=polyfit(erg(:,1),erg(:,2),2); %min suchen
if p(1)>0, k=-p(2)/2/p(1); F=F+df*max(-9,min(9,k)); break
else
    if erg(1,2)>erg(3,2), F=F+2*df;
    else, F=F-2*df;
    end
end
end

for k=1:3, xx=A+da*(k-2);
    z=xx*sin(P+t*F)+Ph; 
    yc=t.*fGWdrift+z; ys=y.*sin(yc); yc=y.*cos(yc); 
    yi=conv(ys,LP1,'same').*stt+conv(yc,LP1,'same').*ctt;
    ys=sync(yi,fZF,Ts); %figure(1),plot(ys)
    y1=abs(hilbert(ys)); y1=y1(50:end-50); yc=(1:length(y1))';
    p=polyfit(yc,y1,1); ys=y1-polyval(p,yc); %plot(ys)
    erg(k,:)=[k-2,sum(abs(ys))];
end, p=polyfit(erg(:,1),erg(:,2),2); %min suchen
if p(1)>0, k=-p(2)/2/p(1); A=A+da*max(-9,min(9,k));
else
    if erg(1,2)>erg(3,2), A=A+4*da;
    else, A=A-4*da;
    end
end
end %A.. ist abgeschlossen
Ph=Ph+A*sin(P+t*F); %nach 3 Iterationen wieder einsetzen
end

%A,F,P so bestimmen, dass der Anstieg linear wird
%[A,F,P]=finde_low(t,fGW,drift,dftt,y,z2,Ts,fZF,stt,ctt,LP2);
function [A,F,P]=finde_low(t,fGW,drift,dftt,y,Ph,Ts,fZF,stt,ctt,LP1)
A=0.1; F=2.3e-9; P=3;
dp=0.1; df=5e-12; da=0.01; erg=zeros(3,2); 
fGWdrift=fGW+drift*dftt;
    yc=t.*fGWdrift+Ph; ys=y.*sin(yc); yc=y.*cos(yc); 
    yc=conv(ys,LP1,'same').*stt+conv(yc,LP1,'same').*ctt;
    ys=sync(yc,fZF,Ts); %plot(ys) %bisheriges Ergebnis
    y1=abs(hilbert(ys)); y1=y1(50:end-50); yc=(1:length(y1))';
    p=polyfit(yc,y1,1); ys=y1-polyval(p,yc); plot(ys)
    [p,~,mu]=polyfit(yc,ys,3); y2=polyval(p,yc,[],mu); plot(y2) %ys=y2?

for j=1:5
    if j==5
        k=k;
    end
while 1
for k=1:3, xx=P+dp*(k-2); z=A*sin(xx+t*F)+Ph; 
    yc=t.*fGWdrift+z; ys=y.*sin(yc); yc=y.*cos(yc); 
    yc=conv(ys,LP1,'same').*stt+conv(yc,LP1,'same').*ctt;
    ys=sync(yc,fZF,Ts); %plot(ys) %bisheriges Ergebnis
    y1=abs(hilbert(ys)); y1=y1(50:end-50); yc=(1:length(y1))';
    p=polyfit(yc,y1,1); ys=y1-polyval(p,yc); plot(ys)
    erg(k,:)=[k-2,sum(abs(ys))];
end, p=polyfit(erg(:,1),erg(:,2),2); %min suchen
if p(1)>0, k=-p(2)/2/p(1); P=P+dp*max(-9,min(9,k)); break
else
    if erg(1,2)>erg(3,2), P=P+2*dp;
    else,  P=P-2*dp;
    end
end
end

while 1
for k=1:3, xx=F+(k-2)*df; z=A*sin(P+t*xx)+Ph;  
    yc=t.*fGWdrift+z; ys=y.*sin(yc); yc=y.*cos(yc); 
    yc=conv(ys,LP1,'same').*stt+conv(yc,LP1,'same').*ctt;
    ys=sync(yc,fZF,Ts); %plot(ys) %bisheriges Ergebnis
    y1=abs(hilbert(ys)); y1=y1(50:end-50); yc=(1:length(y1))';
    p=polyfit(yc,y1,1); ys=y1-polyval(p,yc); plot(ys)
    erg(k,:)=[k-2,sum(abs(ys))];
end, p=polyfit(erg(:,1),erg(:,2),2); %min suchen
if p(1)>0, k=-p(2)/2/p(1); F=F+df*max(-19,min(19,k)); break
else
    if erg(1,2)>erg(3,2), F=F+4*df;
    else, F=F-4*df;
    end
end
end

while 1
for k=1:3, xx=A+da*(k-2); z=xx*sin(P+t*F)+Ph; 
    yc=t.*fGWdrift+z; ys=y.*sin(yc); yc=y.*cos(yc); 
    yc=conv(ys,LP1,'same').*stt+conv(yc,LP1,'same').*ctt;
    ys=sync(yc,fZF,Ts); %plot(ys) %bisheriges Ergebnis
    y1=abs(hilbert(ys)); y1=y1(50:end-50); yc=(1:length(y1))';
    p=polyfit(yc,y1,1); ys=y1-polyval(p,yc); plot(ys)
    erg(k,:)=[k-2,sum(abs(ys))];
end, p=polyfit(erg(:,1),erg(:,2),2); %min suchen
if p(1)>0, k=-p(2)/2/p(1); A=A+da*max(-9,min(9,k)); break
else
    if erg(1,2)>erg(3,2), A=A+4*da;
    else, A=A-4*da;
    end
end
end 
end
end


