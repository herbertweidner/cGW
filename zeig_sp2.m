%Spektrum der Funktion y zeigen; 
%Author=Herbert Weidner
function [erg, f]=zeig_sp2(y,Ts,w) %Ts in s
%w ist die Anzahl der angezeigten Punkte
%die angefügte Dezimale wird zur Erhöhung
%der Auflösung zu NFFT addiert.
%Die effektive Bildbreite bleibt konstant
    %figure(1)
%[b,a] = cheby1(4,0.1,0.95); y=filtfilt(b,a,y);
[z,p,k] = cheby1(6,0.1,0.95); [sos,g] = zp2sos(z,p,k);
y=g(1,1)*sosfilt(sos,y);
%if Ts<0, Ts=-Ts; y=y.*blackman(length(y));
if Ts<0, Ts=-Ts; y=y.*hamming(length(y));
end
if isnumeric(w)
  if w<0, w=-w; a=1; else, a=0; end %kein Graph
    verg=round(10*(w-floor(w)));
    NFFT = 2^(nextpow2(length(y))+verg);% z.B. +3
    Fs=1/Ts; %Ein Wert pro Minute?
    f = Fs/2*linspace(0,1,NFFT/2+1);
    Y = fft(y,NFFT);
    erg=abs(Y(1:NFFT/2+1));
    w=floor(w);
    w=min(2^verg*floor(w),length(f));
    erg=erg(1:w); f=f(1:w);
    if a==0
      plot(1e6*f(1:w),erg(1:w)) %Ausschnitt vergrößern
      %semilogy(1e6*f,erg) %Ausschnitt vergrößern
      title('Spektrum'), xlabel('Frequenz in µHz')
      ylabel('relative Amplitude')
    end
else
    NFFT = 2^(nextpow2(length(y))+2);%Verfeinerung
    Fs=1/Ts; %Ein Wert pro Minute?
    f = Fs/2*linspace(0,1,NFFT/2+1);
    Y = fft(y,NFFT);
    erg=abs(Y(1:NFFT/2+1));
    if strcmp(w,'lin'),plot(f,erg),end %Ausschnitt vergrößern
    if strcmp(w,'log'),semilogy(f,erg),end
    title('Spektrum'), xlabel('Frequenz in Hz')
    ylabel('relative Amplitude')
end
%Yb=abs(Y(1:NFFT/2+1));
%plot(Yb(1:5000,1)) 
%semilogy(Yb(1:5000,1))
