%Authoe=Herbert Weidner
function [output] = window_sinc_filter(M, fc1, fc2, fs, filter_type, window_type) 
% Die Funktion berechnet die Filterkoeffizienten eines Window-Sinc Filters. 
% Hierzu werden die sinc Funktion sin(x)/x und ein wählbares Fenster verwendet. 
% 
% Inputparameter: 
% order        = Filterordnung ->  Anzahl Filterkoeff. = order + 1 
%                Bedingung: order = gerade ganze Zahl 
% 
% fc1 [Hz]     = cutoff frequency 1 = Grenzfrequenz der ersten Filterstufe 
%                -> für Hochpass, Tiefpass, Bandpass und Bandsperre 
%                Bedingung: 0 < fc1 < fs/2 
% 
% fc2 [Hz]     = cutoff frequency 2 = Grenzfrequenz der zweiten Filterstufe 
%                -> neben fc1 zusätzlich nur für Bandpass und Bandsperre 
%                Bedingung: 0 < fc1 < fc2 < fs/2 
% 
% fs [Hz]      = Sampling-frequency = Abtastfrequenz 
% 
% filter_type  = 'low' = Tiefpass, 'high' = Hochpass,, 'bandpass' = Bandpass 
%                und 'bandreject' = Bandsperre 
% 
% window_type  = Fenstertyp: 'Hamming' oder 'Blackman' 
% 
% analyse_plot = Bodediagramm, Frequenzantwort linear und Sprungfunktion(nur Tiefpass) 
%                des Filters darstellen: 'y' oder 'n' 
% 
% Outputparameter: 
% output     = Filterkoeffizienten = Impulsantwort des Filters 
% 
% Quelle des Filter-Algorithmus: 
% The Scientist and Engineer's Guide to Digital Signal Processing 
% By Steven W. Smith, Ph.D. , Chapter 16: Window-Sinc Filter 
% http://www.dspguide.com/ch16.htm 
% ------------------------------------------------------------------------- 

% Überprüfung der Inputparameter: 
% Grenzfrequenz 1 überprüfen 
Fc1 = fc1 / fs ; % Normierung auf fn (Fc1 = 0...0.5 fs) 
if (Fc1 <= 0) || (Fc1 >= 0.5) 
    error(['Wrong cutoff frequency fc1 [Hz]: ', sprintf('%8.1f',fc1),' ...choose: 0 < fc1 < ',sprintf('%8.1f',fs/2),' = fs/2']); 
end 
% Grenzfrequenz 2 überprüfen 
% wird nur für Bandpass- und Bandsperrfilter benötigt 
if (strcmp(filter_type, 'bandpass')) || (strcmp(filter_type, 'bandreject')) 
    Fc2 = fc2 / fs ; % Normierung auf fn (Fc2 = 0...0.5 fs) 
    if ((Fc2 <= 0) || (Fc2 >= 0.5) || (Fc2 <= Fc1)) 
        error(['Wrong cutoff frequency fc2 [Hz]: ', sprintf('%8.1f',fc2),' ...choose: ',sprintf('%8.1f',fc1) ' < fc2 < ',sprintf('%8.1f',fs/2),' = fs/2']); 
    end 
end 

% Auswahl Filtertyp 
switch lower(filter_type) 
    case 'high' % Hochpass - Filter 
        B = zeros(M+1, 1); % Init 
        % Fensterfunktion erstellen 
        window = Fenster(M+1, window_type); 
        for i = 0:M 
            if 2 * i == M 
                B(i+1) = 2*pi*Fc1; 
            else 
                B(i+1) = sin(2*pi*Fc1 * (i-(M/2))) / (i-(M/2)); 
            end 
            B(i+1) = B(i+1) * window(i+1); 
        end                 

        % Verstärkungsfaktor des Filters auf 1 normieren 
        B = B./sum(B); 
        % Tiefpass in Hochpass durch Inversion des Spektrums wandeln 
        output      = - B; 
        output((M/2)+1) = output((M/2)+1) + 1; 
        
    case 'low' % Tiefpass - Filter 
        B = zeros(M+1, 1); % Init 
        window = Fenster(M+1, window_type); 
        for i = 0:M 
            if 2 * i == M   % Multiplikation ist schneller als Division 
                B(i+1) = 2*pi*Fc1; 
            else 
                B(i+1) = sin(2*pi*Fc1 * (i-(M/2))) / (i-(M/2)); 
            end 
            B(i+1) = B(i+1) * window(i+1); 
        end 

        % Verstärkungsfaktor des Filters auf 1 normieren 
        B      = B./sum(B); 
        output = B; 

    case 'bandpass' % Bandpass - Filter 
        A = zeros(M+1, 1); B=A; % Init 
        % Fensterfunktion erstellen 
        window = Fenster(M+1, window_type); 
        
        % 1. Filterstufe: Berechnung der Filterkoeffizienten 
        for i = 0:M 
            if 2 * i == M   % Multiplikation ist schneller als Division 
                A(i+1) = 2*pi*Fc1; 
            else 
                A(i+1) = sin(2*pi*Fc1 * (i-(M/2))) / (i-(M/2)); 
            end 
            A(i+1) = A(i+1) * window(i+1); 
        end 

        % 2. Filterstufe: Berechnung der Filterkoeffizienten 
        for i = 0:M 
            if 2 * i == M 
                B(i+1) = 2*pi*Fc2; 
            else 
                B(i+1) = sin(2*pi*Fc2 * (i-(M/2))) / (i-(M/2)); 
            end 
            B(i+1) = B(i+1) * window(i+1); 
        end 

        % Verstärkungsfaktor des Filters auf 1 normieren 
        A = A./sum(A); B = B./sum(B); 
        % Tiefpass in Hochpass durch Inversion des Spektrums wandeln 
        B          = - B; 
        B((M/2)+1) = B((M/2)+1) + 1; 
        output     = A + B; 
        % Bandsperre in Bandpass durch Inversion des Spektrums wandeln 
        output          = - output; 
        output((M/2)+1) = output((M/2)+1) + 1; 

    case 'bandreject' % Bandsperr - Filter 
        A = zeros(M+1, 1); B=A; % Init 
        % Fensterfunktion erstellen 
        window = Fenster(M+1, window_type); 

        % 1. Filterstufe: Berechnung der Filterkoeffizienten 
        for i = 0:M 
            if 2 * i == M   % Multiplikation ist schneller als Division 
                A(i+1) = 2*pi*Fc1; 
            else 
                A(i+1) = sin(2*pi*Fc1 * (i-(M/2))) / (i-(M/2)); 
            end 
            A(i+1) = A(i+1) * window(i+1); 
        end 

        % 2. Filterstufe: Berechnung der Filterkoeffizienten 
        for i = 0:M 
            if 2 * i == M   
                B(i+1) = 2*pi*Fc2; 
            else 
                B(i+1) = sin(2*pi*Fc2 * (i-(M/2))) / (i-(M/2)); 
            end 
            B(i+1) = B(i+1) * window(i+1); 
        end 

        % Verstärkungsfaktor des Filters auf 1 normieren 
        A = A./sum(A); B = B./sum(B); 
        % Tiefpass in Hochpass durch Inversion des Spektrums wandeln 
        B          = - B; 
        B((M/2)+1) = B((M/2)+1) + 1; 
        output = A + B; 
        
    otherwise % ungültiger Filtertyp 
        error(['Unknown filter type: ', filter_type ' ...choose: high, low, bandpass or bandreject']); 
end % switch 

%-----------------------------------
function output = Fenster(window_size, window_type) 
% Input der Funktion: 
% window_size = Fensterlänge 
% window_typ = 'Hamming', 'Hann', 'Blackman', 'Blackman-Harris' 
N = window_size; 
% Fensterlänge output = zeros(N, 1); 
if mod(N, 2) == 0 % N gerade 
  m = fix(N / 2); n = m; 
else % N ungerade 
  m = fix(N / 2)+1; n = m-1; 
end
switch lower(window_type) 
  case 'hamming' 
    window = 0.54 - 0.46 * cos(2*pi*(0:m) / (N-1)); 
  case 'hann' 
    window = 0.50 - 0.50 * cos(2*pi*(0:m) / (N-1)); 
  case 'blackman' 
    window = 0.42 - 0.50 * cos(2*pi*(0:m) / (N-1)) + 0.08 * cos(4*pi* (0:m) / (N-1)); 
  case 'blackmanharris' 
  window = 0.35875 - 0.48829 * cos(2*pi*(0:m) / (N-1)) + 0.14128 * cos(4*pi* (0:m) / (N-1)) - 0.01168 * cos(6*pi* (0:m) / (N-1)); 
  otherwise, error(['Unknown window type: ', window_type]); 
end % Ergebnisvektor 
output = transpose([window(1:m),window(n:-1:1)]);
