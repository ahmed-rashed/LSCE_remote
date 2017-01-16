function [x,time] = InvFft(H,omega)

% ------------------   This file is part of EasyMod   ----------------------------
%  User function
%
%  Calculation of impulse response from FRF.
%
%  Synthax:
%  [x,time] = InvFft(H,omega) 
%
%  Input data:
%  H: FRF vector that we want to analyse,
%  omega: circular frequency vector.
%
%  Output data:
%  x: impulse time history vector,
%  time: time vector.
%
% Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS, Delphine LUPANT


% Test on vector direction
if size(H,1) == 1
    H = H ;
else
    H = H.' ;
end
clear n HD N deltaw temps x 
n = length(H) ;

% Aliasing effect operation
HD(1) = H(1) ;
HD(2:n-1) = H(2:n-1)/2 ;
HD(n) = H(n) ;
for i = 1:(n-2)
   HD(n+i) = H(n-i)'/2 ;
end
N = length(HD) ;

% Time vector definition
deltaf = (omega(2)-omega(1))/(2*pi) ;
time = [0:(N-1)]/(N*deltaf) ;

% Inverse FFT
x = N*real(ifft(HD)) ;
