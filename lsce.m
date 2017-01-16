function [lsce_res,infoFRF,infoMODE] = lsce(H,freq,infoFRF)

% ------------------   This file is part of EasyMod   ----------------------------
%  User function
%
%  Identification based on the least-square complex exponential (LSCE)
%  method (MDOF method) giving the natural frequency, the loss factor 
%  and the modal constant in all the frequency range.
%
%  Synthax:
%  [lsce_res,infoFRF,infoMODE] = lsce(H,freq,infoFRF) 
%
%  Input data:
%  H: FRF matrix containing all the FRFs (column number = FRF number),
%  freq: frequency vector,
%  infoFRF: structure containing information on FRFs 
%            infoFRF(jj).response = jjth FRF response node
%            infoFRF(jj).dir_response = jjth FRF response direction (1=X, 2=Y,
%             3=Z, 4=RotX, 5=RotY, 6=RotZ)
%            infoFRF(jj).excitation = jjth FRF excitation node
%            infoFRF(jj).dir_excitation = jjth FRF excitation direction
%            (1=X, 2=Y, 3=Z, 4=RotX, 5=RotY, 6=RotZ). 
%
%  Output data:
%  lsce_res: results data in tabular forms
%      lsce_res(:,1) = natural frequency
%      lsce_res(:,2) = damping ratio
%      lsce_res(:,3) = stabilitisation state variable (1:
%      stabilization, 2: non stabilization),
%  infoMODE: structure containing the different identified parameters
%                infoMODE.frequencyk = natural frequency
%                infoMODE.etak = loss factor
%                infoMODE.Bijk = modal constant,
%  infoFRF : the completed input data.
%
% Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS, Delphine LUPANT


%  Necessary functions:
%  -----------------------------------------------------------
%  gen_resp_impul.m
%  AnMatIR.m
%  Modord.m
%  MatSur2.m
%  rec.m
%  releve.m
%  stabdiag.m
%  mode_lsce.m
%  mode_stab.m


% Calculation of the impulse response matrix from frequency data
HH = H ;
[MatIRs] = gen_resp_impul(HH,freq,infoFRF) ;
   
% Analysis of the impulse response matrix
[H,Ni,No,Nt,deltaT] = AnMatIR(MatIRs) ;
FMAX = (Nt/2+1)/(Nt*deltaT) ; % Frequency range is between 0 and FMAX

% Calculation parameters capture
%		- maximum iteration
disp(' ') ;
disp('-----------------------------------------------------------------------------') ;
disp(' ') ;
MaxMod = input('Model size - maximum iteration to analyse: (defaults: 30)\n') ;
if isempty(MaxMod), MaxMod = 30 ; end ;
%		- the tolerance in frequency
disp(' ');
prec1=input('Tolerance (%) in frequency: (defaults: 1%)\n') ;
if isempty(prec1), prec1 = 1 ; end ;
prec1 = prec1/100 ;
%		- the tolerance in damping
disp(' ') ;
prec2 = input('Tolerance (%) in damping: (defaults: 1%)\n') ;
if isempty(prec2), prec2 = 1 ; end ; 
prec2 = prec2/100 ;

% Definition of matrices allowed from the modal parameters at each step
clear WD WN XI LL Z ;
WD = zeros(2*MaxMod,MaxMod-1) ;
WN = zeros(2*MaxMod,MaxMod-1) ;
XI = zeros(2*MaxMod,MaxMod-1) ;
LL = zeros(Ni*(MaxMod-1),MaxMod*2) ;
Z = zeros(2*MaxMod,MaxMod-1) ;

% Definition of matrices allowed from the comparison results at each step
clear FTEMP FNMOD XIMOD XITEMP TESTXI ;
FTEMP = zeros(2*MaxMod,MaxMod-1) ;
FNMOD = zeros(2*MaxMod,MaxMod-1) ;
XIMOD = zeros(2*MaxMod,MaxMod-1) ;
XITEMP = zeros(2*MaxMod,MaxMod-1) ;
TESTXI = zeros(2*MaxMod,MaxMod-1) ;

% Loop begin
check = 0;
while check==0
   clear N ; 
   for N = 1:MaxMod ;
       % Data
       p = modord(N,Ni) ; % p represents the order of the linear differential equations 
       % Overdetermined matrix definition
       [G] = MatSur2(H,No,Ni,Nt,p) ;
       % Overdetermined equation A*x=B formulation
       B = -G(:,1:Ni) ;
       A = G(:,Ni+1:(p+1)*Ni) ;
       % Resolving (using the pseudoinverse function)
       x = pinv(A)*B ; 
       % Erreur calculation 
       %		- Inverse of conditionnning number
       InvCond(N) = 1/cond(A) ;
       %		- Singular normalized  values
       [U,S,V] = svd(A,0) ; 
       SingVals = diag(S) ;
       err(N) = 1/(max(SingVals)/min(SingVals)) ;
       %		- of least squares
       epsilon(N) = norm(B-(A*x));
       % Eigenvalue problem resolving
       clear L z ;
       [L,z] = PbValPp(x,Ni,p) ;
       % Modal parameters extraction
       clear lambda wd delta wn xi i ;
       lambda = log(z)./deltaT ;
       wd = imag(lambda) ;
       delta = real(lambda) ;
       wn = sqrt(wd.^2+delta.^2) ;
       xi = -(delta./wn) ;
       fn = wn/(2*pi) ;
       % Stabilization validation (frequency and damping)
       [FTEMP,XITEMP,TESTXI,FNMOD,XIMOD] = rec(fn,xi,N,FMAX,FTEMP,XITEMP,TESTXI,FNMOD,XIMOD,prec1,prec2) ;
       % Data saving
       WD(1:length(wd),N) = wd/2*pi ;
       WN(1:length(wn),N) = wn ;
       XI(1:length(xi),N) = xi ;
       LL((N-1)*Ni+1:N*Ni,1:size(L,2)) = L ;
       Z(1:length(z),N) = z ;
   end
   check = -1 ;
end

% Stabilization chart
stabdiag(FTEMP,XITEMP,TESTXI,FMAX,MaxMod,HH,freq) ;

% Least squares error chart visualization
subplot(2,2,3) ;
semilogy(epsilon,'-*') ;
xlabel('Number of modes') ;
ylabel('Amplitude') ;
title('Least squares error chart') ;
grid on ;
zoom on ;

% Error chart visualization
subplot(2,2,4) ;
semilogy(InvCond,'-*') ;
xlabel('Number of modes') ;
ylabel('Amplitude') ;
title('Conditioning error chart') ;
grid on ;
zoom on ;

% Selection of the model size
disp(' ') ;
disp('-----------------------------------------------------------------------------') ;
disp(' ') ;
num=input('Select the model size in order to extract the modal parameters:\n') ;
if isempty(num), num = MaxMod ; end ;
warning off

% Results statement
[lsce_res,Y] = releve(FTEMP,XITEMP,TESTXI,num) ;
kk = find(lsce_res(:,1)>0) ; % Only non-null frequencies are considered
[M,N] = size(lsce_res) ;
lsce_res = lsce_res(kk:M,:) ;
disp('     fok        xik   total stabilization') ;
disp(lsce_res) ;

% Mode shape extraction
[psi_phys,fn_phys,xi_phys] = mode_lsce(H,deltaT,Z,num,Nt) ;
infoMODE = mode_stab(infoFRF,psi_phys,fn_phys,xi_phys,lsce_res,prec1,prec2) ;

