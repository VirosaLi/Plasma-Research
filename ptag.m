%ptag.m
%This script calculates the tagging transfer function for a straight
%magnetic field and a quadratic Hamiltonian.
%
%a1,w1,x1,e1 are the spatial and velocity widths, a reference point on the
%beam in cartesian coordinates, and a direction vector for laser one.  The
%second laser also has these variables as does the viewing volume - to start
%with I assume that there is no velocity resolution on the detection ie no
%constraint based on w3.
%
%The plasma is characterized by thermal velocities vtll vtpp and a parallel
%drift vd and its gradient dvd.  The plasma may be rotating (this should be
%made to agree with the Hamiltonian) wr.  The plasma radius is r.  All
%variables are in cgs units !WITH VELOCITIES /wc  ie in cm.
%
%**************************************************************************
function S = ptag(xx,nu,sepz,a,tp,wc)
%Frequency units are in GHz and wavelength in nm (then *100)
% Laser 1--Tag Beam
a1 = a;                        % Beam width (cm)
w1 = .01;                      % Laser bandwidth [GHz?]
x1 = [xx 0 0];                 % Tag beam spatial location
e1 = [0 -1 0];                 % Tag beam orientation
lam1 = 100*585.4;
nu1 = nu;
%**************************************************************************
%Laser 2--Search Beam
a2 = a;                        % Beam width (cm)
w2 = .01;                      % Laser bandwidth
x2 = [0 0 sepz];               % Search beam spatial location
e2 = [1 1 0]/sqrt(2);          % Search beam orientation
lam2 = 100*585.4;
nu2 = nu;
%**************************************************************************
% Viewing volume
a3 = a;                        % Width (cm)
x3 = [0 0 sepz];               % Location
e3 = [1 0 0];                  % Orientation; originally e3=[1 -1 0]/sqrt(2); collinear [1 0 0];
%**************************************************************************
%Plasma
am = 138;
r = 2.5;                       % Plasma radius (cm)
wr = 0;                        % Angular velocity
tpara = .1;                    % Temperature in parallel direction (eV)
tperp = tp;                    % Temperature in perpendicular direction (eV)
vtll = 1e6*sqrt(tpara/am);     % Thermal velocity in parallel direciton
vtpp = 1e6*sqrt(tperp/am);     % Thermal velocity in perpendicular direction
vd = 1.2e5;                    % Drift velocity
vdp = [1 0 0]*3e3;
%**************************************************************************
%   .
%   _
%  |x  =  G*  |x
%  |v         |v		Matrix mapping
%   -
%
%   z(t)=M(t-t')*z(t')             Lie Transform   M=exp(G*(t-t'))
%
%**************
axx = -2*wr/wc;
ayy = axx;
azz = 0;
we = 0;       %originally we=ee; plasma rotation (due to drift?);
axy = we/wc;
ayx = axy;
ayz = 0;
azy = ayz;
azx = 0;
axz = azx;

G = Gmatrix( axx,axy,axz,ayx,ayy,ayz,azx,azy,azz );
%************************
%calculate constraints on the initial conditions
%
u = zeros(size(1:16));
c = zeros(16,6);
ud = sqrt((1-e1.^2)/2)/a1;
u(1:3) = ud.*x1;
c(1,:) = ud(1)*[1 0 0 0 0 0];
c(2,:) = ud(2)*[0 1 0 0 0 0];
c(3,:) = ud(3)*[0 0 1 0 0 0];
c(4,:) = [1 0 0 0 0 0]/(r*sqrt(2));
c(5,:) = [0 1 0 0 0 0]/(r*sqrt(2));
c(6,:) = [0 wr 0 wc 0 0]/(vtpp*sqrt(2));
c(7,:) = [-wr 0 0 0 wc 0]/(vtpp*sqrt(2));
c(8,:) = [vdp 0 0 wc]/(vtll*sqrt(2));
u(8) = vd/wc;
ud = sqrt((1-e2.^2)/2)/a2;
u(9:11) = ud.*x2;
cc9 = ud(1)*[1 0 0 0 0 0];
cc10 = ud(2)*[0 1 0 0 0 0];
cc11 = ud(3)*[0 0 1 0 0 0];
ud = sqrt((1-e3.^2)/2)/a3;
u(12:14) = ud.*x3;
cc12 = ud(1)*[1 0 0 0 0 0];
cc13 = ud(2)*[0 1 0 0 0 0];
cc14 = ud(3)*[0 0 1 0 0 0];
ud = 1/(sqrt(2)*w2);
u(15) = nu2*ud;
cc15 = [0 0 0 e2]*ud*wc/lam2;
ud = 1/(sqrt(2)*w1);
u(16) = nu1*ud;
c(16,:) = [0 0 0 e1]*ud*wc/lam1;
ut = sum(u.^2);
%add a constraint on the viewing velocity here .....
%************************
%constraints on final conditions act on initial conditions through
%the mapping.
N = 30000;               %Discretize time into 30000 intervals
S = zeros(1,N);    
for i = 1:N
    t = i*1e-2;
    M = expm(G*t);
    c(9,:) = cc9*M;
    c(10,:) = cc10*M;
    c(11,:) = cc11*M;
    c(12,:) = cc12*M;
    c(13,:) = cc13*M;
    c(14,:) = cc14*M;
    c(15,:) = cc15*M;
    %*************************************
    %construct the constraint matrix
    zz = u*c;
    Q = zeros(6);
    for j = 1:16
        Q = Q+c(j,:)'*c(j,:);
    end
    S(i) = exp((zz/Q)*zz'-ut)/sqrt(det(Q));
end
