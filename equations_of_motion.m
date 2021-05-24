function [statederivative]=equations_of_motion(t0,state,par)
%
%This function contains the equations of motion of a simplified human 
%on a slackline in 2D, which means considering only the frontal plane,
%as well as a balancing control strategy that generates joint moments in
%hip and shoulder.
%
%Coordinate system: x points to the right, y points up, origin is in
%slackline attachment point. Generalized coordinates q are:
%x- and y-positions of the leg CoM and segment angles phi_i. 
%
%Note that these angles are absolute, with respect to the vertical, 
%and positive when the respective segment com is to the right of 
%the respective base joint, whereby the "base" is the foot for the leg, 
%the hip for the trunk, and the shoulder for the arm.
%
%Internal forces and internal moment acting on a segment are defined 
%positive if they point in positive coordinate directions 
%at the *base* joint of a segment. 
%
%H. Vallery, last update May 2020

%extract parameters:
g = par.g; %gravity (m/s^2)
T = [0:0.01:4]';
CoM = 0.8;
mb = par.mb; %leg mass (kg)
Jb = par.Jb; %leg inertia (kgm^2)
mt = par.mt; %trunk mass (m)
Jt = par.Jt; %trunk inertia (kgm^2)
ma = par.ma; %arm mass (m)
Ja = par.Ja; %arm inertia (kgm^2)

sb = par.sb; %distance of leg center of mass from foot (m)
st = par.st; %distance of trunk center of mass from hip (m)
sa = par.sa; %distance of arm center of mass from shoulder (m)

lb = par.lb; %leg length (m)
lt = par.lt; %trunk length (m)
la = par.la; %arm length (m)

c=par.c; %stiffness of slackline (N/m)
d=par.d; %damping coefficient of slackline (Ns/m)

%state vector:
%xb,yb,phib,phit,phia,dxb,dyb,dphib,dphit,dphia, 
%with:
xb=state(1); %com position leg in x direction (m)
yb=state(2); %com position leg in y direction (m)
phib = state(3); %leg angle (rad)
phit = state(4); %trunk angle (rad)
phia = state(5); %arm angle (rad)
%derivatives:
dxb=state(6); dyb=state(7); dphib=state(8); dphit=state(9); dphia=state(10); 
   
%All variables and parameters must be expressed in SI units (see above). 

%calculate foot position and velocity (needed for slackline force):
    xf=xb-sb*sin(phib);%foot position in x direction
    yf=yb-sb*cos(phib);%foot position in y direction
    dxf=dxb-sb*dphib*cos(phib);%foot velocity in x-direction
    dyf=dyb+sb*dphib*sin(phib);%foot velocity in y-direction
    
%calculate slackline geometry and force:    
    r=sqrt(xf^2+yf^2);%m, slackline radius
    drdt=(xf*dxf+yf*dyf)/sqrt(xf^2+yf^2);%radius derivative wrt time
    delta=atan2(xf,-yf);%slackline angle
    Fs=r*c+drdt*d;%force in radial direction
    FFx=-Fs*sin(delta);%force x-component acting on foot
    FFy=Fs*cos(delta);%force y-component acting on foot
    
    FF = sqrt((FFx^2)+(FFy)^2); %resultant force
    Mf = FF * (lb*cos(phib)); %torque at the hip due to the feet
%----------------------------------------------------------------------
%INPUTS (hip and shoulder torques)
%----------------------------------------------------------------------

%STUDENTS, THIS IS WHERE YOUR CODE MUST GO (you can also create a separate
%function for this):
%hip:


if T < 1
tauS= -FF*(la/2 - CoM); %arm torque strategy
tauH= -Mf;
else
    tauS = 0;
    tauH= 0;
end


%note the sign definitions provided above.

%%
%----------------------------------------------------------------------
%EQUATIONS OF MOTION (DERIVED VIA LAGRANGE):
%----------------------------------------------------------------------
%mass matrix:
massmatrix =[ ma + mb + mt,          0,      lb*ma*cos(phib) + lb*mt*cos(phib) - ma*sb*cos(phib) - mt*sb*cos(phib),                                                                                                                                           lt*ma*cos(phit) + mt*st*cos(phit),                                                                      ma*sa*cos(phia);
     0,       ma + mb + mt,          ma*sb*sin(phib) - lb*mt*sin(phib) - lb*ma*sin(phib) + mt*sb*sin(phib),                                                                                                                                         - lt*ma*sin(phit) - mt*st*sin(phit),                                                                     -ma*sa*sin(phia);
    (ma*cos(phib)*(2*lb - 2*sb))/2 + (mt*cos(phib)*(2*lb - 2*sb))/2, - (ma*sin(phib)*(2*lb - 2*sb))/2 - (mt*sin(phib)*(2*lb - 2*sb))/2, Jb + (ma*(sin(phib)*(lb*sin(phib) - sb*sin(phib))*(2*lb - 2*sb) + cos(phib)*(lb*cos(phib) - sb*cos(phib))*(2*lb - 2*sb)))/2 + (mt*(sin(phib)*(lb*sin(phib) - sb*sin(phib))*(2*lb - 2*sb) + cos(phib)*(lb*cos(phib) - sb*cos(phib))*(2*lb - 2*sb)))/2, (ma*(lt*cos(phib)*cos(phit)*(2*lb - 2*sb) + lt*sin(phib)*sin(phit)*(2*lb - 2*sb)))/2 + (mt*(st*cos(phib)*cos(phit)*(2*lb - 2*sb) + st*sin(phib)*sin(phit)*(2*lb - 2*sb)))/2, (ma*(sa*cos(phia)*cos(phib)*(2*lb - 2*sb) + sa*sin(phia)*sin(phib)*(2*lb - 2*sb)))/2;
    lt*ma*cos(phit) + mt*st*cos(phit),      - lt*ma*sin(phit) - mt*st*sin(phit),        lb*lt*ma*cos(phib - phit) - lt*ma*sb*cos(phib - phit) + lb*mt*st*cos(phib - phit) - mt*sb*st*cos(phib - phit),     ma*lt^2 + mt*st^2 + Jt,        lt*ma*sa*cos(phia - phit);
   ma*sa*cos(phia),      -ma*sa*sin(phia),    lb*ma*sa*cos(phia - phib) - ma*sa*sb*cos(phia - phib),   lt*ma*sa*cos(phia - phit),   ma*sa^2 + Ja];

%further terms (gravity etc.):
restterms =  [ - ma*sa*sin(phia)*dphia^2 - dphit*(dphit*lt*ma*sin(phit) + dphit*mt*st*sin(phit)) - dphib*(dphib*lb*ma*sin(phib) + dphib*lb*mt*sin(phib) - dphib*ma*sb*sin(phib) - dphib*mt*sb*sin(phib));
    - ma*sa*cos(phia)*dphia^2 - dphit*(dphit*lt*ma*cos(phit) + dphit*mt*st*cos(phit)) + g*(ma + mb + mt) - dphib*(dphib*lb*ma*cos(phib) + dphib*lb*mt*cos(phib) - dphib*ma*sb*cos(phib) - dphib*mt*sb*cos(phib));
    (lb - sb)*(dphib*dyb*ma*cos(phib) - g*mt*sin(phib) - g*ma*sin(phib) + dphib*dyb*mt*cos(phib) + dphib*dxb*ma*sin(phib) + dphib*dxb*mt*sin(phib) + dphib*dphit*lt*ma*sin(phib - phit) - dphia*dphib*ma*sa*sin(phia - phib) + dphib*dphit*mt*st*sin(phib - phit)) - dphib*((ma*(sin(phib)*(2*lb - 2*sb)*(dxb + dphib*lb*cos(phib) + dphit*lt*cos(phit) + dphia*sa*cos(phia) - dphib*sb*cos(phib)) - sin(phib)*(dphib*lb*cos(phib) - dphib*sb*cos(phib))*(2*lb - 2*sb) - cos(phib)*(2*lb - 2*sb)*(dphib*lb*sin(phib) - dyb + dphit*lt*sin(phit) + dphia*sa*sin(phia) - dphib*sb*sin(phib)) + cos(phib)*(dphib*lb*sin(phib) - dphib*sb*sin(phib))*(2*lb - 2*sb)))/2 + (mt*(sin(phib)*(2*lb - 2*sb)*(dxb + dphib*lb*cos(phib) - dphib*sb*cos(phib) + dphit*st*cos(phit)) + cos(phib)*(2*lb - 2*sb)*(dyb - dphib*lb*sin(phib) + dphib*sb*sin(phib) - dphit*st*sin(phit)) - sin(phib)*(dphib*lb*cos(phib) - dphib*sb*cos(phib))*(2*lb - 2*sb) + cos(phib)*(dphib*lb*sin(phib) - dphib*sb*sin(phib))*(2*lb - 2*sb)))/2) - dphit*((mt*(dphit*st*cos(phib)*sin(phit)*(2*lb - 2*sb) - dphit*st*cos(phit)*sin(phib)*(2*lb - 2*sb)))/2 + (ma*(dphit*lt*cos(phib)*sin(phit)*(2*lb - 2*sb) - dphit*lt*cos(phit)*sin(phib)*(2*lb - 2*sb)))/2) + (dphia*ma*(dphia*sa*cos(phia)*sin(phib)*(2*lb - 2*sb) - dphia*sa*cos(phib)*sin(phia)*(2*lb - 2*sb)))/2;
    dphit*dyb*lt*ma*cos(phit) - dphit*(dyb*lt*ma*cos(phit) + dyb*mt*st*cos(phit) + dxb*lt*ma*sin(phit) + dxb*mt*st*sin(phit) - dphib*lb*lt*ma*sin(phib - phit) - dphia*lt*ma*sa*sin(phia - phit) + dphib*lt*ma*sb*sin(phib - phit) - dphib*lb*mt*st*sin(phib - phit) + dphib*mt*sb*st*sin(phib - phit)) - g*lt*ma*sin(phit) - g*mt*st*sin(phit) - dphia^2*lt*ma*sa*sin(phia - phit) - dphib*(dphib*lb*lt*ma*sin(phib - phit) - dphib*lt*ma*sb*sin(phib - phit) + dphib*lb*mt*st*sin(phib - phit) - dphib*mt*sb*st*sin(phib - phit)) + dphit*dyb*mt*st*cos(phit) + dphit*dxb*lt*ma*sin(phit) + dphit*dxb*mt*st*sin(phit) - dphib*dphit*lb*lt*ma*sin(phib - phit) - dphia*dphit*lt*ma*sa*sin(phia - phit) + dphib*dphit*lt*ma*sb*sin(phib - phit) - dphib*dphit*lb*mt*st*sin(phib - phit) + dphib*dphit*mt*sb*st*sin(phib - phit);
    dphib*(dphib*lb*ma*sa*sin(phia - phib) - dphib*ma*sa*sb*sin(phia - phib)) - dphia*(dyb*ma*sa*cos(phia) + dxb*ma*sa*sin(phia) + dphib*lb*ma*sa*sin(phia - phib) + dphit*lt*ma*sa*sin(phia - phit) - dphib*ma*sa*sb*sin(phia - phib)) + ma*sa*(dphia*dxb*sin(phia) - g*sin(phia) + dphia*dyb*cos(phia) + dphia*dphib*lb*sin(phia - phib) + dphia*dphit*lt*sin(phia - phit) - dphia*dphib*sb*sin(phia - phib)) + dphit^2*lt*ma*sa*sin(phia - phit)];

%generalized applied forces, from slackline and joint moments:
appliedforces=[FFx;FFy;tauH-FFx*sb*cos(phib)+FFy*sb*sin(phib);-tauH+tauS;-tauS];

%solve equations of motion for accelerations:
accelerations=(massmatrix)\(appliedforces-restterms);

statederivative=[state(6:10);accelerations];




