clear
%This script simulates a three-link headless human figure on a slackline 

%in the frontal plane, 
%as used in assignment 5 of the course "human and robot locomotion", 2020.
%
%Dimensions and weight of segments are inspired by anthropometric data 
%(mass of arm equivalent to two arms, mass of trunk contains swing leg).
%
%The slackline is a parallel assembly of linear spring+damper.
%The controller is to be implemented by students (A benchmark is provided).
%
%Definition of variables:
%g: gravitational constant
%
%mi: mass of a segment i
%Ji: moment of inertia of a segment i about its own CoM
%li: length of segment i
%si: distance of center of mass (CoM) of a segment i from the first joint below it
%r: instantaneous 'radius' of slackline (so distance of foot from anchor)
%sag: static 'radius' of slackline
%
%Indices used (for i):
%b: leg (note that 'b' is the abbreviation for the Dutch 'been')
%t: trunk
%a: arm
%f: foot
%h: hip
%s: shoulder
%
%Inertial coordinate system: x points to the right, y points up, origin is in
%slackline attachment point. Kinematics are described by the x- and y-position of
%the leg's CoM and the angle phi_i of each segment i (See drawing for definition).
%These angles are absolute, with respect to the vertical, and positive
%when the respective segment com is to the right of the segment's base joint.
%
%states used for the simulation:
%xb,yb,phib,phit,phia,
%dxb,dyb,dphib,dphit,dphia
%
%See the defintions also in the drawing.
%Heike Vallery, TU Delft, last update May 2020
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parameters:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%visualization choices:
showanimation=1;%do you want to see an animation? Then set this to 1. 

%simulation parameters:
Ts=.01;%sampling rate (s)
tend=4;%end time of simulation (s)
T = [0:Ts:tend]'; 
%general mechanical constants:
g=9.81;%gravity (m/s^2)

%segment properties:
mb=17;mt=51;ma=4;%mass of segments (kg)
lb=.8;lt=.5;la=.5;%length of segments (m)
sb=lb/2;st=lt/2;sa=la/2;%distance of segment CoM from respective joint below it (m)
Jb=1/12*mb*(lb)^2;Jt=1/12*mt*(lt)^2;Ja=1/12*ma*(la)^2;%segment inertia (kg m^2), assuming cylinders
CoM =0.8;
%((-mb*(lb-sb))+ (mt*st) + (ma*(sa+lt)))/ (ma +mb+mt);
%slackline properties:
sag=.3;%slackline sag in static conditions (m)
c=(mt+mb+ma)*g/sag;%slackline stiffness, as calculated from sag (N/m)
d=2*.5*sqrt(c*(mt+mb+ma));%damping of the slackline (Ns/m)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%initial conditions for the simulation:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initial segment angles:
phib0=0;%(rad)
phit0=0;%(rad)
phia0=0;%(rad)
%initial slackline 'radius':
r0=sag;%m
%define initial leg position via slackline radius and angles:
delta0=pi/1000;%rad
xb0=r0*sin(delta0)+sb*sin(phib0);%initial CoM x-coord of leg
yb0=-r0*cos(delta0)+sb*cos(phib0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%prepare simulation:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%initial time and state:
t0 = 0;
x0 =[xb0; yb0; phib0; phit0;phia0;zeros(5,1)];

%prepare parameter structure to pass to equations of motion (eom):
par.mb = mb; par.Jb=Jb;
par.mt=mt  ; par.Jt = Jt;
par.ma = ma; par.Ja = Ja;

par.sb = sb; par.st = st; par.sa = sa;
par.lb = lb; par.lt = lt; par.la = la;
par.g = g;
par.c=c;par.d=d;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%simulate:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options = odeset('AbsTol',1e-10,'RelTol',1e-8);
%students, uncomment this to test your own controller:
[t,state] = ode45(@equations_of_motion,[0:Ts:tend],x0,options,par);
%for the benchmarking controller (the code of which is hidden):
%[t,state] = ode45(@equations_of_motion_benchmark,[0:Ts:tend],x0,options,par);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%extract and process data from the simulation output:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%extract states from output:
xb=state(:,1);yb=state(:,2);
phib = state(:,3); phit = state(:,4); phia = state(:,5);
dxb=state(:,6);dyb=state(:,7);
dphib = state(:,8); dphit = state(:,9); dphia = state(:,10);

%calculate trunk position and velocity:
xt=xb+(lb-sb)*sin(phib)+st*sin(phit);
yt=yb+(lb-sb)*cos(phib)+st*cos(phit);
dxt=dxb+(lb-sb)*dphib.*cos(phib)+st*dphit.*cos(phit);
dyt=dyb-(lb-sb)*dphib.*sin(phib)-st*dphit.*sin(phit);

%calculate arm position and velocity:
xa=xt+(lt-st)*sin(phit)+sa*sin(phia);
ya=yt+(lt-st)*cos(phit)+sa*cos(phia);
dxa=dxt+(lt-st)*dphit.*cos(phit)+sa*dphia.*cos(phia);
dya=dyt-(lt-st)*dphit.*sin(phit)-sa*dphia.*sin(phia);

%calculate foot position and velocity:
xfvec=xb-sb*sin(phib);
yfvec=yb-sb*cos(phib);
rvec=sqrt(xfvec.^2+yfvec.^2);
deltavec=atan2(xfvec,-yfvec);

%ZMP Calculation
ddxb = diff(deltavec)./diff(T); %acceleration of foot
newRow = zeros(1,size(ddxb,2));                       % row of 0s
acc = [ddxb(1, :);newRow; ddxb(1:end, :)]; %resizing the matrix
acc(2,:) = []; %acceleration of the foot

ZMP = xfvec - ((CoM/g)*acc); % ZMP tracking of the given model

figure(4)
plot(t,ZMP);
legend('ZMP');
xlabel('time in s')
ylabel('ZMP position in m')


%store for visualization:
footx=rvec.*sin(deltavec);
footy=-rvec.*cos(deltavec);

hipx=footx+lb*sin(phib);
hipy=footy+lb*cos(phib);

shoulderx=hipx+lt*sin(phit);
shouldery=hipy+lt*cos(phit);

handx=shoulderx+la*sin(phia);
handy=shouldery+la*cos(phia);


%plot joint angles:
figure(1)
clf
plot(t,phib*180/pi)
hold all
plot(t,phit*180/pi,':bs')
plot(t,phia*180/pi,'+')
legend('phib','phit','phia');
ylabel('angles in deg')
xlabel('time in s')
title(gca,'joint angles')
set(gcf,'name','joint angles')


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%animate, if desired:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if showanimation
    %prepare visualization:
    figure(2)
    clf
    %plot dummies of the three segments and the slackline, to initialize the handles for later editing:
    stanceleghandle=plot([0,0],[-r0,-r0+lb],'r','linewidth',3);
    hold on
    torsohandle=plot([0,0],[0,0],'k','linewidth',3);
    armhandle=plot([0,0],[0,0],'g','linewidth',3);
    slackhandle=plot([0,0],[0,-r0]);
    
    xlim([-2,2])
    ylim([-2,2])
    
    xlabel('x in m')
    ylabel('y in m')
    
    set(gcf,'name','slackliner animation')

    %loop through the data to update the animated figure:
    for index=1:tend/Ts-1
        set(slackhandle,'XData',[0,footx(index)],'YData',[0,footy(index)])
        set(stanceleghandle,'XData',[footx(index),hipx(index)],'YData',[footy(index),hipy(index)])
        set(torsohandle,'XData',[hipx(index),shoulderx(index)],'YData',[hipy(index),shouldery(index)])
        set(armhandle,'XData',[shoulderx(index),handx(index)],'YData',[shouldery(index),handy(index)])
        
        pause(Ts)
        shg
    end
end




