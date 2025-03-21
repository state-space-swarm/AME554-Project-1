function tracking_control
clear all
% Set up the time vector
tf=10;
time=[0:.01:10];
theta0=2*pi/3;thetaf=0;
a1=theta0;
a3=3*(thetaf-theta0)/tf^2;
a4=-2*(thetaf-theta0)/tf^3;

% initial conditions and tolerance
ICs = [theta0,0]; 
tolerance = 1e-12;
options = odeset('RelTol',tolerance,'AbsTol',tolerance);

% Integrate unperturbed optimal control reference trajectory
[timer, stater] = ode45(@reference, [0,tf], ICs, options);
% Integrate perturbed optimal control
[timep1, statep1] = ode45(@perturbed, [0,tf], [2*pi/3+0.3,0], options);
[timep2, statep2] = ode45(@perturbed, [0,tf], [2*pi/3-0.3,0], options);
% Obtain the nominal control
ur = 2*a3.*ones(1,length(time))+6*a4.*time;
% Obtain the tracking control
om=.6;kp=om^2;kd=2*om;
thetar1=a1.*ones(length(timep1),1)+a3.*timep1.^2+a4.*timep1.^3;
thetadotr1=2*a3.*timep1+3*a4.*timep1.^2;
up1 = 2*a3.*ones(length(timep1),1)+6*a4.*timep1-kp*(statep1(:,1)-thetar1)-kd*(statep1(:,2)-thetadotr1);
thetar2=a1.*ones(length(timep2),1)+a3.*timep2.^2+a4.*timep2.^3;
thetadotr2=2*a3.*timep2+3*a4.*timep2.^2;
up2 = 2*a3.*ones(length(timep2),1)+6*a4.*timep2-kp*(statep2(:,1)-thetar2)-kd*(statep2(:,2)-thetadotr2);

figure(1)
plot(timer,stater(:,1),'k',timep1,statep1(:,1),'r',timep2,statep2(:,1),'r')
ylabel('\theta (rad)');
legend('reference','perturbed IC')
figure(2)
plot(time,ur,'k',timep1,up1,'r',timep2,up2,'r')
xlabel('time (sec)');ylabel('Control u(t) (rad/s^2)');
legend('reference','perturbed IC')

end

function dstate = reference(t,state)

theta=state(1);
omega=state(2);
theta0=2*pi/3;thetaf=0;
tf=10;
a1=theta0;
a3=3*(thetaf-theta0)/tf^2;
a4=-2*(thetaf-theta0)/tf^3;
thetar=a1+a3*t^2+a4*t^3;
thetadotr=2*a3*t+3*a4*t^2;
ur = 2*a3+6*a4*t;

kp=.4; kd=1.25;
thetadot=omega;
omegadot=ur;
dstate = [thetadot;omegadot];
end

function dstate = perturbed(t,state)

theta=state(1);
omega=state(2);
theta0=2*pi/3;thetaf=0;
tf=10;
a1=theta0;
a3=3*(thetaf-theta0)/tf^2;
a4=-2*(thetaf-theta0)/tf^3;
thetar=a1+a3*t^2+a4*t^3;
thetadotr=2*a3*t+3*a4*t^2;
ur = 2*a3+6*a4*t;
om=.9;kp=om^2;kd=2*om;
upd = -kp*(theta-thetar)-kd*(omega-thetadotr);

kp=.4; kd=1.25;
thetadot=omega;
omegadot=ur+upd; % for open loop delete "+upd"
dstate = [thetadot;omegadot];
end
