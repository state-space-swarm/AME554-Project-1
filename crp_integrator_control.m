function crp_integrator_control
global K
% Set up the time vector
ti0 = 0;
tif = 30;      
dt = .01;         
time = [ti0:dt:tif];

% Set up initial conditions
ICs = tan(3/2)*[1, -2, 3]/norm([1, -2, 3])
PHI0=acos((1-norm(ICs)^2)/(1+norm(ICs)^2))
x0 = ICs';
crp(1,:) = ICs;
dt2 = dt/2;
dt6 = dt/6; 

for i=2:length(time)
   
   t = time(i-1);
   tt = t + dt2;
   dxdt = deriv_crp(t,x0);
   x0t = x0 + dt2*dxdt;
   dx0t = deriv_crp(tt,x0t);
   x0t = x0 + dt2*dx0t;
   dxmold = deriv_crp(tt,x0t);
   x0t = x0 + dt*dxmold;
   dxm = dx0t + dxmold;
   dx0t = deriv_crp(t + dt, x0t);
   x0 = x0 + dt6*(dxdt + dx0t + 2*dxm);
   crp(i,:) = x0';
%   crp(i,:) = (x0+ref-cross(x0,ref))/(1-dot(x0,ref));
%   wvec(i,:)=-K*crp(i,:);% Linear CRP feedback
    wvec(i,:)=-K/(1+norm(x0)^2)*crp(i,:); % Nonlinear CRP feedback
end

figure()
subplot(1,3,1),plot(time, crp(:,1),time,crp(:,2),time,crp(:,3))
title('CRP error')
ylabel('q1-3')
xlabel('time(s)')
subplot(1,3,2),plot(time, acos((1.-(crp(:,1).^2+crp(:,2).^2+crp(:,3).^2))./(1.+crp(:,1).^2+crp(:,2).^2+crp(:,3).^2)))
axis([ti0, tif, 0, pi])
title('Principal angle of error')
xlabel('time(s)')
ylabel('PHI (rad)')
subplot(1,3,3),plot(time,wvec(:,1),time,wvec(:,2),time,wvec(:,3))
title('Angular velocity')
ylabel('w1-3 (rad/s)')
xlabel('time(s)')

end

function dcrp = deriv_crp(time, crp)
global K
q1 = crp(1);
q2 = crp(2);
q3 = crp(3);
qvec = [q1;q2;q3];

B(1,1) = 1+q1^2;
B(1,2) = q1*q2-q3;
B(1,3) = q1*q3+q2;
B(2,1) = q2*q1+q3;
B(2,2) = 1+q2^2;
B(2,3) = q2*q3-q1;
B(3,1) = q3*q1-q2;
B(3,2) = q3*q2+q1;
B(3,3) = 1+q3^2;

% Control law
K = 0; t1=5;
if (time>=t1)
    K = .5;
end
%omega = -K*qvec; % uncomment for linear CRP feedback
omega = -K/(1+norm(qvec)^2)*qvec; % Nonlinear CRP feedback

% Returns a vector of crp derivatives
dcrp = 0.5*B*omega;

end