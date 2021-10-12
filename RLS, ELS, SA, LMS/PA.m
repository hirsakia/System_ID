clear all
close all;
clc;
%% Tf Generation
m1=1;
m2=1;
k1=(86+1)/15;
b2=k1;
k2=(1+86)/100;
b1=k2;
k3=k2/2;
A=[0 1 0 0;
    -(k1+k3)/m1 -b1/m1 k3/m1 0;
    0 0 0 1;
    k3/m2 0 -(k2+k3)/m2 -b2/m2];
B=[0 ;1/m1 ;0 ;0]; 
C=[1 0 1 0]; 
D=0;
[b,a]= ss2tf(A,B,C,D);
sys=tf(b,a);
fb = bandwidth(sys);
Ts=0.05*2*pi/fb; %seconds
sysd = c2d(sys, Ts,'zoh');
[c,d]=tfdata(sysd,'v');

%% input

N_samples=1000;
t=0:Ts:Ts*N_samples;
tt=size(t,2);
noise=2*rand(N_samples,1)-1;
u=[0 0 0 0 10*ones(1,N_samples-4)];
%%u=u+noise;
%%u=noise
num=c(2:end);                                                % c -> c(2:end)
den=d;
dd=d(1,2:5);
init=[dd,num];                                              %[dd,c] -> [dd,num]
% Sampling
y(1)=0;
y(2)=0;
y(3)=0;
y(4)=0;
%%thetainitial=
for n=5:N_samples
    y(n)=[-(y(n-1:-1:n-4)),(u(n-1:-1:n-4))]*(init')+0.1*noise(n);         %  (u(n:-1:n-4)) -> n-1:-1:n-4
end
 
Y=y(1:end)';

%% RLS
p(:,:,1)=(10^6)*eye(8);
p(:,:,2)=p(:,:,1);
p(:,:,3)=p(:,:,1);
p(:,:,4)=p(:,:,1);

theta_hat(:,1)=zeros(1,8)';
theta_hat(:,2)=theta_hat(:,1);
theta_hat(:,3)=theta_hat(:,1);
theta_hat(:,4)=theta_hat(:,1);

% k(:,1)=zeros(8,1);
% k(:,2)=zeros(8,1);
% k(:,3)=zeros(8,1);
% k(:,4)=zeros(8,1);

theta0=init

theta01=theta0(1)*ones(1,N_samples)';
theta02=theta0(2)*ones(1,N_samples)';
theta03=theta0(3)*ones(1,N_samples)';
theta04=theta0(4)*ones(1,N_samples)';
theta05=theta0(5)*ones(1,N_samples)';
theta06=theta0(6)*ones(1,N_samples)';
theta07=theta0(7)*ones(1,N_samples)';
theta08=theta0(8)*ones(1,N_samples)';

landa=1.5
alpha=0.1

for n=5:N_samples
    
    phi(:,n) = [-y(n-1) -y(n-2) -y(n-3) -y(n-4) u(n-1) u(n-2) u(n-3) u(n-4)]';
    
    theta_hat(:,n) = theta_hat(:,n-1)+landa*phi(:,n)*(y(n)-phi(:,n)'*theta_hat(:,n-1))/(alpha+phi(:,n)'*phi(:,n));
    y_hat(n)=phi(:,n)'*theta_hat(:,n-1);
end

thetafinal=theta_hat(:,N_samples)'

errorRLS=norm(thetafinal-theta0)
    
RMSDRLS=sqrt(sum((thetafinal-theta0).^2)/8)

errorPAfunc=norm(y_hat-y)

plot(theta_hat')


sysdd=tf([0,thetafinal(5:end)],[1,thetafinal(1:4)],Ts);

figure
step(sysdd,t,'r+')
hold on
step(sysd,t)


figure
plot(Y)
hold on
plot(phi'*theta_hat,'r*')
xlabel('Sample Number')
ylabel('Output')

figure
subplot(2,4,1)
plot(theta_hat(1,:),'r')
hold on
plot(theta01,'b')
title('parameter1')

subplot(2,4,2)
plot(theta_hat(2,:),'r')
hold on
plot(theta02,'b')
title('parameter2')

subplot(2,4,3)
plot(theta_hat(3,:),'r')
hold on
plot(theta03,'b')
title('parameter3')

subplot(2,4,4)
plot(theta_hat(4,:),'r')
hold on
plot(theta04,'b')
title('parameter4')

subplot(2,4,5)
plot(theta_hat(5,:),'r')
hold on
plot(theta05,'b')
title('parameter5')

subplot(2,4,6)
plot(theta_hat(6,:),'r')
hold on
plot(theta06,'b')
title('parameter6')
hold off

subplot(2,4,7)
plot(theta_hat(7,:),'r')
hold on
plot(theta07,'b')
title('parameter6')
hold off

subplot(2,4,8)
plot(theta_hat(8,:),'r')
hold on
plot(theta08,'b')
title('parameter6')
hold off

