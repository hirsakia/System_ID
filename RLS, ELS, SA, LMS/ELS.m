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

N_samples=5000;
t=0:Ts:(Ts*N_samples);
t2=t(1:end-1);
tt=size(t,2);
u=[0 0 0 0 0 0 0 10*ones(1,tt-8)];
distrb=normrnd(0,1,1,N_samples);

%%u=noise-mean(noise);
%%u=[0 0 0 0 normrnd(0,1,1,tt-5)];
u=u+distrb;
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
    y(n)=[-(y(n-1:-1:n-4)),(u(n-1:-1:n-4))]*(init')+0.2*distrb(n-1)-0.5*distrb(n-2);         %  (u(n:-1:n-4)) -> n-1:-1:n-4
end
 
Y=y(1:end)';

%% RLS
p(:,:,1)=10^8*eye(10);
p(:,:,2)=p(:,:,1);
p(:,:,3)=p(:,:,1);
p(:,:,4)=p(:,:,1);
%%p(:,:,5)=p(:,:,1);
%%p(:,:,6)=p(:,:,1);


theta_hat(:,1)=zeros(1,10);
theta_hat(:,2)=theta_hat(:,1);
theta_hat(:,3)=theta_hat(:,1);
theta_hat(:,4)=theta_hat(:,1);
%%theta_hat(:,5)=theta_hat(:,1);
%%theta_hat(:,6)=theta_hat(:,1);

phi(:,1)=zeros(1,10);
phi(:,2)=phi(:,1);
phi(:,3)=phi(:,1);
phi(:,4)=phi(:,1);
%%phi(:,5)=[0 0 0 0 0 0 0 0 ,0 0]';
%%phi(:,6)=[0 0 0 0 0 0 0 0 ,0 0]';

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

for n=5:N_samples
    
    phi(:,n) = [-y(n-1) -y(n-2) -y(n-3) -y(n-4) u(n-1) u(n-2) u(n-3) u(n-4) (y(n-1)-phi(:,n-1)'*theta_hat(:,n-1)) (y(n-2)-phi(:,n-2)'*theta_hat(:,n-2))]';
    
    k(:,n) = p(:,:,n-1)*phi(:,n)/(1+phi(:,n)'*p(:,:,n-1)*phi(:,n));
    
    p(:,:,n) = p(:,:,n-1)-p(:,:,n-1)*phi(:,n)*phi(:,n)'*p(:,:,n-1)/(1+phi(:,n)'*p(:,:,n-1)*phi(:,n));
    
    theta_hat(:,n) = theta_hat(:,n-1)+k(:,n)*(y(n)-phi(:,n)'*theta_hat(:,n-1));
    
    y_hat(n)=phi(:,n)'*theta_hat(:,n);


end

thetafinal=theta_hat(:,N_samples)';
theta_hat(:,N_samples)
theta00=[theta0 0.2 -0.5];

denRMS=size(theta00);

errorELS=norm(thetafinal-theta00)
    
errorEMS=sqrt(sum((thetafinal-theta00).^2)/8)

errorElSfunc=norm(y_hat-y)

plot(theta_hat')

sysdd=tf([0,thetafinal(5:8)],[1,thetafinal(1:4)],Ts);


figure
step(sysdd,t,'r+')
hold on
step(sysd,t)

figure
plot(Y)
hold on
plot(y_hat','r*')
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

















