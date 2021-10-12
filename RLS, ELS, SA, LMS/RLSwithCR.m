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

N_samples=800;
t=0:Ts:Ts*N_samples;
tt=size(t,2);

noise=normrnd(0,1,1,N_samples);
%%u=[0 0 0 0 100*ones(1,N_samples-4)];
%%u=u+noise;
;u=noise;
num=c(2:end);
den=d;
dd=d(1,2:5);
init=[dd,num];
% Sampling
y(1)=0;
y(2)=0;
y(3)=0;
y(4)=0;

changetime=200;
endtime=500;

for n=5:changetime
    y(n)=[-(y(n-1:-1:n-4)),(u(n-1:-1:n-4))]*(init');
    d2(n,:)=dd;
    n2(n,:)=num;
end

for n=changetime:endtime
    
    d2(n,:)=(1-0.01*(n-200)/300)*dd;
    n2(n,:)=(1-0.01*(n-200)/300)*num;
    initedit(n,:)=[d2(n,:),n2(n,:)];
   
    y(n)=[-(y(n-1:-1:n-4)),(u(n-1:-1:n-4))]*(initedit(n,:)');
    
end

for n=endtime:N_samples
    y(n)=[-(y(n-1:-1:n-4)),(u(n-1:-1:n-4))]*(initedit(500,:)');
    d2(n,:)=d2(500,:);
    n2(n,:)=n2(500,:);
end


%%d3=[1 d2]
%%c3=[0 n2]
%%sys2=tf(c3,d3,Ts)
%%B = isstable(sys2)



Y=y(1:end)';

%% RLS
p(:,:,1)=(10^6)*eye(8);
p(:,:,2)=p(:,:,1);
p(:,:,3)=p(:,:,1);
p(:,:,4)=p(:,:,1);

theta_hat(:,1)=zeros(1,8);
theta_hat(:,2)=theta_hat(:,1);
theta_hat(:,3)=theta_hat(:,1);
theta_hat(:,4)=theta_hat(:,1);

% k(:,1)=zeros(8,1);
% k(:,2)=zeros(8,1);
% k(:,3)=zeros(8,1);
% k(:,4)=zeros(8,1);

theta0=init;

theta01=theta0(1)*ones(1,N_samples)';
theta02=theta0(2)*ones(1,N_samples)';
theta03=theta0(3)*ones(1,N_samples)';
theta04=theta0(4)*ones(1,N_samples)';
theta05=theta0(5)*ones(1,N_samples)';
theta06=theta0(6)*ones(1,N_samples)';
theta07=theta0(7)*ones(1,N_samples)';
theta08=theta0(8)*ones(1,N_samples)';

landa=0.5;

for n=5:N_samples
    
    phi(:,n) = [-y(n-1) -y(n-2) -y(n-3) -y(n-4) u(n-1) u(n-2) u(n-3) u(n-4)]';
    
    p(:,:,n) = p(:,:,n-1)-p(:,:,n-1)*phi(:,n)*phi(:,n)'*p(:,:,n-1)/(1+phi(:,n)'*p(:,:,n-1)*phi(:,n));
  
        p(:,:,50)=p(:,:,1);
        p(:,:,100)=p(:,:,1);
        p(:,:,150)=p(:,:,1);
        p(:,:,200)=p(:,:,1);
        p(:,:,250)=p(:,:,1);
        p(:,:,300)=p(:,:,1);
        p(:,:,350)=p(:,:,1);
        p(:,:,400)=p(:,:,1);
        p(:,:,450)=p(:,:,1);
        p(:,:,500)=p(:,:,1);
        p(:,:,550)=p(:,:,1);
        p(:,:,600)=p(:,:,1);
        p(:,:,650)=p(:,:,1);
        p(:,:,700)=p(:,:,1);
        p(:,:,750)=p(:,:,1);
        p(:,:,800)=p(:,:,1);
        p(:,:,850)=p(:,:,1);
        p(:,:,900)=p(:,:,1);
        p(:,:,950)=p(:,:,1);
        
        
    k(:,n) = p(:,:,n-1)*phi(:,n)/(1+phi(:,n)'*p(:,:,n-1)*phi(:,n));
    
    theta_hat(:,n) = theta_hat(:,n-1)+k(:,n)*(y(n)-phi(:,n)'*theta_hat(:,n-1));
    
    y_hat(n)=phi(:,n)'*theta_hat(:,n-1);

end

thetafinal=theta_hat(:,N_samples)'

errorRLS=norm(thetafinal-theta0)
    
RMSDRLS=sqrt(sum((thetafinal-theta0).^2)/8)

plot(theta_hat')

sysdd=tf([0,thetafinal(5:end)],[1,thetafinal(1:4)],Ts);

%%disp([den(2:end),num(2:end)]')

%%figure
%%step(sysdd,t,'r+')
%%hold on
%%step(sysd,t)

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

figure
plot(d2(5:end,1))