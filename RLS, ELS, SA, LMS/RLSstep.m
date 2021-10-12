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
u=[0 0 0 0 20*ones(1,tt-5)];
%%u=normrnd(0,100,1,N_samples+5);
%%u=noise-mean(noise);
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
    y(n)=[-(y(n-1:-1:n-4)),(u(n-1:-1:n-4))]*(init');         %  (u(n:-1:n-4)) -> n-1:-1:n-4
end
 
Y=y(1:end)';

%% RLS
p(:,:,1)=100000000*eye(8);
p(:,:,2)=p(:,:,1);
p(:,:,3)=p(:,:,1);
p(:,:,4)=p(:,:,1);

theta_hat(:,1)=[0 0 0 0 0 0 0 0]';
theta_hat(:,2)=theta_hat(:,1);
theta_hat(:,3)=theta_hat(:,1);
theta_hat(:,4)=theta_hat(:,1);

% k(:,1)=zeros(8,1);
% k(:,2)=zeros(8,1);
% k(:,3)=zeros(8,1);
% k(:,4)=zeros(8,1);

theta0=init



for n=5:N_samples
    
     phi(:,n) = [-y(n-1) -y(n-2) -y(n-3) -y(n-4) u(n-1) u(n-2) u(n-3) u(n-4)]';
    
    p(:,:,n) = p(:,:,n-1)-p(:,:,n-1)*phi(:,n)*phi(:,n)'*p(:,:,n-1)/(1+phi(:,n)'*p(:,:,n-1)*phi(:,n));
    
    k(:,n) = p(:,:,n-1)*phi(:,n)/(1+phi(:,n)'*p(:,:,n-1)*phi(:,n));
    
    theta_hat(:,n) = theta_hat(:,n-1)+k(:,n)*(y(n)-phi(:,n)'*theta_hat(:,n-1));
    
     y_hat(n)=phi(:,n)'*theta_hat(:,n-1);

end

thetafinal=theta_hat(:,N_samples)'

errorRLS=norm(thetafinal-theta0)
    
errorRMS=sqrt(sum((thetafinal-theta0).^2)/8)

errorRlSfunc=norm(y_hat-y)

plot(theta_hat')


sysdd=tf([0,thetafinal(5:end)],[1,thetafinal(1:4)],Ts);

disp([den(2:end),num(2:end)]')
figure
step(sysdd,t)
hold on
step(sysd,t,'r+')
figure
plot(Y)
hold on
plot(phi'*theta_hat,'r*')
xlabel('Sample Number')
ylabel('Output')












