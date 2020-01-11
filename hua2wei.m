function hua2wei

clear all
clc,close all
Nu = 4;
R = diag([25,25]);
P = diag([25,25]);
Rinv = eye(size(R))/R;



x = [];
y = [];
z = [];
for i = -20:1:20
    for j = -20:1:20
        x = [x i];
        y = [y j];
        m = [i,j]';
 
Psi = Rinv * (m*m' + P);
t1 = trace(Psi);
t2 = (Nu+2)/(Nu+t1);

vp = trace(Psi);
lq = -0.5*t2*vp+((Nu+2)/2-1)*log(t2)-t2;
lq2 = -0.5*t2*m'*Rinv*m;%-0.5*(x-xm)'*Pinv*(x-xm);
lq = exp(lq);
lq2 = exp(lq2);
qx = lq*lq2;
z = [z qx];
    end
end
figure(21)
plot3(x,y,z)