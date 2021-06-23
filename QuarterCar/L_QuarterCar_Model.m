function [M,C,K]=L_QuarterCar_Model(n)
%This code gives the Mass, Damping and Stiffness Matrix for a Forced linear
%Quarter Car model with stochastic forcing applied to the lower car. It
%also gives a forcing PSD information. The beam is forced in terms of an
%displacement imposed on the last node. The equations of motion are
%M*u''+C*u'+K*u=f_stochastic.

m1=2;
m2=3;
c=160;
k1=500;
k2=640;
Ce=[c -c;-c c];
Ke=[k2 -k2;-k2 k2];
L=ones(n,1)*m1;  L(end)=m2;
M=diag(L); C=zeros(n,n); K=zeros(n,n);
for i=1:n-1
    C(i:i+1,i:i+1)=C(i:i+1,i:i+1)+Ce;
    K(i:i+1,i:i+1)=K(i:i+1,i:i+1)+Ke;
end
K(end,end)=K(end,end)+k1;

end