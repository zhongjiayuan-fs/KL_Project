clear;
clc;
close all;


e=exp(1);
for i=1:8
    X(i,1)=rand()*(5-1)+1;
end
delta_t=1;
L=2100;       
N=round(L/delta_t);  
ts=zeros(N,1);
sigma=0.01;

es=0.000000001;

s=zeros(1,25);

s(1)=-0.5;
s(2)=-0.475;
s(3)=-0.45;
s(4)=-0.4;
s(5)=-0.375;
s(6)=-0.35;
s(7)=-0.325;
s(8)=-0.3;
s(9)=-0.275;
s(10)=-0.25;
s(11)=-0.225;
s(12)=-0.2;
s(13)=-0.175;
s(14)=-0.15;
s(15)=-0.125;
s(16)=-0.1;
s(17)=-0.08;
s(18)=-0.05;
s(19)=-0.02;
%p(20)=-0.0001;
s(20)=-0.001;
s(21)=0.02;
s(22)=0.05;
s(23)=0.08;
s(24)=0.1;
s(25)=0.15;
D= [-1 1 0 0 0 0 0 0 ;...
    -1 -1 0 0 0 0 0 0 ;...
    1 0 -1 0 0 0 0 0 ;...
    1 0 0 -1 0 0 0 0 ;...
    1 0 0 0 -1 0 0 0 ;...
    0 0 0 0 1 -1 0 1 ;...
    0 0 0 0 0 0 -1 1 ;...
    0 0 0 0 0 0 -1 -1];
TT=100;
kld=zeros(1,24);

sample_num=1;
CC=zeros(8,4,sample_num);

CC1=zeros(8,4,20);
ss=-0.5:5/1900:-0.45;


for ll=1:20
    qq(ll)=0.96^(1/abs(ss(ll)));
    E=[-3/5*qq(ll) 0 0 0 0 0 0 0 ;...
            0 -4/5 0 0 0 0 0 0 ;...
            0 0 -5/5 0 0 0 0 0 ;...
            0 0 0 -6/5 0 0 0 0 ;...
            0 0 0 0 -7/5 0 0 0 ;...
            0 0 0 0 0 -8/5 0 0 ;...
            0 0 0 0 0 0 -9/5 0 ;...
            0 0 0 0 0 0 0 -10/5];
    J=D*E*inv(D);
    for i=1:N-1;
        ts(i+1)=ts(i)+delta_t;
        eJ=e^(J*delta_t);
        for jj=1:8;
            X(jj,i+1)=eJ(jj,:)*X(:,i)+sigma*normrnd(0,1)*delta_t;
        end
    end
    CC1(:,1,ll)=X(:,2000); 
end


pre_TC=reshape(CC1(:,1,:),8,20);
for f=1:TT
    for l=2:25
        q(l)=0.96^(1/abs(s(l)));
        E=[-3/5*q(l) 0 0 0 0 0 0 0 ;...
            0 -4/5 0 0 0 0 0 0 ;...
            0 0 -5/5 0 0 0 0 0 ;...
            0 0 0 -6/5 0 0 0 0 ;...
            0 0 0 0 -7/5 0 0 0 ;...
            0 0 0 0 0 -8/5 0 0 ;...
            0 0 0 0 0 0 -9/5 0 ;...
            0 0 0 0 0 0 0 -10/5];
        J=D*E*inv(D);
        
        for k=1:sample_num
            for i=1:N-1
                ts(i+1)=ts(i)+delta_t;
                eJ=e^(J*delta_t);
                for jj=1:8
                    X(jj,i+1)=eJ(jj,:)*X(:,i)+sigma*normrnd(0,1)*delta_t;
                end
            end
            CC(:,l,k)=X(:,2000);
        end
        TC=reshape(CC(:,l,:),8,sample_num);
        
        for i=1:8
             mu0=mean(pre_TC(i,1:20));
             sigma0=std(pre_TC(i,1:20));
             tumor(i)= normpdf(TC(i,1),mu0,sigma0);
             pvalue(i,1:20)=normpdf(pre_TC(i,1:20),mu0,sigma0);
             casecdf(i)=normcdf(TC(i,1),mu0,sigma0)+es;
             concdf(i,1:20)=normcdf(pre_TC(i,1:20),mu0,sigma0)+es;
        end
        [tmp_com_idx,index]=sort(tumor);
        for k=1:8
            Q(k)=casecdf(index(k));
            P(k)=mean(concdf(index(k),1:20));
        end
        Q=Q/sum(Q);
        P=P/sum(P);
        kld(l-1)= kld(l-1)+0.5*(sum(P .* log(P./Q))+sum(Q.* log(Q./P)));
    end
    kld
        
end

kld=kld/TT;

%%
%subplot(1,2,1);
save kl.mat;
%plot(p(2:25),kld,'k-','LineWidth',2);
%title('KL');
%subplot(1,2,2);

plot(s(2:25),kld,'r-*','LineWidth',2.5);
%xlabel('Parameter');
%ylabel('sKLD');

































