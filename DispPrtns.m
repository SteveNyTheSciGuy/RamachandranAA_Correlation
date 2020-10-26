clear all; clc; close all
% cp=[19;31;21;14]; % alpha chosen protiens
% prtns={'4de6','2frf','1k01','10cc'};

cp=[9;6;41;28]; % beta chosen protiens
prtns={'2lve','1bfg','1ba7','4iba'};
PS=dlmread('protein.txt'); % Input protein spectra
nu=PS(:,1);
for i=1:length(cp)
    pn=cp(i); %protien number
    spec=PS(:,pn+1);
    spec1=spec./max(spec);
    Spec(:,i)=spec1;
end

alphaN=dlmread('AlphaNorm.txt');
betaN=dlmread('BetaNorm.txt');
aBetaN=.9.*betaN+.3.*alphaN;
figure(1)
plot(nu,Spec(:,1),'m',nu,Spec(:,2),'k',nu,Spec(:,3),'b',nu,Spec(:,4),'y')
hold on
%plot(nu,alphaN,'--r','linewidth',2)
plot(nu,aBetaN,'--r','linewidth',2)
axis([1200,2000,0,1.25]) %min(NormSpec(:,i))
% legend(prtns,'alpha')
legend(prtns,'beta')
set(gca,'xdir','reverse')
