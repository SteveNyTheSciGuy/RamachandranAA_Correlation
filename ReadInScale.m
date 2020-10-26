clc
AAT={'Ala';'Arg';'Asn';'Asp';'Cys';'Glu';'Gln';'Gly';'His';'Ile';'Leu';'Lys';'Met';'Phe';'Pro';'Ser';'Thr';'Trp';'Tyr';'Val'};
AAT=lower(AAT);
plt5=[1];
for V=1:length(AAT)
    AAt=AAT(V); %AA type, enter 3 letter code
%     str1=sprintf('%sScaled.txt',AAt{1});
%     Data1=dlmread(str1);
%     nu=Data1(:,1);
%     data1(:,V)=Data1(:,2);
    str2=sprintf('%sRaw.txt',AAt{1});
    Data2=dlmread(str2);
    data2(:,V)=Data2(:,2);
    nu=Data2(:,1);

end
NormSpec=zeros(length(nu),length(AAT),1);
FracNormSpec=zeros(length(nu),length(AAT),1);
FRAC=dlmread('Frac_Alpha.txt');
for i=1:length(AAT)
    data=data2(:,i);
    Max=max(data); Min=min(data);
    normspecmax=data./Max; normspecmin=data./Min;
    NormSpecMax(:,i)=normspecmax; 
    NormSpecMin(:,i)=normspecmin;
    FracNormSpecMax(:,i)=normspecmax.*FRAC(i);
    FracNormSpecMin(:,i)=normspecmin.*FRAC(i);
end

for i=1:length(nu)
    SumNormMax(i)=sum(NormSpecMax(i,:));
    SumSpecMax(i)=sum(FracNormSpecMax(i,:));
    SumNormMin(i)=sum(NormSpecMin(i,:));
    SumSpecMin(i)=sum(FracNormSpecMin(i,:));
end
% SumNorm=SumNorm/20;
% SumNormRaw=SumNormRaw/20;

alpha=dlmread('AlphaHelix.txt');
beta=dlmread('BetaSheet.txt');
alphaN=dlmread('AlphaNorm.txt');
betaN=dlmread('BetaNorm.txt');
if plt5
    fig5=figure(5);clf
    plot(nu,SumNormMax,'b-')
    hold on
    plot(nu,SumNormMin,'r-')
    plot(nu,SumSpecMax,'b.')
    plot(nu,SumSpecMin,'r.')
end
for i=1%:length(AAT)
    fig6=figure(6);clf
    plot(nu,NormSpecMax(:,i),'b-','linewidth',2)
    plot(nu,NormSpecMin(:,i),'r-','linewidth',2)
    %plot(nu,FracNormSpecMax(:,i),'b--','linewidth',1)
    %plot(nu,FracNormSpecMin(:,i),'r--','linewidth',1)
    %plot(nu,alphaN,'k--',nu,betaN,'k--','linewidth',1)
    axis([1200,2000,-1,1]) %min(NormSpec(:,i))
    legend('NormSpecMax','NormSpecMin','FracNormSpecMin','FracNormSpecMin','Alpha','Beta')
    set(gca,'xdir','reverse')
    hold off
    sv1=sprintf('%s_Alpha.jpg',AAt{i});saveas(fig6,sv1)

end
