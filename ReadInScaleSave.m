clc
AAT={'Ala';'Arg';'Asn';'Asp';'Cys';'Gln';'Glu';'Gly';'His';'Ile';'Leu';'Lys';'Met';'Phe';'Pro';'Ser';'Thr';'Trp';'Tyr';'Val'};
AAT=lower(AAT);
plt5=[];
for V=1:length(AAT)
    AAt=AAT(V); %AA type, enter 3 letter code
%     str1=sprintf('%sScaled.txt',AAt{1});
%     Data1=dlmread(str1);
%     nu=Data1(:,1);
%     data1(:,V)=Data1(:,2);
    str2=sprintf('%sRaw.txt',AAt{1});
    Data2=dlmread(str2);
    data2(:,V)=Data2(:,2);
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
for i=2
    disp(AAT{i}), disp(i)
    fig6=figure(6);clf
    plot(nu,NormSpecMax(:,i),'b-','linewidth',2)
    hold on
    %plot(nu,FracNormSpecMax(:,i),'b--','linewidth',1)
    plot(nu,alphaN,'r--','linewidth',2)
    %plot(nu,betaN,'k--','linewidth',1)
    axis([1200,2000,-1,1]) %min(NormSpec(:,i))
    legend(AAT{i},'Alpha')%'location','southeast')
    set(gca,'xdir','reverse')
    Str1=sprintf('Divided\nby Max\n\n%6.3f%%',FRAC(i)*100);
    text(1900,.5,Str1);
    hold off
    sv1=sprintf('%s_AlphaMax.jpg',AAT{i});saveas(fig6,sv1)

    fig7=figure(7);clf
    plot(nu,NormSpecMin(:,i),'b-','linewidth',2)
    hold on
    %plot(nu,FracNormSpecMin(:,i),'r--','linewidth',1)
    %plot(nu,alphaN,'k--','linewidth',1)
    plot(nu,betaN,'r--','linewidth',2)
    axis([1200,2000,-.75,1]) %min(NormSpec(:,i))
    legend(AAT{i},'Beta')
    set(gca,'xdir','reverse')
    Str2=sprintf('Divided\nby Min\n\n%6.3f%%',FRAC(i)*100);
    text(1900,.5,Str2);
    hold off
    sv1=sprintf('%s_AlphaMin.jpg',AAT{i});saveas(fig7,sv1)
end
close figure(6)
close figure(7)