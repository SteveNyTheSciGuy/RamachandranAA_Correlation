%Program written by Steven Nystrom & Dr. Coe 6/5/14
clear all; clc;
% % %
% How to use:
% Lines: Enter the ellipse parameters
%    The Phi and Psi centers, the major/minor axis length, the ellipse tilt
% Line: Select plots to be shown
% Line: The proteins that are to be included in the plot are listed by their 4
%    character PDB file identifier
% The Outputs:
%   The command window reveals the fraction inside and outside of the
%   ellipse, the Phi/Psi centers that were chosen, and the averaged Phi/Psi
%   value
%   The best way to ensure all positive going spectra is to make the
%   average value be approximately equal to the ellipse center.
%   The non-linear least squares approach can lead to the IR spectra having
%   positive and negative going peaks; this requires some finesse with the
%   ellipse parameters, by moving it up/down a few degrees and changing the
%   major and minor axis to get good spectra
% Figure1:
%   The isolated spectra appears as a red line with the standard deviation
%   at every 2 wave numbers appearing in blue
% Figure2:
%   The scaled spectra of inside/outside the ellipse and the total spectra
%   Change the parameters of the ellipse so that the entire spectra is
%   positive or negative
% Figure3:
%   The contoured Ramachandran plot with the ellipse overlaid.
%   The more AA in a 5x5 bin the greater the number
% % %

% Select you region for the ellipse
% beta sheet
PhiC=-113.16; PsiC=134.42; % Centers of the ellipse
ae=39; be=23;      % Major & Minor Axis respectfully
phip1=-35*2*pi/360;    % Tilt angle of major axis from Psi
%Alpha Helix
% PhiC=-62.79; PsiC=-40.94; % Centers of the ellipse
% ae=20; be=13;      % Major & Minor Axis respectfully
% phip1=-45*2*pi/360;    % Tilt angle of major axis from Psi
upw=3.15; % Upper omega bound
low=-3.15; % Lower omega bound
AAt={'tyr'}; %AA type, enter 3 letter code
regn={'alpha'};
% Select your features by putting a 1 in the box if you want that feature
plt1=[1]; plt2=[1]; plt3=[1]; plt4=[]; savv=[1];

% List of 4 character protein data bank file identifiers that label the
% angle files that were extracted from the PDB files using the Ramachandran function
Prtns=['1KCT';'4F5S';'1E7I';'6ADH';'2HCY';'1BFG';'1V9E';'3CNA';'1YPH';'2GIW';'1AKK';'3CYT';'2LIR';'1OCC';'1DNK';'2V35';'3ENL';'1F13';'4DE6';'3GHG';'1K0Y';'1NS6';'1GB1';'1EKU';'1F6S';'3H3F';'1CJ5';'4IBA';'2LVE';'1AZF';'2FRF';'9PAP';'2QCA';'1SBC';'1CB4';'1SXN';'1NUC';'1R2S';'4I8L';'1TGN';'1BA7'];
[NumPro sizb]=size(Prtns); % Gets the number of proteins in variable NumPro
AAc=zeros(NumPro,2,1);AAf=zeros(NumPro,2,1);% Initialize arrays
% AllAngle=[]; AllName=[];
[npro cmax]=size(AAc); % Gets number of groups

% reads in all the angles and the AA for each corrosponding angle
% It isolates the AA chosen in variable AAt
AllName1=fopen('All_AA.txt'); %reads in the txt file of the AA 3 letter code
AllName=textscan(AllName1,'%s');
AllName=AllName{1};
fclose(AllName1);
Angs=dlmread('angles_all.txt');
aaFlag1 = strcmpi(AAt{1},AllName); % if the chosen AA is in the list it becomes a 1, else 0
aaFlag1=ones(length(AllName),1,1);
cnt=0; Flg1=zeros(1,3,1); %intialize array
for i=1:length(AllName)
    if aaFlag1(i)
        cnt=cnt+1;
        Flg1(cnt,:)=Angs(i,:); % Grabs the position of the certian AA for plotting
    end
end
Flg1=ones(1,3,1); %ignore specific aa
% Count the number of amino acid residues inside the ellipse
cnt=0; cnt11=0; %intialize arrays
for i=1:NumPro
    prtns=sprintf('a_%s.txt',Prtns(i,:)); % Completes the name of the angle files for each protein
    Proteins{i,:}=dlmread(prtns); % Cell array of the amino acid angles indexed by ii
    [nAA(i) nc]=size(Proteins{i,:}); % Number of amino acids in each protein
    for j=1:nAA(i) % Loop over amino acid residues in each protein
        cnt11=cnt11+1;
        if Proteins{i,1}(j,3)<=upw && Proteins{i,1}(j,3)>=low % Sorts desired omega range
            if aaFlag1(cnt11) % Sorts for chosen AA, if the AA is present its a 1 in the array
                k1=((cos(phip1)*(Proteins{i,1}(j,1)-PhiC)+sin(phip1)*(Proteins{i,1}(j,2)-PsiC))/ae)^2;
                k2=((sin(phip1)*(Proteins{i,1}(j,1)-PhiC)-cos(phip1)*(Proteins{i,1}(j,2)-PsiC))/be)^2;
                if ((k1+k2)<=1)
                    AAc(i,1)=AAc(i,1)+1;
                    cnt=cnt+1; % Counts number of AA in ellipse per protein
                    avgphi(cnt)=Proteins{i,1}(j,1); % Sum all phi in ellipse for an average
                    avgpsi(cnt)=Proteins{i,1}(j,2); % Sum all psi in ellipse for an average
                end
            end
        end
    end
    AAc(i,2)=nAA(i)-AAc(i,1); % Counts amino acid residues outside of groups
    AAf(i,:)=AAc(i,:)/nAA(i); % Converts the counts to fractions
end
% Ends program if no AA are in the elliptical cylinder
if cnt==0
    Str=sprintf('\n No Amino Acids in range.\n Re-enter lines 5-9');
    error(Str)
end
nAAsum=sum(nAA); % Total AA counts
% Get an amino acid weighted fraction (proteins have different lengths)
frac1=0.0; frac2=0.0;
for j=1:NumPro
    frac1=frac1+AAf(j,1)*nAA(j)/nAAsum; % Fraction of AA residues inside ellipse
    frac2=frac2+AAf(j,2)*nAA(j)/nAAsum; % Fraction of AA residues outside ellipse
end
AvgPhi=mean(avgphi);
AvgPsi=mean(avgpsi);
str=sprintf('AA fraction inside ellipse %5.3f \nAA fraction outside ellipse %5.3f \nPhiCenter: %7.2f  PhiAvg: %7.2f \nPsiCenter: %7.2f  PsiAvg: %7.2f',frac1,frac2,PhiC,AvgPhi,PsiC,AvgPsi);
disp(str) % Displays quantitative results in command window
str2=sprintf('Total AA:  %4i \tAA In: %4i',sum(aaFlag1),sum(cnt));
disp(str2)
% Input the protein IR spectra, wavenumbers are in the first column, IR spectra are in each successive column
PS=dlmread('protein.txt'); % Input protein spectra
[nsteps nr]=size(PS);
nu=PS(:,1); % Pick off wavenumbers for plotting

% Loop through wavenumber steps
for i=1:nsteps
    R=PS(i,2:NumPro+1)'; % IR intensity of each protein at one wavenumber
    B=inv(AAf'*AAf)*AAf'*R; % This is the fit
    E=R-AAf*B; % Deviations of experimental IR absorbances from fit
    s(i)=((E'*E)/(NumPro-2)); % Est. stand. dev. of protein spectra fits
    a1(i)=B(1,1); % Spectrum inside the ellipse
    ot(i)=B(2,1);
    % Elements of nXn variance/covariance matrix without weighting
    for j=1:cmax
        for k=1:cmax
            A(j,k)=sum(AAf(:,j).*AAf(:,k));
        end
    end
    Ainv=inv(A);
    Er1(i)=Ainv(1,1)*s(i); % Error inside
    Er2(i)=Ainv(2,2)*s(i); % Error outside
end

if plt1==1
    fig1=figure(1); clf% Figure 1: Spectra w/ error
    errorbar(nu,a1,Er1,'.')
    hold on
    plot(nu,a1,'r-','Linewidth',2)
    pnts=[1750 1200;0 0];
    plot(pnts(1,:),pnts(2,:),'k-')
    title('Isolated Spectra With Error')
    axis ( [ 1200 1750 min(a1)*1.1 max(a1)*1.1] )
    xlabel('\nu (cm^{-1})','fontsize',16)
    ylabel('Absorbance','fontsize',16)
    set(gca,'FontSize',14)
    set(gca,'xdir','reverse')
end
A1=a1;
% Integrate the spectra within and outside and scale by fraction of AA
as=sum(a1);
ots=sum(ot);
sall=as+ots;
a1=a1*sall/as;
ot=ot*sall/ots;
a1=frac1*a1;
ot=frac2*ot;

% Figure 2: Plot scaled spectra for inside and outside of box
if plt2==1
    fig2=figure(2); clf
    plot(nu,a1*1,'b',nu,ot,'k',nu,a1+ot,'g')
    title('Spectra In/Out/Toltal')
    %axis ( [ 1200 1750 min(a1) max(a1+ot)*1.05] )
    %legend('Scaled Inner','Scaled Outer','Total','Location','best')
    xlabel('\nu (cm^{-1})','fontsize',16)
    ylabel('Absorbance','fontsize',16)
    set(gca,'FontSize',14)
    set(gca,'xdir','reverse')
end

% Make a contour diagram of angles
[na nc]=size(Angs);
cnt=0; % Intialize count
for i=1:na
    if Angs(i,3)<=upw && Angs(i,3)>=low
        cnt=cnt+1;
        angs(cnt,:)=Angs(i,:);
    end
end

dx=5; % Bin width in angle (in degree)
nsa=365/dx; % Number of steps in angle, nsa must be an integer
% Define bin edges
x=zeros(nsa,1); y=zeros(nsa,1);
for j=1:nsa
    x(j)=-180+dx*(j-1);
    y(j)=-180+dx*(j-1);
end
xc=x+dx/2; yc=y+dx/2;% Define bin centers
% Count angles within 2d angle bins
for j=1:nsa-1
    for k=1:nsa-1
        ch=0;
        for i=1:na
            if ( angs(i,1)>=x(j) && angs(i,1)<x(j+1) && angs(i,2)>=y(k) && angs(i,2)<y(k+1) )
                ch=ch+1;
            end
        end
        ic(j,k)=ch;
    end
end
% Output in irregular form for psi-plot
o2=zeros((nsa-1)*(nsa-1),3);
for j=1:nsa-1
    for k=1:nsa-1
        m=(k-1)*(nsa-1)+j;
        o2(m,1)=xc(j);
        o2(m,2)=yc(k);
        o2(m,3)=ic(j,k);
    end
end

if plt3==1
    fig3=figure(3); clf% Figure 3: Make a contour diagram
    v=0:1.5:30; % Defines contours
    contourf(ic',v,'k')
    hold on
    axis square
    colormap(flipud(autumn))
    colorbar
    % Plotting ellipse
    for i=1:100
        t=(i/100)*2*pi;
        xe(i)=PhiC+ae*cos(t)*cos(phip1)-be*sin(t)*sin(phip1);
        ye(i)=PsiC+ae*cos(t)*sin(phip1)+be*sin(t)*cos(phip1);
    end
    % Plotting psi&phi=0 line
    linn=0:1:72;
    xej=1+(xe-dx/2+180)/dx;
    yej=1+(ye-dx/2+180)/dx;
    plot(36,linn,'.r',linn,36,'.r','MarkerSize',4)
    if ~isempty(AAt)
        plot(Flg1(:,1)/5+36,Flg1(:,2)/5+36,'b<','MarkerSize',4)
    end
    plot(xej,yej,'w','linewidth',2)
    set(gca,'XTickLabel',{'-130','-80','-30','20','70','120','170'})
    set(gca,'YTickLabel',{'-130','-80','-30','20','70','120','170'})
    xlabel('\Phi','fontsize',14)
    ylabel('\Psi','fontsize',14)
    title('Ramachandran Plot')
end

if plt4==1
    fig4=figure(4);
    plot3(angs(:,1),angs(:,2),angs(:,3),'.','MarkerSize',5), grid on
    hold on
    upw1=ones(1,length(xe),1)*upw;
    low1=ones(1,length(xe),1)*low;
    plot3(Flg1(:,1),Flg1(:,2),Flg1(:,3),'r<','MarkerSize',6)
    plot3(xe,ye,upw1,'*g',xe,ye,low1,'*g')
    xlabel('\Phi','fontsize',14)
    ylabel('\Psi','fontsize',14)
    zlabel('\omega','fontsize',14)
    title('3-D Ramachandran Plot')
    legend('AA',AAt{1},'Elliptical Cylinder','Location','NorthEast')
    axis([-180 180 -180 180 low upw])
    hold off
end

if savv==1
    %         sv1=sprintf('%s_SpecError %s.jpg',AAt{1},regn{1});saveas(fig1,sv1)
    %         %sv2=sprintf('All_Spectra.jpg');saveas(fig2,sv2)
    %         sv3=sprintf('%s_RamaPlot %s.jpg',AAt{1},regn{1});saveas(fig3,sv3)
    %     %     sv4=sprintf('Ramachandran3D.jpg');saveas(fig4,sv4)
    %     %Output spectra to txt files
    ac=cat(2,nu,A1');
    %    Sv1=sprintf('%s_Spec_%s.txt',AAt{1},regn{1});
    %     dlmwrite(Sv1,ac,'delimiter','\t');
    %     otc=cat(2,nu,ot');
    %     dlmwrite('result_other_a1_1.txt',otc,'delimiter','\t');
    dlmwrite('BetaSheet.txt',ac);
end