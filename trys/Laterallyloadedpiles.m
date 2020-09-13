%% Comments
clc; clear all;

% This analysis is to find the influence of long and short piles.
% The value of hh is standardised at 40/200 = 0.2


%Loop for range of pile diameter
        LP = 40;                                   %Length of Pile in metres
        DP = 1;                                    %Diameter of Pile in metres
        NN = 2000;                                  %Number of elements underground
        hh = LP/NN;                                %Length of elements (epsilon/m)
        tt = 0.06;                                 %Thickness of pile - 6cm (0.06m)
        EE = 210*10^6;                             %Young modulus of steel (kPa-kN/m2)
        II = (pi/4)*(((DP/2)^4)-((DP/2)-tt)^4);    %Second moment of area (m4)
        EI = EE*II;                                %EI for pile (kNm2)

        %Input Data
        FF = 1000;                               %UDL at the top of pile (kN/m)

        %Stiffness calculation
        ES = 50000;                               %Young modulus for dense sand (kN/m2)
        Ef = 1;                                    %Efact (0.8 for clay/ for sand)
        KS = ES*Ef;                                %Soil stiffness (kN/m2)


    %% WITH ROTATIONAL STIFFNESS
    N = 100; %no of slices
    c = 100/0.0100;  %slope of t-z curve
    KR = 0;
        for ii=1:N
            x(ii)=DP/(4*N)*(1+(ii-1)*2);
            alpha(ii)= acos((ii-1)/N)-acos(ii/N);
            C(ii)=alpha(ii)*DP/2;
            KR=KR+2*c*x(ii)*x(ii)*C(ii);
        end

    AA = zeros(NN+5,NN+5);                     
        for ii=3:NN+3                              %For the 3rd to NN+4th nodes
            if ii ==3
                AA(ii-2,ii-2:ii+2)=[1,-4,6,-4,1] + KS*hh^4/EI*[0,0,1,0,0];
            else
                AA(ii-2,ii-2:ii+2)=[1,-4,6,-4,1] + KS*hh^4/EI*[0,0,1,0,0] - ...
                    KR*hh^2/EI*[0,1,-2,1,0];  
            end
        end

    %Boundary conditions
    AA(NN+2,1:5)=[0,(-hh*KR+2*EI),-4*EI,(hh*KR+2*EI),0];                                       %Top rotational stiffness relationship
    AA(NN+3,1:5)=1/2*[-EI/hh^3, 2*EI/hh^3+KR/hh,0, -2*EI/hh^3-KR/hh, EI/hh^3];                 %Top shear
    AA(NN+4,NN+1:NN+5)=[-1/12,4/3,-5/2,4/3,-1/12];                                             %Bottom  moment (=0)
    AA(NN+5,NN+1:NN+5)=1/2*[-EI/hh^3, 2*EI/hh^3+KR/hh,0, -2*EI/hh^3-KR/hh, EI/hh^3];           %Bottom shear   (=0)

    %bb terms
    bb = zeros(NN+5,1);
    bb(1,1) = FF*hh^4/EI;                    %Top(node 3) load  
    bb(NN+2,1) = 0;                           %Top(node 3) moment
    bb(NN+3,1) = FF;                      %Top(node 3) shear

    %Displacement/deflection (metres)
    ww = AA\bb; 

%% WITHOUT ROTATIONAL STIFFNESS
    % AA matrix term
    AA = zeros(NN+5,NN+5);                     
    for ii=3:NN+3                              %For the 3rd to NN+4th nodes
        CS = KS*(hh^4)/EI;                     %Stiffness of soil
        AA(ii-2,ii-2:ii+2)=[1,-4,6,-4,1] + CS*[0,0,1,0,0];
    end

    % Boundary conditions
    AA(NN+2,1:5)=[-1/12,4/3,-5/2,4/3,-1/12];        %Top  moment 
    AA(NN+3,1:5)=[-1/2,1,0,-1,1/2];                 %Top shear
    AA(NN+4,NN+1:NN+5)=[-1/12,4/3,-5/2,4/3,-1/12];  %Bottom  moment (=0)
    AA(NN+5,NN+1:NN+5)=[-1/2,1,0,-1,1/2];           %Bottom shear  (=0)

    % bb terms 
    bb = zeros(NN+5,1);
    bb(1,1) = FF*hh^4/EI;                        %Top(node 3) load  
    bb(NN+2,1)= (0)*hh^2/EI;                     %Top(node 3) moment
    bb(NN+3,1)= (FF) *hh^3/EI;                   %Top(node 3) shear


%%
    
hFig = figure(1);
set(hFig, 'Position', [0 0 700 900])

xx = linspace(0,-LP,NN+1);
w1 = ww*10^3;

subplot (5,1,1); ylim([-LP,0]); 
hold all;  
plot(w1(3:end-2),xx,'b-');
legend('With Rotational Resistance','Without Rotational Resistance','location','Southeast')
set(gca, 'Fontsize', 12)
xlabel('Displacement,mm'); ylabel('depth,m'); 
  
% subplot (5,1,2); ylim([-LP,0]);
% hold all; 
% plot(MM(3:end-2),xx,'b-');
% legend('With Rotational Resistance','Without Rotational Resistance','location','southwest')
% set(gca, 'Fontsize', 12)
% xlabel('Bending Moment,kNm'); ylabel('depth,m');
% 
% subplot (5,1,3); ylim([-LP,0]);
% hold all;  
% plot(SS(3:end-2),xx,'b-');
% legend('With Rotational Resistance','Without Rotational Resistance','location','Southeast')
% set(gca, 'Fontsize', 12)
% xlabel('Shear Force,kN'); ylabel('depth,m');
% 
% subplot (5,1,4); ylim([-LP,0]);
% hold all; 
% plot(RR(3:end-2),xx,'b-');
% legend('With Rotational Resistance','Without Rotational Resistance','location','Southeast')
% set(gca, 'Fontsize', 12)
% xlabel('Rotation,Rad'); ylabel('depth,m');
%  
% subplot (5,1,5); ylim([-LP,0]);
% hold all;
% plot(PP(3:end-2),xx,'b-');
% legend('With Rotational Resistance','Without Rotational Resistance','location','Southeast')
% set(gca, 'Fontsize', 12) 
% xlabel('Pressure,kPa'); ylabel('depth,m')