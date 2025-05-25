clear;clc;
%Calculate the value of phase reestimation
%Determine the number of bins by yourself
numbin=256;
%Import the standard template of pulsar contours
load('profile_value.mat')

%Calculate the period estimation verification of the actual XPNAV-1 observation data
%Set the number of observation groups according to the actual situation
Z=121;
for j= 1:Z
    %Import the flag bit information required for preprocessing data
    load('flag.mat')
    %Import the actual period of the JB Observatory to facilitate the subsequent calculation of error values
    load('JB.mat')
    %Initialize the matrix variables
    TOA1=[];
    TOA2=[];
    TOA=[];
    %Automatically import the data corrected by Einstein one by one
    str1='tdb_time11.txt';
    num=flag(j,1);
    str11='11';
    str12=num2str(num);
    str1=strrep(str1, str11, str12);
    a1=importdata(str1);
    %Perform MJD time conversion
    DATA=MJD(a1);
    %Automatically import the shapiro delay correction information one by one
    str2='r_shapiroDE432s11.mat';
    str21='11';
    str22=num2str(num);
    str2=strrep(str2, str21, str22);
    load(str2);
    %Automatically import the roemer delay correction information one by one
    str3='r_roemerDE432s11.mat';
    str31='11';
    str32=num2str(num);
    str3=strrep(str3, str31, str32);
    load(str3);
    %Perform delay information correction
    str41='r_roemerDE432s';
    str42=num2str(num);
    str4=strcat(str41, str42);
    str51='r_shapiroDE432s';
    str52=num2str(num);
    str5=strcat(str51, str52);
    DATA=DATA+eval(str4)+eval(str5);
    str61='';
    str62=num2str(flag(j,2));
    str63='.txt';
    str6=strcat(+str61, str62,str63);
    a2=importdata(str6);
    %Extract the energy carried by photons
    DATA(:,2)=a2(flag(j,3):flag(j,4),2);
    %Carry out energy screening
    B = DATA(DATA(:,2) >= 500 & DATA(:,2) <= 4000, :);
    DATA=[];
    DATA=B;
    [m,n] = size(DATA);
    TOA(:,1)=DATA(:,1)-DATA(1,1)*ones(m,1);
    %The first-order adaptive HTO value is designed based on the relationship between the number of photons and HTO
    rank = HTO(m);
    bin_num = BIN;
    chi_range = 1;
    %Lower bound of the search cycle
    eP_lowerbound = P1; 
    %Upper bound of the search cycle
    eP_upperbound = P2;
    %Search step size
    eP_step = 0.0000000001;

    for P = eP_lowerbound:eP_step:eP_upperbound
        Tb = P/bin_num;
        %Perform the photon epoch folding operation
        profile = Epoch_folding(bin_num,P,Tb,m,TOA);

        %Calculate the distance between the standard profile and the template profile, and move the template profile to correspond to the standard profile
        profile_std=move_phase(profile_value,profile,numbin);

        %Calculate Δfai for each track based on the slope of the template
        xx=linspace(0, 2*pi, (numbin+1))';
        for w=1:numbin
            faib(w,1)=sovle_k(xx(w,1),xx(w+1,1),profile_std(w,1),profile_std(w+1,1));
        end
        %Calculate the value of Z
        N=m;
        chi2(chi_range,1) = Z2(rank,N,bin_num,profile,faib);
        chi_range =  chi_range+1;
    end
    %Determine the maximum value of Z and its position index
    [maxchi,p_pos] = max(chi2);
    P = (p_pos-1)*eP_step+eP_lowerbound;
    error(j,1)=P;
    error(j,2)=1e9*abs(P-JB(j,1))
end