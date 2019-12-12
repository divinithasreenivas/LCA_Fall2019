close all; 
clear all; 
clc;
format short e
% Technology+Intervention matrix columns: Incandescent (by process)
% Linear space with units
% lamp(i)mfg;electricity(MJ);used lamp (i); fuel in kg; copper
% in kg; secondary copper in kg; heat in MJ; glass in kg; 
% waste residue in kg; copper ore in kg; CO2 in kg; SO2 in kg; sand in kg;
% crude oil in kg; copper to soil in kg; 
i_lamp = [1000;-1000;0;0;-5;0;0;-10;0;0;0;0;0;0;0];
i_elec = [0;100;0;-7.5;0;0;200;0;0;0;25;0.15;0;0;0];
i_eol = [0;0;-100;-30;0;0.425;0;0;1.05;0;100;0;0;0;0.025];
i_glass = [0;-100;0;0;0;0;0;1000;0;0;0;0;-1000;0;0];
i_Cu = [0;-10000;0;0;100;0;0;0;0;-1000;0;0;0;0;0];
i_fuel = [0;0;0;1000;0;0;0;0;0;0;200;5;0;-1200;0];

% note we need a unit process for waste residue to landfill since landfill 
% is still an engineered system. But here for simplicity we assume landfill
% = environment.
% Technology+Intervention Matrix: Incandescent
i_Matrix = [i_lamp i_elec i_eol  i_glass i_Cu i_fuel];
i_A = i_Matrix(1:8,:); % Technology Matrix (incandescent)
i_B = i_Matrix(9:15,:); % Intervention Matrix (incandescent)
i_f = [1;1800;-1;0;0;0;0;0]; % final demand vector

% check dimesions of technology matrix
size(i_A);
% note the original technology matrix is not square due to
% 1) for electricity generation we have heat as co-product
% 2) secondary copper is a product not used by the current but apparently
% it is not waste or emission.
% for heat, we will do allocation/partitioning based on dollar value.
% assuming heat is sold at the 1/3 price of electricty per MJ basis
% ratio of electricity revenue vs. heat revenue
r=100/(200*1/3);
% partition coefficient
lmd=r/(1+r);
i_elec_n=i_elec*lmd;
i_elec_n(2)=i_elec(2);
i_elec_n(7)=0;
i_elec_n;
i_heat=i_elec*(1-lmd);
i_heat(2)=0;
i_heat(7)=i_elec(7);
i_heat;
% now for secondary copper, it is assumed that it can be purified. This
% process consumes only 20% of energy as in the primary copper production.
% It is also assumed the yield is 90%.
i_2nd_Cu=[0;-2000;0;0;90;-100;0;0;10;0;0;0;0;0;0];
% now reassemble the matrix
names_ls = {'lamp(i)';'electricity(MJ)';'used lamp (i)'; 'fuel in kg'; 'copper in kg'; 'secondary copper in kg'; 'heat in MJ'; 'glass in kg'};
names_p ={'lamp(i)mfg';'electricity(MJ)';'treat used lamp (i)';'fuel production';'copper production';'processing secondary Cu';'heat production';'glass production'};
i_Matrix_n=[i_lamp i_elec_n i_eol i_fuel i_Cu  i_2nd_Cu i_heat i_glass];
i_A_n = i_Matrix_n(1:8,:); % Technology Matrix (incandescent)
i_B_n = i_Matrix_n(9:15,:); % Intervention Matrix (incandescent)
% now we have a squre matrix
[mA,nA]=size(i_A_n);
[mB,nB]=size(i_B_n);
% check if technology is invertible
det(i_A_n);
% scaling vector
i_s=inv(i_A_n)*i_f;
% inventory
i_g=i_B_n*i_s;



%===============================================================

% ++++++++++++  Sensitivity Analysis ++++++++++++++

% using analytical equations
% sensitivity to entries in technology matrix
i_lambda = i_B_n*(inv(i_A_n)); % Incandescent uncertainty parameter


for k = 1:nB;
    for i = 1:mA
        for j = 1:nA
            if k == 1; % Entry 1 corresponds to Waste Residue
              i_waste(i,j) = -(i_A_n(i,j))/i_g(k) * i_lambda(k,i)*i_s(j);
                           
            elseif k == 2; % Entry 2 corresponds to Cu Ore
              i_CuOre(i,j) = -(i_A_n(i,j))/i_g(k) * i_lambda(k,i)*i_s(j);
                
             elseif k == 3; % Entry 3 corresponds to CO2
                i_CO2(i,j) = -(i_A_n(i,j))/i_g(k) * i_lambda(k,i)*i_s(j);
             
            elseif k == 4; % Entry 4 corresponds to SO2
               i_SO2(i,j) = -(i_A_n(i,j))/i_g(k) * i_lambda(k,i)*i_s(j);
  
            elseif k == 5; % Entry 5 corresponds to sand
              i_sand(i,j) = -(i_A_n(i,j))/i_g(k) * i_lambda(k,i)*i_s(j);
 
            elseif k == 6; % Entry 6 corresponds to crude oil
              i_crude(i,j) = -(i_A_n(i,j))/i_g(k) * i_lambda(k,i)*i_s(j);
                          
            elseif k == 7; % Entry 7 corresponds to Cu to Soil
              i_CuSoil(i,j) = -(i_A_n(i,j))/i_g(k) * i_lambda(k,i)*i_s(j);

            end
        end
    end
end


% sensitivity to entries in intervention matrix


for k = 1:nB;
    for i = 1:mA
        for j = 1:nA
            if k == 1; % Entry 1 corresponds to Waste Residue
              if i== k
                  i_waste_B(i,j) = i_B_n(i,j)/i_g(k)*i_s(j);
              else
                  i_waste_B(i,j) = 0;
              end
                           
            elseif k == 2; % Entry 2 corresponds to Cu Ore
              if i== k
                  i_CuOre_B(i,j) = i_B_n(i,j)/i_g(k)*i_s(j);
              else
                  i_CuOre_B(i,j) = 0;
              end
                
             elseif k == 3; % Entry 3 corresponds to CO2
              if i== k
                  i_CO2_B(i,j) = i_B_n(i,j)/i_g(k)*i_s(j);
              else
                  i_CO2_B(i,j) = 0;
              end
             
            elseif k == 4; % Entry 4 corresponds to SO2
              if i== k
                  i_SO2_B(i,j) = i_B_n(i,j)/i_g(k)*i_s(j);
              else
                  i_SO2_B(i,j) = 0;
              end
  
            elseif k == 5; % Entry 5 corresponds to sand
              if i== k
                  i_sand_B(i,j) = i_B_n(i,j)/i_g(k)*i_s(j);
              else
                  i_sand_B(i,j) = 0;
              end
 
            elseif k == 6; % Entry 6 corresponds to crude oil
              if i== k
                  i_crude_B(i,j) = i_B_n(i,j)/i_g(k)*i_s(j);
              else
                  i_crude_B(i,j) = 0;
              end                          
            elseif k == 7; % Entry 7 corresponds to Cu to Soil
              if i== k
                  i_CuSoil_B(i,j) = i_B_n(i,j)/i_g(k)*i_s(j);
              else
                  i_CuSoil_B(i,j) = 0;
              end

            end
        end
    end
end




% using numerical approach
% for each entriy in technology matrix A, we introduce +/- delta% change
delta=0.01;

% +delta% disturbance:


for i = 1:mA
   for j = 1:nA
                i_A_dlt=i_A_n;
                i_A_dlt(i,j)=i_A_n(i,j)*(1+delta);
                i_g_dlt=i_B_n*inv(i_A_dlt)*i_f;
        for k = 1:nB;
            if k == 1; % Entry 1 corresponds to Waste Residue
                i_waste_p(i,j) = (i_g_dlt(k)-i_g(k))/i_g(k)/delta;
                           
            elseif k == 2; % Entry 2 corresponds to Cu Ore
              i_CuOre_p(i,j) = (i_g_dlt(k)-i_g(k))/i_g(k)/delta;
                
            elseif k == 3; % Entry 3 corresponds to CO2
                i_CO2_p(i,j) = (i_g_dlt(k)-i_g(k))/i_g(k)/delta;
             
            elseif k == 4; % Entry 4 corresponds to SO2
               i_SO2_p(i,j) = (i_g_dlt(k)-i_g(k))/i_g(k)/delta;
  
            elseif k == 5; % Entry 5 corresponds to sand
              i_sand_p(i,j) = (i_g_dlt(k)-i_g(k))/i_g(k)/delta;
 
            elseif k == 6; % Entry 6 corresponds to crude oil
              i_crude_p(i,j) = (i_g_dlt(k)-i_g(k))/i_g(k)/delta;
                          
            elseif k == 7; % Entry 7 corresponds to Cu to Soil
              i_CuSoil_p(i,j) = (i_g_dlt(k)-i_g(k))/i_g(k)/delta;

            end
        end
    end
end

% -delta% disturbance:


for i = 1:mA
   for j = 1:nA
                i_A_dlt=i_A_n;
                i_A_dlt(i,j)=i_A_n(i,j)*(1-delta);
                i_g_dlt=i_B_n*inv(i_A_dlt)*i_f;
        for k = 1:nB;
            if k == 1; % Entry 1 corresponds to Waste Residue
                i_waste_m(i,j) = (i_g_dlt(k)-i_g(k))/i_g(k)/delta;
                           
            elseif k == 2; % Entry 2 corresponds to Cu Ore
              i_CuOre_m(i,j) = (i_g_dlt(k)-i_g(k))/i_g(k)/delta;
                
            elseif k == 3; % Entry 3 corresponds to CO2
                i_CO2_m(i,j) = (i_g_dlt(k)-i_g(k))/i_g(k)/delta;
             
            elseif k == 4; % Entry 4 corresponds to SO2
               i_SO2_m(i,j) = (i_g_dlt(k)-i_g(k))/i_g(k)/delta;
  
            elseif k == 5; % Entry 5 corresponds to sand
              i_sand_m(i,j) = (i_g_dlt(k)-i_g(k))/i_g(k)/delta;
 
            elseif k == 6; % Entry 6 corresponds to crude oil
              i_crude_m(i,j) = (i_g_dlt(k)-i_g(k))/i_g(k)/delta;
                          
            elseif k == 7; % Entry 7 corresponds to Cu to Soil
              i_CuSoil_m(i,j) = (i_g_dlt(k)-i_g(k))/i_g(k)/delta;

            end
        end
    end
end

% Now let's make a tornado plot to show sensitivity
% Take CO2 emission as the impact of interest
% note we need to remove entry corresponding the product coming out a unit
% process
% the technology is set up in a way that all these entries are on the
% diagonal

for i=1:mA
    i_CO2_p(i,i)=0;
    i_CO2_m(i,i)=0;
end;

nargin=length(names_p);
% convert i_CO2_p matrix to a column vector so we can sort
v_i_CO2_p=i_CO2_p(:);
v_i_CO2_m=i_CO2_m(:);
[v_i_CO2_p_sort,ind_CO2_p]=sort(abs(v_i_CO2_p/delta),'descend');
v_i_CO2_p_sort=v_i_CO2_p(ind_CO2_p);
v_i_CO2_m_sort=v_i_CO2_m(ind_CO2_p);
% here we consider the top 5 entries in technology matrix which have
% largest effects on CO2 emission
n_s=5;
v_i_CO2_ps=v_i_CO2_p_sort(1:n_s)
v_i_CO2_ms=v_i_CO2_m_sort(1:n_s)

% finding names of the top n_s entries
for i=1:n_s
    ind_t=ind_CO2_p(i);
     np=floor(ind_t/mA)+1; % index of unit process
     nr=mod(ind_t,mA); % position of entry in the unit process
     % names of variable has two parts: name of unit process + entry in
     % linear space
     name1=names_p(np);
     if nr==0
         name2=names_ls(mA);
     else
         name2=names_ls(nr);
     end;
     
     names(i)=strcat('amount of',{' '},name2, ' in',{' '},name1);
   
end;  
Objective_high_sum = v_i_CO2_ps';
Objective_low_sum = v_i_CO2_ms';
Objective_base_value = 0;
[Objective_low_sum_t,ind]=sort(abs(Objective_low_sum),'ascend');
Objective_low_sum=Objective_low_sum(ind);
Objective_high_sum=Objective_high_sum(ind);
names_Objective=names;
figure (1)
h = barh(Objective_high_sum);
hold on
xmin=min([min(Objective_low_sum),min(Objective_high_sum)]);
xmax=max([max(Objective_low_sum),max(Objective_high_sum)]);
xlim([0.975*xmin 1.025*xmax])
barh(Objective_low_sum,'r')
bh = get(h,'BaseLine');
set(bh,'BaseValue',Objective_base_value);
title('Sensitivities')
if nargin > 1
    set(gca,'yticklabel',names)
    set(gca,'Ytick',[1:length(names)],'YTickLabel',[1:length(names)])
    set(gca,'yticklabel',names_Objective)
end
xlabel('Percentage Changes on output due to -% red and +% blue change in input')


% Monte-Carlo Simualtion

% Here we select CO2 emission as the emission of interest. 
% Note all the entries in the technology and intervention matrix can be
% stochastic. Here we just pick one entry to deomnstrate the approach.
% for incandescent lamp manufacturing, it was originally given that 5 kg of
% copper is needed to make 1000 lamps. Now let's assume that copper
% consumption follows a certain distribution with 5 as the mean. In
% Ecoinvent database, there are many inputs that follow lognormal
% distribution. Here we assume copper consumption also follows lognormal
% distribution.
% specify parameters for the distribution
% two are needed for lognormal

% now we need a random number generator. Fortunately Matlab provides
% functions for many commonly used distribution. And you can always search
% online.
% For Monte Carlo simulation, we need to specify sample size. Usually the
% larger the better.
N_s=15000;
% number of bins for histgram
M_b=25;
% Note in matlab the "lognrnd" command requires input parameters as
% the mean and standard deviation of the corresponding normal distribution.
% While in Ecoinvent the median and geometric variance are given.
% Assume from Ecoinvent we have
med=5;
var_g=3;
% then for the corresponding normal distribution, we have:
mu = log(med);
sigma = log(sqrt(var_g));
% generating random numbers
R=lognrnd(mu,sigma,N_s,1);
% for each sample, we need to calculate a new CO2 emission
% the entry we want to change is A(5,1).
i_A_mc=i_A_n;
for i=1:N_s
    i_A_mc(5,1)=-R(i);
    i_g_mc=i_B_n*inv(i_A_mc)*i_f;
    carbon(i)=i_g_mc(3);
end;

mu_c=mean(carbon);
var_c=var(carbon);
i_g(3);
figure(2);
hist(carbon,M_b);

% another commonly used distribution is triangular distribution
% For triangular distribution,we need three parameters, the range
% [min,max], and the mode (most probable value).
% assume for electricity generation, fuel needed follows a triangular
% distribution.
min=6;
max=8;
md=7.5;
% note in the original case 7.5 kg fuel is needed to generate 100 MJ
% electricity and 200 MJ heat, while emitting 25 kg CO2 and 0.15 kg SO2. As
% we change the fuel needed, the CO2 emission and SO2 emission should
% change accordingly. Here we assume there exists a linear relationship.

% below we show how to use Matlab's function for generating uniformly distributed
% pseudorandom numbers to realize a random numer generator for any given
% distribution

R_t=zeros(N_s,1);
for i=1:N_s
z=rand;
if sqrt(z*(max-min)*(md-min))+min < md
    R_t(i)=sqrt(z*(max-min)*(md-min))+min;
else
    R_t(i)=max-sqrt((1-z)*(max-min)*(max-md));
end
end 
%hist(R_t);

% for each sample, we need to calculate a new CO2 emission
% the entries we want to change are A(4,2) and A(4,7); B(3,2), B(4,2), B(3,7), and B(4,7).
i_A_mc_t=i_A_n;
i_B_mc_t=i_B_n;
for i=1:N_s
    i_A_mc_t(4,2)=-R_t(i)*lmd;
    i_A_mc_t(4,7)=-R_t(i)*(1-lmd);
    i_B_mc_t(3,2)=i_B_n(3,2)*R_t(i)/md;
    i_B_mc_t(4,2)=i_B_n(4,2)*R_t(i)/md;
    i_B_mc_t(3,7)=i_B_n(3,7)*R_t(i)/md;
    i_B_mc_t(4,7)=i_B_n(4,7)*R_t(i)/md;   
    i_g_mc_t=i_B_mc_t*inv(i_A_mc_t)*i_f;
    carbon_t(i)=i_g_mc_t(3);
end;
figure(3);
hist(carbon_t,M_b);

% now we check the combined effects on carbon footprint. Assume the two
% random variables are independent. Note this may not always the case.
% Joint distributions are needed but usaully these are difficult to find.
% Also, we may have many random varibales invovled.
i_A_mc_cm=i_A_n;
i_B_mc_cm=i_B_n;
for i=1:N_s
    i_A_mc_cm(5,1)=-R(i);
    i_A_mc_cm(4,2)=-R_t(i)*lmd;
    i_A_mc_cm(4,7)=-R_t(i)*(1-lmd);
    i_B_mc_cm(3,2)=i_B_n(3,2)*R_t(i)/md;
    i_B_mc_cm(4,2)=i_B_n(4,2)*R_t(i)/md;
    i_B_mc_cm(3,7)=i_B_n(3,7)*R_t(i)/md;
    i_B_mc_cm(4,7)=i_B_n(4,7)*R_t(i)/md; 
    
    i_g_mc_cm=i_B_mc_cm*inv(i_A_mc_cm)*i_f;
    
    carbon_cm(i)=i_g_mc_cm(3);
end;
figure(4);
hist(carbon_cm,M_b);
mean_cm=mean(carbon_cm)
var_cm=var(carbon_cm);

% determine the x-th percentile of the data. 
% usually 25% and 75% or 5% and 95% are used.
y25=prctile(carbon_cm, 25)
y75=prctile(carbon_cm, 75)

% determine the confidence interval
pd = fitdist(carbon_cm','Normal');
ci = paramci(pd,'Alpha',.05)
