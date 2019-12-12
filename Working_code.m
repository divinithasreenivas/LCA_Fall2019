close all; 
clear all; 
clc;
format short e
%% Setup to Import B
opts = spreadsheetImportOptions("NumVariables", 41);

% Specify sheet and range
opts.Sheet = "Invertible Matrix";
opts.DataRange = "C45:AQ47";

% Specify column names and types
opts.VariableNames = ["Transportschoolbusdieselpowered", "Dieselatrefinery", "Bitumenatrefinery", "Gasolineatrefinery", "Keroseneatrefinery", "Liquefiedpetroleumgasatrefinery", "Petroleumcokeatrefinery", "Petroleumrefiningcoproductatrefinery", "Petroleumrefiningatrefinery", "Refinerygasatrefinery", "Residualfueloilatrefinery", "Transportbargeaveragefuelmix", "Transportcombinationtruckaveragefuelmix", "Transportoceanfreighteraveragefuelmix", "Transportpipelineunspecifiedpetroleumproducts", "Transporttraindieselpowered", "Crudeoilextracted", "Dieselcombustedinindustrialequipment", "Gasolinecombustedinequipment", "Naturalgascombustedinindustrialboiler", "Residualfueloilcombustedinindustrialboiler", "Liquefiedpetroleumgascombustedinindustrialboiler", "Transportpipelinenaturalgas", "Transportcombinationtruckdieselpowered", "Naturalgascombustedinindustrialequipment", "Naturalgasprocessedatplant", "Transportbargedieselpowered", "Transportbargeresidualfueloilpowered", "Transportoceanfreighterdieselpowered", "Transportoceanfreighterresidualfueloilpowered", "Electricityatgrid", "Electricitybiomassatpowerplant", "Electricitybituminouscoalatpowerplant", "Bituminouscoalcombustedinindustrialboiler", "Bituminouscoalatmine", "Dieselcombustedinindustrialboiler", "Electricitynaturalgasatpowerplant", "Electricitynuclearatpowerplant", "Fuelgradeuraniumatregionalstorage", "Electricityresidualfueloilatpowerplant", "Naturalgasextracted"];
opts.SelectedVariableNames = ["Transportschoolbusdieselpowered", "Dieselatrefinery", "Bitumenatrefinery", "Gasolineatrefinery", "Keroseneatrefinery", "Liquefiedpetroleumgasatrefinery", "Petroleumcokeatrefinery", "Petroleumrefiningcoproductatrefinery", "Petroleumrefiningatrefinery", "Refinerygasatrefinery", "Residualfueloilatrefinery", "Transportbargeaveragefuelmix", "Transportcombinationtruckaveragefuelmix", "Transportoceanfreighteraveragefuelmix", "Transportpipelineunspecifiedpetroleumproducts", "Transporttraindieselpowered", "Crudeoilextracted", "Dieselcombustedinindustrialequipment", "Gasolinecombustedinequipment", "Naturalgascombustedinindustrialboiler", "Residualfueloilcombustedinindustrialboiler", "Liquefiedpetroleumgascombustedinindustrialboiler", "Transportpipelinenaturalgas", "Transportcombinationtruckdieselpowered", "Naturalgascombustedinindustrialequipment", "Naturalgasprocessedatplant", "Transportbargedieselpowered", "Transportbargeresidualfueloilpowered", "Transportoceanfreighterdieselpowered", "Transportoceanfreighterresidualfueloilpowered", "Electricityatgrid", "Electricitybiomassatpowerplant", "Electricitybituminouscoalatpowerplant", "Bituminouscoalcombustedinindustrialboiler", "Bituminouscoalatmine", "Dieselcombustedinindustrialboiler", "Electricitynaturalgasatpowerplant", "Electricitynuclearatpowerplant", "Fuelgradeuraniumatregionalstorage", "Electricityresidualfueloilatpowerplant", "Naturalgasextracted"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Import the data from path 
mat_b = readtable("\\nas01.itap.purdue.edu\puhome\desktop\Copy of DIV_HW#6_newchanges - backup.xlsx", opts, "UseExcel", false);

%% Convert to output type
mat_b = table2array(mat_b);

%% Clear temporary variables
clear opts


%% Clear temporary variables
clear opts

opts = spreadsheetImportOptions("NumVariables", 1);

% Specify sheet and range
opts.Sheet = "Invertible Matrix";
opts.DataRange = "AS2:AS42";

% Specify column names and types
opts.VariableNames = "FinalDemand";
opts.SelectedVariableNames = "FinalDemand";
opts.VariableTypes = "double";

% Import the data
mat_f = readtable("\\nas01.itap.purdue.edu\puhome\desktop\Copy of DIV_HW#6_newchanges - backup.xlsx", opts, "UseExcel", false);

%% Convert to output type
mat_f = table2array(mat_f);% importing final demand data

%% Clear temporary variables
clear opts

opts = spreadsheetImportOptions("NumVariables", 41);

% Specify sheet and range
opts.Sheet = "Invertible Matrix";
opts.DataRange = "C2:AQ42";

% Specify column names and types
opts.VariableNames = ["Transportschoolbusdieselpowered", "Dieselatrefinery", "Bitumenatrefinery", "Gasolineatrefinery", "Keroseneatrefinery", "Liquefiedpetroleumgasatrefinery", "Petroleumcokeatrefinery", "Petroleumrefiningcoproductatrefinery", "Petroleumrefiningatrefinery", "Refinerygasatrefinery", "Residualfueloilatrefinery", "Transportbargeaveragefuelmix", "Transportcombinationtruckaveragefuelmix", "Transportoceanfreighteraveragefuelmix", "Transportpipelineunspecifiedpetroleumproducts", "Transporttraindieselpowered", "Crudeoilextracted", "Dieselcombustedinindustrialequipment", "Gasolinecombustedinequipment", "Naturalgascombustedinindustrialboiler", "Residualfueloilcombustedinindustrialboiler", "Liquefiedpetroleumgascombustedinindustrialboiler", "Transportpipelinenaturalgas", "Transportcombinationtruckdieselpowered", "Naturalgascombustedinindustrialequipment", "Naturalgasprocessedatplant", "Transportbargedieselpowered", "Transportbargeresidualfueloilpowered", "Transportoceanfreighterdieselpowered", "Transportoceanfreighterresidualfueloilpowered", "Electricityatgrid", "Electricitybiomassatpowerplant", "Electricitybituminouscoalatpowerplant", "Bituminouscoalcombustedinindustrialboiler", "Bituminouscoalatmine", "Dieselcombustedinindustrialboiler", "Electricitynaturalgasatpowerplant", "Electricitynuclearatpowerplant", "Fuelgradeuraniumatregionalstorage", "Electricityresidualfueloilatpowerplant", "Naturalgasextracted"];
opts.SelectedVariableNames = ["Transportschoolbusdieselpowered", "Dieselatrefinery", "Bitumenatrefinery", "Gasolineatrefinery", "Keroseneatrefinery", "Liquefiedpetroleumgasatrefinery", "Petroleumcokeatrefinery", "Petroleumrefiningcoproductatrefinery", "Petroleumrefiningatrefinery", "Refinerygasatrefinery", "Residualfueloilatrefinery", "Transportbargeaveragefuelmix", "Transportcombinationtruckaveragefuelmix", "Transportoceanfreighteraveragefuelmix", "Transportpipelineunspecifiedpetroleumproducts", "Transporttraindieselpowered", "Crudeoilextracted", "Dieselcombustedinindustrialequipment", "Gasolinecombustedinequipment", "Naturalgascombustedinindustrialboiler", "Residualfueloilcombustedinindustrialboiler", "Liquefiedpetroleumgascombustedinindustrialboiler", "Transportpipelinenaturalgas", "Transportcombinationtruckdieselpowered", "Naturalgascombustedinindustrialequipment", "Naturalgasprocessedatplant", "Transportbargedieselpowered", "Transportbargeresidualfueloilpowered", "Transportoceanfreighterdieselpowered", "Transportoceanfreighterresidualfueloilpowered", "Electricityatgrid", "Electricitybiomassatpowerplant", "Electricitybituminouscoalatpowerplant", "Bituminouscoalcombustedinindustrialboiler", "Bituminouscoalatmine", "Dieselcombustedinindustrialboiler", "Electricitynaturalgasatpowerplant", "Electricitynuclearatpowerplant", "Fuelgradeuraniumatregionalstorage", "Electricityresidualfueloilatpowerplant", "Naturalgasextracted"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Import the data
mat_a = readtable("\\nas01.itap.purdue.edu\puhome\desktop\Copy of DIV_HW#6_newchanges - backup.xlsx", opts, "UseExcel", false);

%% Convert to output type
mat_a = table2array(mat_a);

%% Clear temporary variables
clear opts


%% Setup the Import Options
opts = spreadsheetImportOptions("NumVariables", 1);

% Specify sheet and range
opts.Sheet = "Invertible Matrix";
opts.DataRange = "A2:A42";

% Specify column names and types
opts.VariableNames = "LinearSpace";
opts.VariableTypes = "string";
opts = setvaropts(opts, 1, "WhitespaceRule", "preserve");
opts = setvaropts(opts, 1, "EmptyFieldRule", "auto");

% Import the data
linear_space = readtable("\\nas01.itap.purdue.edu\puhome\desktop\Copy of DIV_HW#6_newchanges - backup.xlsx", opts, "UseExcel", false);

%% Convert to output type
linear_space = table2array(linear_space);

%% Clear temporary variables
clear  opts




%% Setup the Import Options
opts = spreadsheetImportOptions("NumVariables", 41);

% Specify sheet and range
opts.Sheet = "Invertible Matrix";
opts.DataRange = "C1:AQ1";

% Specify column names and types
opts.VariableNames = ["Transportschoolbusdieselpowered", "Dieselatrefinery", "Bitumenatrefinery", "Gasolineatrefinery", "Keroseneatrefinery", "Liquefiedpetroleumgasatrefinery", "Petroleumcokeatrefinery", "Petroleumrefiningcoproductatrefinery", "Petroleumrefiningatrefinery", "Refinerygasatrefinery", "Residualfueloilatrefinery", "Transportbargeaveragefuelmix", "Transportcombinationtruckaveragefuelmix", "Transportoceanfreighteraveragefuelmix", "Transportpipelineunspecifiedpetroleumproducts", "Transporttraindieselpowered", "Crudeoilextracted", "Dieselcombustedinindustrialequipment", "Gasolinecombustedinequipment", "Naturalgascombustedinindustrialboiler", "Residualfueloilcombustedinindustrialboiler", "Liquefiedpetroleumgascombustedinindustrialboiler", "Transportpipelinenaturalgas", "Transportcombinationtruckdieselpowered", "Naturalgascombustedinindustrialequipment", "Naturalgasprocessedatplant", "Transportbargedieselpowered", "Transportbargeresidualfueloilpowered", "Transportoceanfreighterdieselpowered", "Transportoceanfreighterresidualfueloilpowered", "Electricityatgrid", "Electricitybiomassatpowerplant", "Electricitybituminouscoalatpowerplant", "Bituminouscoalcombustedinindustrialboiler", "Bituminouscoalatmine", "Dieselcombustedinindustrialboiler", "Electricitynaturalgasatpowerplant", "Electricitynuclearatpowerplant", "Fuelgradeuraniumatregionalstorage", "Electricityresidualfueloilatpowerplant", "Naturalgasextracted"];
opts.VariableTypes = ["string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string"];
opts = setvaropts(opts, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41], "WhitespaceRule", "preserve");
opts = setvaropts(opts, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41], "EmptyFieldRule", "auto");

% Import the data
process = readtable("\\nas01.itap.purdue.edu\puhome\desktop\Copy of DIV_HW#6_newchanges - backup.xlsx", opts, "UseExcel", false);

%% Convert to output type
process = table2array(process);% importing process names

%% Clear temporary variables
clear  opts


%assign imported data to variables used in the code
i_A = mat_a; % Technology Matrix (incandescent)
i_B = mat_b; % Intervention Matrix (incandescent)
i_f = mat_f; % final demand vector
names_ls = linear_space ;
names_p = process;

% Since th matrix has aldedy been allocted and converted to an intervtible
% matrix, we asssign it to the variable 
i_A_n =i_A; % Technology Matrix (incandescent)
i_B_n = i_B; % Intervention Matrix (incandescent)
% now we have a squre matrix
[mA,nA]=size(i_A_n);
[mB,nB]=size(i_B_n);
% scaling vector
i_s=inv(i_A_n)*i_f;
% inventory
i_g=i_B_n*i_s;

inverseA= inv(i_A_n);

%===============================================================

% ++++++++++++  Sensitivity Analysis ++++++++++++++

% using analytical equations
% sensitivity to entries in technology matrix
i_lambda = i_B_n*(inv(i_A_n)); % Incandescent uncertainty parameter


for k = 1:nB;
    for i = 1:mA
        for j = 1:nA
           if k == 1; % Entry 1 corresponds Dinitrogen dioxide
              i_dinitogen(i,j) = -(i_A_n(i,j))/i_g(k) * i_lambda(k,i)*i_s(j);
                           
            elseif k == 2; % Entry 2 corresponds Methane
              i_methane(i,j) = -(i_A_n(i,j))/i_g(k) * i_lambda(k,i)*i_s(j);
                
             elseif k == 3; % Entry 3 corresponds to CO2
                i_carbon(i,j) = -(i_A_n(i,j))/i_g(k) * i_lambda(k,i)*i_s(j);
             
           end
        end
     end
  end 


% sensitivity to entries in intervention matrix


for k = 1:nB;
    for i = 1:mA
        for j = 1:nA
            if k == 1; %Entry 1 corresponds Dinitrogen dioxide
              if i== k
                  i_dinitogen_B(i,j) = i_B_n(i,j)/i_g(k)*i_s(j);
              else
                  i_dinitogen_B(i,j) = 0;
              end
                           
            elseif k == 2; % Entry 2 corresponds Methane
              if i== k
                  i_methane_B(i,j) = i_B_n(i,j)/i_g(k)*i_s(j);
              else
                  i_methane_B(i,j) = 0;
              end
                
             elseif k == 3; % Entry 3 corresponds to CO2
              if i== k
                  i_CO2_B(i,j) = i_B_n(i,j)/i_g(k)*i_s(j);
              else
                  i_CO2_B(i,j) = 0;
             
            end
        end
    end
    end
end

% -delta% disturbance:
delta=0.01;

for i = 1:mA
   for j = 1:nA
                i_A_dlt=i_A_n;
                i_A_dlt(i,j)=i_A_n(i,j)*(1-delta);
                i_g_dlt=i_B_n*inv(i_A_dlt)*i_f;
        for k = 1:nB;
            if k == 1; %Entry 1 corresponds Dinitrogen dioxide
                i_dinitogen_m(i,j) = (i_g_dlt(k)-i_g(k))/i_g(k)/delta;
                           
            elseif k == 2; % Entry 2 corresponds Methane
              i_methane_m(i,j) = (i_g_dlt(k)-i_g(k))/i_g(k)/delta;
                
            elseif k == 3; % Entry 3 corresponds to CO2
                i_CO2_m(i,j) = (i_g_dlt(k)-i_g(k))/i_g(k)/delta;
           
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
      %names of variable has two parts: name of unit process + entry in
      %linear space
     name1=names_p(np);
     if nr==0
         name2=names_ls(mA);
     else
         name2=names_ls(nr);
     end;
     
     names(i)=strcat('  Amount of',{' '},name2, ' in',{' '},name1);
   
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



% triangular distribution
% For triangular distribution,we need three parameters, the range
% [min,max], and the mode (most probable value).

min=(0.442*0.3)/0.42;
max=(0.442*0.3)/0.25;
md=(0.442*0.3)/0.33;

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

% for each sample, we need to calculate a new CO2 emission
% the entries we want to change are for Electricity, bituminous coal, at power plant, which corresponds to column 33 

i_A_mc_t=i_A_n;
i_B_mc_t=i_B_n;
i_A_mc_t(41,33)=(33/30)*i_A_mc_t(41,33) % changing the value of electricity output to account for efficency change.
for i=1:N_s
    i_A_mc_t(38,33)=-R_t(i);
   
    i_B_mc_t(1,33)=-(i_B_n(1,33)*(R_t(i)/0.442));
    i_B_mc_t(2,33)=-(i_B_n(2,33)*(R_t(i)/0.442));
    i_B_mc_t(3,33)=-(i_B_n(3,33)*(R_t(i)/0.442));
   
    i_g_mc_t=i_B_mc_t*inv(i_A_mc_t)*i_f;
    carbon_t(i)=i_g_mc_t(3);
end;
figure(3);
hist(carbon_t,M_b);

% now we check the combined effects on carbon footprint. Assume the two
% random variables are independent. Note this may not always the case.
% Joint distributions are needed but usaully these are difficult to find.
% Also, we may have many random varibales invovled.


% determine the x-th percentile of the data. 
% usually 25% and 75% or 5% and 95% are used.

y5=prctile(carbon_t, 5)
y95=prctile(carbon_t, 95) % we used 5th and 95th percentile according to the question

% determine the confidence interval
pd = fitdist(carbon_t','Normal');
ci = paramci(pd,'Alpha',.05)
