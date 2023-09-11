function SEIR_Omicron_Fitting
  %This function fits the first set of BSfludata to an SIR model
  clear all
  close all
  clc
  
  format long

  %Define array with t-coordinates of the data
  tdata = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,  ...
      25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47, ...
      48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70, ...
      71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93, ...
      94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,...
      113,114,115,116,117,118,119,120,121,122,123,124,125,126]; 

  %Define array with y-coordinates of the data
  qdata = 36*[4579,4602,4641,4642,4714,4806,4902,5073,5205,5387,5637,5896,6164, ...
      6497,6923,7431,7975,8609,9248,9919,10603,11301,11905,12523,13088,13599, ...
      14094,14558,14873,15180,15476,15668,15816,15933,15990,15989,15914, ...
      15711,15502,15330,15021,14750,14540,14329,14144,13900,13653,13443, ...
      13206,12959,12760,12525,12355,12073,11810,11502,11161,10739,10308, ...
      9875,9447,9032,8665,8262,7932,7579,7244,6981,6722,6457,6218,5978, ...
      5773,5580,5340,5156,4968,4822,4622,4427,4229,4088,3921,3747,3564, ...
      3423,3310,3197,3091,2984,2906,2812,2735,2664,2601,2538,2465,2385, ...
      2319,2212,2117,2027,1944,1893,1818,1761,1755,1736,1745,1748,1751, ...
      1771,1777,1769,1761,1731,1723,1709,1738,1777,1827,1898,1974,2029, ...
      2094,2147]; 
  
  n = length(qdata);
  %t-mesh for the solution of the differential equation
  tforward = [1:0.01:126];

  %Selects the points in the solution corresponding to the t values of tdata
  tmeasure = [1:100:12501]';
 
  %Initial values of parameters to be fitted
  a = 0.01; %time sick
  b = 0.000009; %transmissibility
  c = 0.01; %time exposed
  suscep_init = 5000000.0;
  expos_init = 10000.0;
  infec_init = qdata(1);

  function dy = model_1(t,y,k) %DE
    a = k(1); % Assignes the parameters in the DE the current value of the parameters
    b = k(2);
    c = k(3);
    dy = zeros(3,1); %Initializes dy array
    dy(1) = - b * y(1) * y(3); %Susceptible equation
    dy(2) = b * y(1) * y(3) - c * y(2); %Exposed equation
    dy(3) = c * y(2) - a * y(3); %Infected equation
  end

  function error_in_data = moder(k) % computing the error in the data
    [T, Y] = ode23s(@(t,y)(model_1(t,y,k)),tforward,[suscep_init, expos_init, infec_init]);
    % solves the DE; output is written in T and Y

    q=Y(tmeasure(:),3); %Assigns the y-coordinates of the solution at the 
                        %t-coordinates of tdata
    
    error_in_data = sum((q' - qdata).^2); %Computes SSE
    AIC = n*log(error_in_data/n)+2*(length(k0)+1) + (2*(length(k0)+1)*...
     (length(k0)+2))/(n-(length(k0)+1)-1) % computes AIC
  end

  k0 = [a, b, c]; % main routine; assigns initial values of parameters
  [T, Y] = ode23s(@(t,y)(model_1(t,y,k0)),tforward,[suscep_init, expos_init, infec_init]);
  %Solves the DE with the initial values of the parameters
  yint = Y(tmeasure(:),3);
  %Assigns the y-coordinates of the solution at tdata to yint
  
  %Generate plot with data and regression based only on initial parameter
  %values. 
  figure(1)
  subplot(1,2,1);
  plot(tdata,qdata,'r.');
  hold on
  plot(tdata,yint,'b-'); 
  xlabel('time in days');
  ylabel('Number of cases');
  axis([0, 130, 0, 36*20000]);

[k,fval] = fminsearch(@moder,k0); %Minimization routine; assigns the new
                                  %values of parameters k and the SSE to fval
disp('optimal values ');
disp(k);

[T, Y] = ode23s(@(t,y)(model_1(t,y,k)),tforward,[suscep_init, expos_init, infec_init]);
    %Solving the DE with the final values of the parameters

yint = Y(tmeasure(:),3); %Computing the y-coordinates corresponding to 
                         %the tdata

%Plots data along with final regression for the data.
subplot(1,2,2)
plot(tdata,qdata,'r.');
hold on
  plot(tdata,yint,'b-');
  xlabel('time in days');
  ylabel('Number of cases');
  axis([0, 130, 0, 36*20000]);
end



