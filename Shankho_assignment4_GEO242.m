%% MATLAB script to model the Slip Rate and Locking Depth of a fault from GPS data
% Author: Shankho Niyogi, 27th Feb, 6.30 pm.

clear

d = 2:4:20;% Matrix containing range of Locking Depth values 

x1 = -100:100;%index value ranges for modeling

%defining zero matrices to store values beforehand to prevent run time allocation
data = zeros(3,11);
vel = zeros(length(d),length(data(1,:)));
error_est = zeros(length(d),1);

%Assigning values manually from the pdf file
data(1,:) = [-22 -67 4 15 2 -55 6 20 -3 -17 1];%Distance from fault in km
data(2,:) = [37.5 42.9 25.8 23.4 29.3 42.5 25.4 19.8 32.0 36.0 29.7];%Slip rate observed in mm/yr
data(3,:) = [0.6 1.4 0.5 0.3 0.4 0.5 0.8 0.5 1.2 0.6 0.8];%Slip rate uncertainty in mm/yr
% data(4,:) = [3.0 4.1 2.9 1.9 2.5 4.4 3.9 3.4 0.9 1.6 3.5];
% data(5,:) = [0.4 0.5 0.5 0.1 0.3 0.5 0.7 0.5 0.3 0.4 0.4];
%% Matrix Inversion Section using the form A*x = b
% The matrix A contains the equation (1/pi)*atan(x/D) as first column and
% 1s as second column
% Note the matrix x contains the slip rate value and constant which we need to find out for different locking Depth values.
% Matrix b contains the velocity values.

for j=1:length(d)% Loop to run through different Locking Depth values
    
    data_vel = -(1/pi)*atan(data(1,:)/d(j));%creating the first column of A matrix
    
    A = [data_vel ; ones(1,length(data_vel))]';%Creating A matrix.
    R = qr(A);%storing the result of QR factorization.
    b = data(2,:)';%creating the b matrix
    
    x = R\(R'\(A'*b));%Finding the 2 values of x, slip rate and constant c
    residual =  b - (A*x);% calculating the residual
    error = R\(R'\(A'*residual));% calculating the error
    x = x + error;%updating the values of x to reduce the error between modelled and actual data
    
    %% Creating datasets for plotting
    
    data_vel2 = -(x(1,1)/pi)*atan(data(1,:)/d(j)) + x(2,1);% Calculating new velocities using the new slip rate and constant value
    
    error_est(j) = sum(sqrt(mean(((data(2,:) - data_vel2)/data(3,:))^2)));% calculating the total error estimation between the actual and modelled result
    
    for i=1:length(x1)
        vel(j,i) = -(x(1,1)/pi)*atan(x1(i)/d(j)) + x(2,1);% generating a set of values to be plotted as a line along with modelled results
    end
end

%% Data plotting section.

figure(1)
errorbar(data(1,:),data(2,:),data(3,:),'LineStyle','none');grid on;grid minor;legend('GPS data observations');%Displaying actual datavalues with uncertainty
hold on
for j = 1:length(d)%Loop to plot different modelled results  
    txt = sprintf('Modeled Slip rate for %d km Locking Depth, RMSE error: %f mm/yr',d(j),error_est(j));%Legend entry for each modelled data
    plot(x1,vel(j,:),'DisplayName',txt);% Plotting modelled values and custom legend text
end
xlim([-100 100])% Defining the x axis limits
xlabel('Distance from fault (km)');
ylabel('Modeled Slip Rate (mm/yr)');
title(sprintf('Modelled results for Slip rate of %f mm/yr for different Locking Depths',x(1,1)));
hold off
legend show
set(gca,'FontSize',12);