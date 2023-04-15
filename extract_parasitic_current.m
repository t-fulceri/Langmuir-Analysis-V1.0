function [ fitted_parasitic ] = extract_parasitic_current( meas_directory, filename, current_sign )
%extract_parasitic_current Summary of this function goes here
%   Detailed explanation goes here

full_file_path = strcat(meas_directory,filename);

%Check metadata
metadata = h5info(full_file_path);
%Extract timestep from metadata [s]
time_step = metadata.Groups.Attributes(2,1).Value(1);
%Extract date and time from metadata
date = metadata.Datasets.Attributes(1,1).Value
time = metadata.Datasets.Attributes(2,1).Value

%Read data
data = h5read(full_file_path,'/dataset');
%Transpose
data = data';
%Extract current
current = data(:,2);
%Optional inversion of current signal (depends on the data source)
current = current_sign.*current;

[xData, yData] = prepareCurveData( [], current );

% Set up fittype and options.
ft = fittype( 'fourier1' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0 0 0 0.031403365189822];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
cv = coeffvalues(fitresult);
a0 = cv(1);
a1 = cv(2);
b1 = cv(3);
w  = cv(4);

t_step_vec = 1:1:length(current);
fitted_parasitic = (a0 + a1*cos(t_step_vec*w) + b1*sin(t_step_vec*w))';

end

