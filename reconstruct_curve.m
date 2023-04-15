function [ output_args ] = reconstruct_curve( x_of_t, y_of_t , interpolate_nan)
%RECONSTRUCT_CURVE Reconstructs the y = y(x) curve starting from the x(t) and
%y(t) signals, taking x as the independent variable
%   Warning: x and y must be of the same length and must represent
%   timeseries of measurements taken simultaneously

%Basically follow the scheme
%x(t) -> t -> (y_1(t),...,y_n(t)) -> y(t) -> y(x)

if nargin == 2
    interpolate_nan = true;
end

%Extract the x axis (independent variable)
x = signal_to_axis(x_of_t);

%Prepare the dependent variable
L = length(x);
y_of_x = zeros(1,L);
%Calculate resolution of y
dy = diff(y_of_t);
delta_y= ones(1,L)*min(abs(dy(dy~=0)))/2;

%Scan the x values range
for k = 1:L
    %Look for the time indices where the x value is present
    x_value = x(k);
    ind = find(x_of_t == x_value);
    %Look for the corresponding y value(s)
    y_spread = y_of_t(ind);
    %Average the the y values to get a single value
    if(~isempty(y_spread))
        y_of_x(k) = mean(y_spread);
    else
        y_of_x(k) = nan;
    end
end



% %Keep only the array elements corresponding to non-Nan values
% ind_good = find(~isnan(y_of_x));
% x = x(ind_good);
% y_of_x = y_of_x(ind_good);

%Order by increasing value of the x elements
[x, x_order] = sort(x);
y_of_x = y_of_x(x_order);


if interpolate_nan
    %Interpolate y_of_x where it is Nan

    %Take note of the array elements which are not Nans
    number_indices = find(~isnan(y_of_x));

    %Count how many indices correspond to numbered values
    L_N = length(number_indices);
    %Prepare two arrays of indices:
    %One for the Left Boundary Indices (LBI)
    LBI = nan(L_N,1);
    %One for the Right BOundary Indices (RBI)
    RBI = nan(L_N,1);

    %Between a LBI position and a RBI position there is at least one Nan value,
    %which must be interpolated

    %Fill the LBI and RBI according to therir description
    for k = 1:L_N
        %If the next element to the right is a nan...
        if k < L_N && isnan(y_of_x(number_indices(k)+1))
            %...than it is a Left Boundary Index
            LBI(k) = number_indices(k);
        end
        %If the previous element to the left is a nan...
        if k > 1 && isnan(y_of_x(number_indices(k)-1))
            %...than it is a Right Boundary Index
            RBI(k) = number_indices(k);
        end
    end

    %Keep only the necessary elements of the LBI and RBI arrays
    LBI = LBI(~isnan(LBI));
    RBI = RBI(~isnan(RBI));

    %Take note of the common length
    common_length = min([length(LBI),length(RBI)]);

    %Bridge all the Nan gaps with a linear interpolation
    for k = 1:common_length
        closed_interval = linspace(y_of_x(LBI(k)),y_of_x(RBI(k)),RBI(k)-LBI(k)+1);
        open_interval = closed_interval(2:end-1);
        y_of_x(LBI(k)+1:RBI(k)-1) = open_interval;
    end
end

%Copy the curve and the uncertainty into the output structure
output_args.x = x;
output_args.y_of_x = y_of_x;
output_args.delta_y = delta_y;

end
