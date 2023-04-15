clear all
close all

full_file_path_00 = 'C:\Users\tzf\Desktop\VINETA II Lab\MATLAB\Development\IV Timeseries Reader\IV_readings_00.mat';
load(full_file_path_00);
delta_V = IV_readings.delta_voltage;
delta_I = IV_readings.delta_current;
IV_readings_00 = IV_readings;

full_file_path_01 = 'C:\Users\tzf\Desktop\VINETA II Lab\MATLAB\Development\IV Timeseries Reader\IV_readings_01.mat';
load(full_file_path_01);
IV_readings_01 = IV_readings;

full_file_path_02 = 'C:\Users\tzf\Desktop\VINETA II Lab\MATLAB\Development\IV Timeseries Reader\IV_readings_02.mat';
load(full_file_path_02);
IV_readings_02 = IV_readings;

full_file_path_03 = 'C:\Users\tzf\Desktop\VINETA II Lab\MATLAB\Development\IV Timeseries Reader\IV_readings_03.mat';
load(full_file_path_03);
IV_readings_03 = IV_readings;

full_file_path_04 = 'C:\Users\tzf\Desktop\VINETA II Lab\MATLAB\Development\IV Timeseries Reader\IV_readings_04.mat';
load(full_file_path_04);
IV_readings_04 = IV_readings;

clear IV_readings;

nc = 8;

V_0 = IV_readings_00.voltage_cycles_rise(:,nc);
I_0 = IV_readings_00.current_cycles_rise(:,nc);

V_1 = IV_readings_01.voltage_cycles_rise(:,nc);
I_1 = IV_readings_01.current_cycles_rise(:,nc);

V_2 = IV_readings_02.voltage_cycles_rise(:,nc);
I_2 = IV_readings_02.current_cycles_rise(:,nc);

V_3 = IV_readings_03.voltage_cycles_rise(:,nc);
I_3 = IV_readings_03.current_cycles_rise(:,nc);

V_4 = IV_readings_04.voltage_cycles_rise(:,nc);
I_4 = IV_readings_04.current_cycles_rise(:,nc);


V = [V_0;V_1;V_2;V_3;V_4];
I = [I_0;I_1;I_2;I_3;I_4];

[V,sort_ind] = sort(V);
I = I(sort_ind);




[out] = reconstruct_curve_improved(V,I);
V_axis = out.x;
I_of_V = out.y_of_x;

L = length(V_axis/2);
n_sgolay = 5;
frame_length = round(L/2);
if mod(frame_length,2) == 0
    frame_length = frame_length + 1;
end

%Apply Savitzky–Golay filter
I_of_V_sgolay = sgolayfilt(I_of_V,n_sgolay,frame_length);
fprintf('Savitzky–Golay filter applied...\n');

%Focus on a meaningful range of indices (exclude margins)
margin_percent = 0.1;
margin = round(margin_percent*length(V_axis));
margin_range = margin:length(V_axis)-margin;

figure
set(gcf,'Renderer','painters');
ax_IV = axes;
movegui('north');
vec = ones(length(V_axis),1);
errorbar(ax_IV,V_axis,I_of_V,delta_I.*vec,'ko');
hold(ax_IV,'on');
plot(ax_IV,V_axis,I_of_V_sgolay,'b-','LineWidth',2);
title(ax_IV,['IV-characteristic at cycle number ',num2str(nc)]);
xlabel(ax_IV,'Voltage [V]');
ylabel(ax_IV,'Current [A]');
xlim_temp = xlim;
hold(ax_IV,'off');

figure
set(gcf,'Renderer','painters');
movegui('south');
ax_logI = subplot(2,1,1);
if min(I_of_V_sgolay) < 0
    shift = -(1.1)*min(I_of_V_sgolay);
else
    shift = 0;
end
logI = log(I_of_V_sgolay + shift);
plot(ax_logI,V_axis,logI);
xlim(xlim_temp);
hold(ax_logI,'on');

n_poly = 10

p_logI = polyfit(V_axis,logI,n_poly);
logI_poly = polyval(p_logI,V_axis);
plot(ax_logI,V_axis,logI_poly,'-r');
hold(ax_logI,'off');

ax_DlogI = subplot(2,1,2);
DlogI = diff(logI_poly)./diff(V_axis);
plot(ax_DlogI,V_axis(1:end-1),DlogI);

%Find the peaks od DlogI in the meaningful range (exclude the margins), take the
%maximum peak as a marker for a fitting region
[pks,locs] = findpeaks(DlogI);
locs = locs(locs > margin_range(1) & locs < margin_range(end))
pks = DlogI(locs);
[mp,mp_ind_in_locs] = max(pks);
%First approximation of temperature, NOT fundamental: sometimes it is
%completely out of this world (but we keep it)
T_1st_approx = 1./mp
xlim(xlim_temp);


%The location of the greatest peak is an indication of where to apply the
%fitting procedure
delta_left  = 1000;
delta_right = 10;
%fit_range = intersect(find(DlogI > mp*0.8),locs(mp_ind_in_locs)-delta_left:locs(mp_ind_in_locs)+delta_right);
data_range = locs(mp_ind_in_locs)-delta_left:locs(mp_ind_in_locs)+delta_right;
fit_range = intersect(margin_range,data_range);
if isempty(fit_range)
    fprintf('Fit range does not contain any indices \n');
    return
end

V_tf = V_axis(fit_range);
I_tf = I_of_V(fit_range);
%WARNING: PROVISIONAL
I_tf = I_of_V_sgolay(fit_range);


[xData, yData] = prepareCurveData( V_tf, I_tf );

ft = fittype( 'I_0 + a*V + b*exp(c*V)', 'independent', 'V', 'dependent', 'I' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf 0 0.001 0.01];
opts.MaxFunEvals = 1000;
opts.MaxIter = 1000;
opts.Robust = 'Bisquare';
opts.StartPoint = [0 0.0001 0.1 mp];
opts.TolFun = 0.0001;
opts.TolX = 0.0001;
opts.Upper = [0 Inf 1 Inf];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
cv = coeffvalues(fitresult);
hold(ax_IV,'on');
plot_range = (fit_range);

I_0 = cv(1)
a = cv(2)
b = cv(3)
c = cv(4)

I_of_V_fitted = I_0 + a.*V_axis + b.*exp(c.*V_axis);
plot(ax_IV,V_axis(plot_range),I_of_V_fitted(plot_range),'r-','LineWidth',2);

I_i = I_0 + a.*V_axis;
plot(ax_IV,V_axis,I_i,'r--','LineWidth',2);

T_fit = 1/c


%Extract floating potential and draw a vertical line at its position
[ d, ind_V_float ] = min(abs(I_of_V_sgolay));
if ismember(ind_V_float,1:length(V_axis))
    V_float = V_axis(ind_V_float)
else
    fprintf('Floating potential is out of range\n');
    V_float = nan;
end
axes(ax_IV);
x_bound = [V_float V_float];
y_bound = get(ax_IV,'YLim');
line(x_bound,y_bound,'Color',[0 1 0],'LineWidth',2)

%Extract plasma potential and to draw a vertical line at its position
V_plasma = V_float + T_fit
axes(ax_IV);
x_bound = [V_plasma V_plasma];
y_bound = get(ax_IV,'YLim');
line(x_bound,y_bound,'Color',[1 0 1],'LineWidth',2)

h_leg = legend(ax_IV,'Data points','Savitzky-Golay filtered data','Fit with model I(V) = I_0 + aV + be^{cV}','Ion current I_i(V) = I_0 + aV','Floating potential','Plasma potential');
set(h_leg,'FontSize',10);

I_sat_i_fit = I_0 + a*V_float

