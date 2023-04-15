function [ return_status, output_args ] = extract_plasma_param_IV( voltage_cycle, current_cycle, probe_tip_area, ion_mass_number)

    %USEFUL PHYSICAL CONSTANTS----------------------------------
    %Electron_charge [C]
    e = -1.6021766208 * 1e-19;
    %Boltzmann constant [J/K]
    k_B = 1.38064852 * 1e-23;
    %Atomic Mass Unit [kg]
    AMU = 1.660539040*1e-27;
    %-----------------------------------------------------------
    
    %INTERNAL PARAMETERS----------------------------------------
    %Polynomial degree for fitting
    n = 20;
    %Margin percent from borders of V axis where polynomial approximation
    %is considered valid
    margin_percent = 0.2;
    %Minimum allowed temperature [eV]
    T_MIN = 0;
    %Maximum allowed temperature [eV]
    T_MAX = 20;
    %-----------------------------------------------------------
    
    [ reconstruct_curve_out ] = reconstruct_curve( voltage_cycle, current_cycle, false );

    V_axis = reconstruct_curve_out.x;
    I_of_V = reconstruct_curve_out.y_of_x;
    delta_I = reconstruct_curve_out.delta_y;

    output_args.V_axis = V_axis;
    output_args.I_of_V = I_of_V;
    output_args.delta_I = delta_I;

    %TO DO BEFORE EVERYTHING ELSE:
    %CHECK IF THERE IS PLASMA AT ALL!!!
    if max(I_of_V(I_of_V > 0))/abs(min(I_of_V)) <= 1
        fprintf('No plasma present\n');
        return_status = 1;
        return
    end

    %Extract V_float (voltage at which the current crosses zero)
    [ d, ind_V_float ] = min(abs(I_of_V));
    if ismember(ind_V_float,1:length(V_axis))
        V_float = V_axis(ind_V_float);
    else
        fprintf('Floating potential is out of range\n');
        V_float = nan;
    end
    output_args.V_float = V_float;

    %Extract a first estimate of ion saturation current
    I_sat_i_1st = mean(I_of_V(I_of_V < 0));
    output_args.I_sat_i_1st = I_sat_i_1st;
    
    %Approximate the dataset with a
    %Piecewise Cubic Hermite Interpolating Polynomial (PCHIP)
    I_of_V_pchip = pchip(V_axis(~isnan(I_of_V)),I_of_V(~isnan(I_of_V)),V_axis);
    fprintf('IV_curve fitted with PCHIP\n');

    %Fit the PCHIP with a normal polynomial
    warning('off','all');
    p_IV = polyfit(V_axis,I_of_V_pchip,n);
    fprintf('PCHIP fitted with polynomial\n');
    warning('on','all');

    %Keep note of the polynomial-approximate IV-curve
    I_of_V_poly = polyval(p_IV,V_axis);
    output_args.I_of_V_poly = I_of_V_poly;
    fprintf('I(V) values extracted from polynomial model\n');

%     f = figure
%     movegui(f,'center');
%     plot(V_axis,I_of_V,'*');
%     hold on
%     plot(V_axis,I_of_V_pchip,'-');
%     plot(V_axis,I_of_V_poly,'-');
%     xlabel('Voltage');
%     ylabel('Current');
%     title('I(V) curve, PCHIP approximation (red), and polynomial approximation (yellow)');

    %Polynomial approximation is not reliable near the extremes of the
    %V_axis interval considered:
    L = length(V_axis);
    
    delta_restricted = floor(margin_percent*L);
    ind_restricted = (1+delta_restricted):1:(L-delta_restricted);
    output_args.ind_restricted = ind_restricted;

    %Normalize it
    N_V = max(V_axis) - min(V_axis);
    N_I = max(I_of_V) - min(I_of_V);
    output_args.N_V = N_V;
    output_args.N_I = N_I;

    
    
    %Shift upwards if the minimum is negative or zero
    minI = min(I_of_V_poly);
    shift = 0;
    if minI < 0
        shift = - minI;
        
    else
        if minI == 0
            shift = delta_I;
        end
    end
    output_args.shift = shift;
    
    %Keep note of log(I + shift)
    log_I_of_V_poly = log(I_of_V_poly + shift);
    output_args.log_I_of_V_poly = log_I_of_V_poly;
    
    %Take first derivative of log(I)
    D_log_I_of_V_poly = diff(log_I_of_V_poly)./diff(V_axis);
    output_args.D_log_I_of_V_poly = D_log_I_of_V_poly;

    %Take the second derivative of log(I)
    d2_log_I = diff(D_log_I_of_V_poly);
    dV2 = diff(V_axis).^2;
    D2_log_I_of_V_poly = d2_log_I./dV2(1:end-1);
    output_args.D2_log_I_of_V_poly = D2_log_I_of_V_poly;
    
    %Look for interesting regions
    epsilon = 1e-3;
    exit = false;
    DELTA = 10;
    MAX_ITERATIONS = 1000;
    n_iter = 1;
    while(~exit)
        %Look only in region where polynomial approximation is valid
        ind0 = ind_restricted;
        %Only positive slope
        ind1 = intersect(ind0,find( D_log_I_of_V_poly > 0));
        %Almost vanishing second derivative
        ind2 = intersect(ind1, find( abs(D2_log_I_of_V_poly) < epsilon));
        
        %For each local maximum, construct a range of indices
        [pks, locs] = findpeaks(D_log_I_of_V_poly);
        peak_ranges = NaN(2*DELTA+1,length(locs));
        for nr = 1:length(locs)
            %Select only those in the ind_restricted range
            if ~isempty(find(ind_restricted == locs(nr)))
                peak_ranges(:,nr) = locs(nr)-DELTA:locs(nr)+DELTA;
            end
        end
        %Discard the NaN columns
        peak_ranges = peak_ranges(:,find(~isnan(peak_ranges(1,:))));
        
        %Intersect each index range with the previous conditions and
        %produce a list of candidates for the fit range
        N_cand = size(peak_ranges,2);
        ind_exp_candidates = cell(1,N_cand);
        for nie = 1:N_cand
            ind_exp_candidates{nie} = intersect(ind2,peak_ranges(:,nie));
        end
        
        %Take note of the length of the various candidates
        L_array = 1:N_cand;
        for i = 1:N_cand
            L_array(i) = length(ind_exp_candidates{i});
        end
        
        %Are there enough points for the fit?
        [m,m_l] = min(L_array);
        %if sum(~isnan(I_of_V(ind_exp_candidates{m_l}))) >= 6
        if min(L_array) >= 4
            exit = true;
        else
            %fprintf('Increasing epsilon...\n');
            epsilon = epsilon + 1e-3;
            n_iter = n_iter + 1;
            if n_iter >= MAX_ITERATIONS
                exit = true;
            end
        end

    end
    
    %Here: cycle thorugh the different fit range candidates and produce a
    %number of different [T,C] readings
    
    %Prepare output temperature array:
    T_e_lin_fit_array = NaN(N_cand,1);
    %Prepare output multiplicative constant array:
    C_lin_fit_array = NaN(N_cand,1);
    %Prepare fitresult and gof array
    fitresult_array = cell(N_cand,1);
    gof_array = cell(N_cand,1);
    
    for i_fit = 1:N_cand
        
        %Try to extract temperature from linear fit of log(I)
        V_tf = V_axis(ind_exp_candidates{i_fit});
        log_I_tf = log_I_of_V_poly(ind_exp_candidates{i_fit});

        N_log_I = 1;

        V_tf_norm = V_tf/N_V;
        log_I_tf_norm = log_I_tf/N_log_I;

        [xData, yData] = prepareCurveData( V_tf_norm, log_I_tf_norm );
        % Set up fittype and options.
        ft = fittype( 'poly1' );
        % Fit model to data.
        try
            [fitresult, gof] = fit( xData, yData, ft );
            fitresult_array{i_fit} = fitresult;
            gof_array{i_fit} = gof;
            cv = coeffvalues(fitresult);
            %Temperature = inverse of slope of log(I)
            T_e_lin_fit_array(i_fit) = 1/(cv(1)/N_V);
            %Additive constant
            C_lin_fit_array(i_fit) = cv(2)*N_log_I;
        catch
            fitresult_array{i_fit} = [];
            gof_array{i_fit} = [];
            cv = nan;
            %Temperature = inverse of slope of log(I)
            T_e_lin_fit_array(i_fit) = nan;
            %Additive constant
            C_lin_fit_array(i_fit) = nan;
        end

        
    end
    %Select only the realistic readings
    n_good_cand = find(and(T_e_lin_fit_array >= T_MIN,T_e_lin_fit_array <= T_MAX));
    
    T_e_lin_fit = T_e_lin_fit_array(n_good_cand);
    C_lin_fit = C_lin_fit_array(n_good_cand);
    ind_exp_cell = ind_exp_candidates(n_good_cand);
    fitresult = fitresult_array(n_good_cand);
    gof = gof_array(n_good_cand);

    %Convert electron temperature(s) to SI units (Kelvin)
    T_e_SI = T_e_lin_fit.* ((-e)/k_B);
    %Convert ion mass from AMU to kg
    m_i = ion_mass_number*AMU;
    %Calculate electron charge density/ies [C*(m^-3)]
    rho_e = (I_sat_i_1st/probe_tip_area) .* sqrt(2*pi*m_i./(k_B*T_e_SI));
    %Calculate electron density [m^-3]
    n_e_derived = rho_e/e;
    %Remove negative values
    n_e_derived = n_e_derived(n_e_derived >= 0);
    
    output_args.T_e_lin_fit = T_e_lin_fit;
    output_args.C_lin_fit = C_lin_fit;
    output_args.ind_exp_cell = ind_exp_cell;
    output_args.fitresult = fitresult;
    output_args.gof = gof;
    
    output_args.n_e_derived = n_e_derived;
    
    return_status = 0;

end

