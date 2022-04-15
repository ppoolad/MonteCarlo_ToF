%Sim deadtime effect on a lambda
function final_pdf = analytical_model(lambda_t, system,transient_active)
    n_spad = system.n_spad;
    n_pix = system.n_pix;
    n_tdc = system.n_tdc;
    sim_res = system.sim_res;
    arrival = system.arrival;
    tdc_dead = system.tdc_dead;
    spad_dead = system.spad_dead;
    adjust_factor = system.adjust_cd;

    %% Photon Arrival times have exponential distribution:
    if (arrival == "First")
        lambdas_sum = cumsum(lambda_t).*sim_res;
        r_t = lambda_t .* exp(-1*lambdas_sum);        
    else
        r_t = lambda_t;% .* sim_res; 
    end

    %% Apply deadtime effect for SPADS

    %transient: dead time
    if(transient_active)
        r_t_dd = transient(r_t,spad_dead,sim_res)/sim_res;
    else
        r_t_dd =  r_t ./ ( 1 + r_t .* spad_dead);
    end
    %% Apply CD

    c_f = system.c_f;
    c_t = system.c_t;
    if(c_f > 1)
        
        f_c_t = coincidence_detection(r_t_dd, c_f, c_t, n_spad, sim_res); %analytical_lidar(r_t_dd, c_f, c_t, n_spad, sim_res,adjust_factor);    
        %an older less accurate model 
        %f_c_t = analytical_lidar(r_t_dd, c_f, c_t, n_spad, sim_res,1); <-  
    else
        f_c_t = n_spad .* r_t_dd;
    end
    %% Apply deadtime effect for TDC
    if(tdc_dead > sim_res)
        h = tdc_dd_pn(f_c_t, tdc_dead, n_tdc, n_pix);
        r_eff = f_c_t .* (1-h);
    else
        r_eff = f_c_t;
    end
    final_pdf = r_eff;
    %% Voila
    
end