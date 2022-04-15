%% analytical only histogram %%
function hst = histModelScript_vect(distance, const, refl)
    %%%N = 500000;
     if ~exist('refl','var')
     % third parameter does not exist, so default it to something
      refl = 255*const.reflectivity;
    end
    %% lambda_sig = E_rx/(sqrt(2*pi)*sigma) * exp(-(t-t_tof)^2/(2*sigma^2)) * eta_pde/E_photon 
    %%%% CONSTANTS %%%%
    
    dist = distance;
    if(size(dist,2)>=1 && length(dist(:)) < 1000)
       dist =  reshape(dist',1,length(dist(:)));
    elseif(length(dist(:))>=1000)
        dist = tall(reshape(dist',1,length(dist(:)))); %%<<-------- think about this
    end
    reflectivity = single(refl)/255;
    %E_rx = reflectivity .* (const.E_30cm * (30e-2)^2 ./ dist.^2 * const.pulse_length);%%1e-12 * 5e-9; %% <5
    coverage_area = 4 .* dist.^2 .* tand(const.laser_div/2).^2; %laser illuminates this area
    lambertian_factor = const.Lens_D .^2 ./ (const.Lens_D .^2 + 4 .* dist.^2); %lens collects this much
    %maximum power at target 
    p_max_target = const.laser_max_p ./ coverage_area; %w/m2;
    p_max_reflected = reflectivity .* p_max_target; %reflected
    p_max_lens_collected = lambertian_factor .* p_max_reflected .* const.lens_opacity;
    pixel_projection = const.pix_area .* (dist ./ const.focal) .^ 2;
    p_max_pixel = p_max_lens_collected .* pixel_projection;
    %create the envelope
    
    
    %E_rx = reflectivity .* (const.laser_output_E ./4 ./ pi ./ dist.^2 ./ (tand(const.fov).^2)) .* const.pix_area ; % <- This should be distance in meters.
    sim_res = const.sim_res; %%1e-12;
    e_photon = const.e_photon;%%53.8e-19;
    etha_pde = const.etha_pde;  %%%0.5;
    background_rate = const.backgrount_rate;%%%33.4e6;
    c = const.c;%%physconst('LightSpeed');
    histogram_bins = const.histogram_bins;
    %histogram_res = histogram_bins * sim_res;%%%  50e-12;
    tdc_dead = const.tdc_dead;%% 5e-9;
    n_tdc = const.n_tdc;%%%1;
    n_pixels_to_tdc = const.n_pixels_to_tdc;%%% 3;
    fps = const.fps;%%% 200;
    rpt_freq = const.rpt_freq;%% 100e6;
    time_range = 1/rpt_freq;
    end_time = const.end_time;
    duration = end_time - const.start_time;
    N = rpt_freq/fps;
    coinc = const.coinc;
    spad_per_pix = const.spadPerPix;
    CT_Time = const.ct_time; % Not supported for now.
    
    
    %%%miss = false;
    %%%%%%%%%%%%%%%%%%%  
    if(isfield(const,'start_time'))
        start_time = const.start_time;
        %%disp('Delayed Start');
    else
        start_time = 0;
    end
    
    if(dist(1) < 1e-3)
        t_tof = dist - start_time;
    else
        t_tof = 2*dist/c - start_time;
    end
    miss = (false(size(distance)));
    miss(t_tof<0) = 1;

    laser_sigma = const.pulse_length / (2*sqrt(2*log(2)));
    sigma = sqrt((laser_sigma)^2 +(30e-12)^2 + (10e-12)^2);
    dt = 0:sim_res:0+time_range-sim_res;
    P_sig = p_max_pixel .* exp(-((dt' - t_tof).^2)/ (2 .* sigma.^2));
    %P_sig = E_rx .* normpdf(dt',t_tof,sigma); %%Power Distribution of Signal %
    lambda_sig = P_sig ./ e_photon .* etha_pde; %% Rate of signal

    %% R_eff = (lambda_bg + lambda_sig) * exp^(-1*integral(lambda_bg + lambda_sig,0,t))

    lambda_bg =  background_rate* etha_pde;%% Rate of Background
    lambdas = lambda_bg + lambda_sig; %Rate of recieved photons
    lambdas_sum = cumsum(lambdas).*sim_res; %Integrate
    
    %% %%%%%%%%%%% Supressed lambdaBG with coincidence detection %%%%%%%%%%
    % Considering there is only background:
    % we can calculate expected value of our lambda
    %lambda_sup = cd_supression(spad_per_pix,lambda_bg,lambda_sig,coinc,dt,CT_Time,sim_res,tdc_dead);
    %lambdas_cd = lambda_sup.bg + lambda_sup.sig;
    %lambdas_cd_sum = cumsum(lambdas_cd).*sim_res;
    %% Photon Arrival times have exponential distribution:
    if (const.arrival == "all")
        r_t = lambdas;
    else
        r_t = lambdas .* exp(-1*lambdas_sum);% .* sim_res; 
    end
    %r_ts = spad_per_pix.* lambdas .* exp(-1.*lambdas_sum);% .* sim_res; 
    %r_t_CD = lambdas_cd .* exp(-1.*lambdas_sum);
    %% Window Conv for Coincidence
    %% Experimental
%     nerl = spad_per_pix;
%     r_bg = nerl.*lambdas;
     CT_window = ones(uint32(CT_Time/sim_res),1);
%     k_erl = coinc - 1;
%     r_tn = nerl*lambdas .* exp(-nerl*lambdas_sum);% .* sim_res; 
%     Erlang_C = r_tn; %%nerl*lambda_bg.^k_erl .* dt.^(k_erl-1) .* exp(-nerl*lambda_bg.*dt)/factorial(k_erl-1);
%     erl_conv = conv2(CT_window,1,r_tn,'full').*sim_res;
%     rterl = r_bg(:,12) .* erl_conv(1:5000,12);
%     %%This is more like a PMF %%We need it padded
      r_t_ctf = conv2(CT_window,1,r_t,'full').* sim_res;
%       r_t_ctfa = conv2(CT_window,1,lambdas./(1+lambdas.*tdc_dead),'full').* sim_res;
      r_t_ct = r_t_ctf(1:length(r_t),:);
      r_t_CD = zeros(size(r_t_ct));
%      r_t_cd_n = zeros(size(r_t_ct));
%      %r_t_cd = r_t_ct .^ coinc;
%      %r_t_cd_null = 1-r_t_cd;
% %      
%       for i=coinc:coinc
       r_t_CD = spad_per_pix .* nchoosek(spad_per_pix-1,coinc-1) .*  r_t_ct.^(coinc-1) .* (1-r_t_ct).^(spad_per_pix - coinc-1);
%               %r_t_CD = r_t_CD + factorial(spad_per_pix)/factorial(spad_per_pix-i) .*  r_t_ct.^i .* (1-r_t_ct).^(spad_per_pix - i);
% %              r_t_cd_n = r_t_cd_n + nchoosek(spad_per_pix,coinc) .*  (1-(r_t_ct)).^i;
%       end
       %pad_r_t = padarray(r_t,0,0,'pre');
              
       r_t_CDx = (r_t.*sim_res).*((1-r_t.*sim_res).^(spad_per_pix-1)) .* r_t_CD;%%conv2(CT_window,1,r_t_CD,'full').*sim_res;
       %r_t_CDx = r_t_CD;
       r_t_CD = r_t_CDx(1:length(r_t),:);
%       r_t_CD = nchoosek(spad_per_pix-1,coinc-1) .* r_t_ct.^(coinc-1) .* (1-r_t_ct).^(spad_per_pix-1-coinc-1);
      %r_t_shift = padarray(r_t,length(CT_window));
%       r_t_CD = spad_per_pix.*(r_t.*r_t_CD);
      %figure; plot(rfa);
      %r_t_CD2 = nchoosek(spad_per_pix,coinc) .*  r_t.^coinc .* (1-r_t).^(spad_per_pix - coinc);
      %r_t_CDc = conv2(CT_window,1,r_t_CD2,'full').*sim_res;
      %r_t_CDcs = r_t_CDc(1:length(r_t),:);
      %%r_t_CD = r_t_CD .* r_t .* sim_res .* spad_per_pix .* (1 - r_t.*sim_res).*(spad_per_pix-1);
%     r_t_CD_All  = (coinc * nchoosek(spad_per_pix,coinc) .*  r_t_ctfa.^coinc .* (1-r_t_ctfa).^(spad_per_pix - coinc))./CT_Time;
% %    r_t_CD = nchoosek(spad_per_pix,coinc) .*  r_t_ct.^coinc;% .* (1-r_t_ct).^(spad_per_pix - coinc);
%     figure; plot(r_t_CD_All);
%     lambdas_sum_cd = cumsum(r_t_CD_All).*sim_res; %Integrate
%     r_t_cd_all_f = r_t_CD_All .* exp(-1*lambdas_sum_cd);
%     figure; plot(r_t_cd_all_f);

    %r_t_CD = r_t_CD/CT_Time;
    pn = tdc_dd_pn(r_t,tdc_dead,n_tdc,n_pixels_to_tdc);
    pn_cd = tdc_dd_pn(r_t_CD/sim_res,tdc_dead,n_tdc,n_pixels_to_tdc);
    %pn_cdc = tdc_dd_pn(r_t_CDcs,tdc_dead,n_tdc,n_pixels_to_tdc);
    %% Final
    %%figure; plot(pn_cd);
    rat = 1;%uint32(histogram_res/sim_res); %reshape factor
    histogram_res = sim_res;
    r_eff = zeros([size(r_t),2]);
    r_eff(:,:,1) = r_t .* (1-pn);
    r_eff(:,:,2) = 0.89 * r_t_CD .* (1-pn_cd);
    %reff_split = r_eff; %reshape(r_eff,size(r_eff,1)/rat,rat,size(r_eff,2),2); %split
    %reff_split_cd = reshape(r_eff_cd,rat,size(r_eff_cd,1)/rat,size(r_eff_cd,2));
    h = r_eff .* N; %reshape(sum(reff_split,1).*sim_res,size(reff_split,2),size(reff_split,3),2); %sum
    h(:,:,1) = (h(:,:,1) .* sim_res); %N times
    %h(:,:,2) = (h(:,:,2));
    if(length(duration) == 1)
       h_select = h(1:uint32(duration/histogram_res),:,:);
       h_select = squeeze(sum(reshape(h_select,size(h_select,1)/histogram_bins,histogram_bins,size(h,2),2)));
       tau = (0:duration/histogram_bins:duration-sim_res)' + zeros(histogram_bins, size(h,2));
    else
    h_select = zeros(histogram_bins,size(h,2),'single');
    tau = zeros(histogram_bins,size(h,2));
    for i = 1:size(h,2)
        temp = h(1:uint32(duration(i)/histogram_res),i);
        h_select(:,i) = sum(reshape(temp,size(temp,1)/histogram_bins,histogram_bins,size(temp,2)));
        
        tau(:,i) = 0:duration(i)/histogram_bins:duration(i)-histogram_res;
    end
    end
    %%hcd2 = sum(reshape(h_select,size(h_select,2)/histogram_bins,histogram_bins),1).^coinc;
    %plot(tau,h,'r','LineWidth',3)
     %ylim([0 1000])
     %xlim([0 10e-9])
    %%sum(h)
    
    hst.time = tau + start_time;
    %%%q= sum(reshape(h_select,size(h_select,1)/histogram_bins,histogram_bins,size(h_select,2)));
    hst.h = (h_select);%%reshape(q,size(q,2),size(q,3));
    %h_norm = normalize(hst.h,'norm',1);
    %%h_cd = h_norm.^coinc;
    %hst.h_norm = h_norm;
    %% new HCD %%
%     h_cd = zeros(size(h_norm));
%     for i=coinc:spad_per_pix
%         h_cd = h_cd + nchoosek(spad_per_pix,coinc) *  h_norm.^i .* (1-h_norm).^i;
%     end
%     
%     ph_ratio = sum(h_cd)*spad_per_pix;
%     hst.hst_cd = h_cd .* sum(hst.h).*ph_ratio;
     hst.miss = miss;
%     
%     if(const.plot)
%         figure;
%         bar(tau+start_time,hst.h);
%         title("Normal Histogram");
%         xlabel("Time (ns)");
%         ylabel("Photon Count");
%         figure;
%         bar(tau+start_time,hst.hst_cd);
%         title("Histogram with Coincidence Detection of n=2");
%         xlabel("Time (ns)");
%         ylabel("Photon Count");
%     end
%     if(const.save)
%         f = figure('visible','off');
%         bar(tau+start_time,hst.h);
%         title("Normal Histogram");
%         xlabel("Time (ns)");
%         ylabel("Photon Count");
%         figname = sprintf(const.save_loc + "bc_hist_%.1f_S%.2f_E%.2f.png",distance,start_time*1e9,end_time*1e9);
%         saveas(f,figname);
%         close(f);
%         f=figure('visible','off');
%         bar(tau+start_time,hst.hst_cd);
%         title("Histogram with Coincidence Detection of n=2");
%         xlabel("Time (ns)");
%         ylabel("Photon Count");
%         figname = sprintf(const.save_loc + "cd2_hist_%.1f_S%.2f_E%.2f.png",distance,start_time*1e9,end_time*1e9);
%         saveas(f,figname);   
%         close(f);        
%     end
    
end
