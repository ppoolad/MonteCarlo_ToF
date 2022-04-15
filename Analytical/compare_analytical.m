close all
%coincidece time and factors
c_t = 3e-9;
c_f = 16;
%how many tdc bins?
TDC_bins = 512;
%number of iterations 
N = 20;
%transient??
transient_active = 1;

%how many spad in a sipm?
n_spad_per_pix = 133;

%large numbers throw warning for nchoosek so supress them
warning('off','MATLAB:nchoosek:LargeCoefficient')
sim_res = 10e-12;
%where to read lambda vector and montecarlo histogram matrix generated from
%python
sim_path = "./outputs/milestone/sensl/100k_20it/cd16/";
load(sim_path+"lambdas_strt100_stp500.mat");
load(sim_path+"N20_strt100_stp500.mat");
%list of variables:
lambda_name_list = ["pixel0_alambdas100", "pixel1_alambdas200", "pixel2_alambdas300", "pixel3_alambdas399", "pixel4_alambdas500"];
lvar_name_list = ["pixel0_all_ar_tof100","pixel1_all_ar_tof200","pixel2_all_ar_tof300","pixel3_all_ar_tof399","pixel4_all_ar_tof500"];

%true tofs for SNR calc
true_tof_list  = uint32([100e-9, 200e-9, 300e-9, 400e-9, 500e-9]/sim_res);

SNR_A = zeros(size(true_tof_list));
SNR_MC = zeros(size(true_tof_list));
for i=1:length(lambda_name_list)

    %read lambda
    lambda_t = eval(lambda_name_list(i)+"(1,:)");%pixel0_alambdas100(1,:);
    tdc_width = uint32(length(lambda_t)/TDC_bins);
    hist_t = eval(lvar_name_list(i)+"(1,:)");%pixel0_all_ar_tof100(1,:);
    %creat system struct
    system.n_spad = n_spad_per_pix;
    system.n_pix = 5;
    system.n_tdc = 10;
    system.sim_res = 10e-12;
    system.arrival = "all";
    system.tdc_dead = 1e-12;
    system.spad_dead = 23e-9;
    system.c_f = c_f;
    system.c_t = c_t;
    system.adjust_cd = 4;
    
    %run analytical
    analitycal_pdf = analytical_model(lambda_t,system,transient_active);
    
    %we should match it with TDCs bin width
    %So here we do slice and sum
    analytical_hist = zeros(1,TDC_bins);
    for tau=1:TDC_bins
        idx = (tau-1)*tdc_width +1 ;
        idx_end = idx + tdc_width;
        if(idx_end>length(analitycal_pdf)) 
                idx_end = length(analitycal_pdf);   
        end
        analytical_hist(tau) = sim_res*sum(analitycal_pdf(idx:idx_end));
    end
    % 
    timeline = downsample(0:sim_res:(length(lambda_t)-1)*sim_res,tdc_width);
    %also the actual hist is actually upsampled in motecarlo by repeating
    %samples. -> revert it
    norm_hist = downsample(double(hist_t),tdc_width); %n
    %Quantize analytical
    norm_analytical =  double(uint32(N*analytical_hist));%sum(norm_hist)%.*normalize(analytical_hist,'norm',1); %N.*analytical_hist; %
    
    %plot :)
    f = figure; bar(timeline,norm_hist)
    hold on; plot(timeline, norm_analytical, LineWidth=2,LineStyle="--");
    fontsize(f,14,"points")
    xlim([0, 666e-9])
    xlabel("time (s)");
    ylabel("Normalized histogram")
    legend('Histogram', 'Analycital')
    savefig(sim_path+lvar_name_list(i)+"MATLAB.fig")
    set(gcf,'Units','inches');
    screenposition = get(gcf,'Position');
    set(gcf,...
        'PaperPosition',[0 0 screenposition(3:4)],...
        'PaperSize',[screenposition(3:4)]);
    saveas(gcf, sim_path+lvar_name_list(i)+"MATLAB.pdf")
    saveas(gcf, sim_path+lvar_name_list(i)+"MATLAB.png")

    ds_tofs = uint32(true_tof_list/tdc_width);
    
    %SNR Calc
    peak_A = max(norm_analytical(ds_tofs(i)-2:ds_tofs(i)+2));
    valley_A = mean(norm_analytical(1:ds_tofs(i)-3));
    SNR_A(i) = (peak_A-valley_A)./ sqrt(peak_A);

    peak_MC = max(norm_hist(ds_tofs(i)-2:ds_tofs(i)+2));
    valley_MC = mean(norm_hist(1:ds_tofs(i)-3));
    SNR_MC(i) = double(peak_MC-valley_MC)./ sqrt(double(peak_MC));


    
    
end
    %SNR plot
    f=figure; semilogy(double(true_tof_list)*sim_res,SNR_MC,double(true_tof_list)*sim_res,SNR_A, LineWidth=2)
    fontsize(f,14,"points")
    xlabel("ToF (s)");
    ylabel("SNR")
    legend('MonteCarlo', 'Analycital')
    saveas(gcf, sim_path+"SNR_CMP_"+"MATLAB.pdf")
    saveas(gcf, sim_path+"SNR_CMP_"+"MATLAB.png")

