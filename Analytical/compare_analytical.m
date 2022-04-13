close all
c_t = 3e-9;
c_f = 16;
transient_active = 1;
n_spad_per_pix = 133;
warning('off','MATLAB:nchoosek:LargeCoefficient')
sim_res = 10e-12;
sim_path = "./outputs/milestone/sensl/100k/cd16/";
load(sim_path+"lambdas_strt100_stp500.mat");
load(sim_path+"N1000_strt100_stp500.mat");
lambda_name_list = ["pixel0_alambdas100", "pixel1_alambdas200", "pixel2_alambdas300", "pixel3_alambdas399", "pixel4_alambdas500"];
lvar_name_list = ["pixel0_all_ar_tof100","pixel1_all_ar_tof200","pixel2_all_ar_tof300","pixel3_all_ar_tof399","pixel4_all_ar_tof500"];
true_tof_list  = uint32([100e-9, 200e-9, 300e-9, 400e-9, 500e-9]/sim_res);
SNR_A = zeros(size(true_tof_list));
SNR_MC = zeros(size(true_tof_list));
for i=1:length(lambda_name_list)

    lambda_t = eval(lambda_name_list(i)+"(1,:)");%pixel0_alambdas100(1,:);
    hist_t = eval(lvar_name_list(i)+"(1,:)");%pixel0_all_ar_tof100(1,:);
    
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
    
    analitycal_pdf = analytical_model(lambda_t,system,transient_active);
    
    % analitycal_pdf  = analytical_lidar(lambda_t,c_f,c_t,n_spad_per_pix,sim_res);
    % 
    timeline = 0:sim_res:(length(lambda_t)-1)*sim_res;
    norm_analytical = normalize(analitycal_pdf,'norm',1.00);
    norm_hist = normalize(double(hist_t), 'norm',1);
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

    peak_A = max(norm_analytical(true_tof_list(i)-300:true_tof_list(i)+300));
    valley_A = mean(norm_analytical(1:true_tof_list(i)-600));
    SNR_A(i) = (peak_A-valley_A)./ sqrt(peak_A);

    peak_MC = max(norm_hist(true_tof_list(i)-300:true_tof_list(i)+300));
    valley_MC = mean(norm_hist(1:true_tof_list(i)-600));
    SNR_MC(i) = (peak_MC-valley_MC)./ sqrt(peak_MC);


    
    
end

    f=figure; semilogy(double(true_tof_list)*sim_res,SNR_MC,double(true_tof_list)*sim_res,1*SNR_A, LineWidth=2)
    fontsize(f,14,"points")
    xlabel("ToF (s)");
    ylabel("SNR")
    legend('MonteCarlo', 'Analycital')
    saveas(gcf, sim_path+"SNR_CMP_"+"MATLAB.pdf")
    saveas(gcf, sim_path+"SNR_CMP_"+"MATLAB.png")

