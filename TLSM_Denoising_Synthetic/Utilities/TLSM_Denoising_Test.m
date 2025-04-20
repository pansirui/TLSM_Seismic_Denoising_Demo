function [PSNR_Final,FSIM_Final,SSIM_Final,ERGAS_Final, SAM_Final, iter, Time_s] = TLSM_Denoising_Test(ori, Y, tau, alpha, beta, gamma, lambda1, lambda2)
 
randn ('seed',0);

% parameters setting
par = Parset_TLSM(tau, alpha, beta, gamma, lambda1, lambda2);

time0 = clock;

% TLSM for seismic noise suppression
[Denoising, iter] = TLSM_Denoising(ori, Y, par);  

Time_s = (etime(clock,time0));  

Xnew = Denoising{iter};   

% normalization
output_image = ((Xnew - min(Xnew(:))) / (max(Xnew(:)) - min(Xnew(:))))*255;
Ori_Image = ((ori - min(ori(:))) / (max(ori(:)) - min(ori(:))))*255;

% quality assess
[PSNR_Final, SSIM_Final, FSIM_Final, ERGAS_Final, SAM_Final] = MSIQA(Ori_Image, output_image);

end
