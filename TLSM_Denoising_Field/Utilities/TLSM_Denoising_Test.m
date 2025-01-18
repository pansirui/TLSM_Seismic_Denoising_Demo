function [PSNR_Final,FSIM_Final,SSIM_Final,ERGAS_Final, SAM_Final, iter, Time_s, SNR_Final] = TLSM_Denoising_Test(ori, Y, tau, alpha, beta, gamma, lambda1, lambda2)
 
randn ('seed',0);

par = Parset_TLSM(tau, alpha, beta, gamma, lambda1, lambda2);

time0 = clock;

[Denoising, iter] = TLSM_Denoising(ori, Y, par);  

Time_s = (etime(clock,time0));  

Xnew = Denoising{iter}; 

output_image = ((Xnew - min(Xnew(:))) / (max(Xnew(:)) - min(Xnew(:))))*255;
Ori_Image = ((ori - min(ori(:))) / (max(ori(:)) - min(ori(:))))*255;

[PSNR_Final, SSIM_Final, FSIM_Final, ERGAS_Final, SAM_Final] = MSIQA(Ori_Image, output_image);

end
