function [data_cat,IKK_measure,FWHM,counts,spot_position_bp,spot_position_gauss] = IKK_master(protein_dir,filename_bi,reference_spot,red_dist_lim,green_dist_lim)

[data_cat] = IKK_matfile_combine(protein_dir,filename_bi);

[IKK_measure,FWHM,counts] = IKK_run_ver3(data_cat,reference_spot,red_dist_lim,green_dist_lim);

[spot_position_bp] = IKK_dist_ver3(IKK_measure,'bp',reference_spot);

[spot_position_gauss] = IKK_dist_ver3(IKK_measure,'gauss',reference_spot);