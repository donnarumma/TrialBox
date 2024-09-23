function   DCM_RESULTS_LFP_KS(DCM)
% function DCM_RESULTS_LFP_KS(DCM)
% ~/TESTS/SAPIENZA/DCM/
% LFP_M3_SK009_1  dir 1
% LFP_M3_SK009_9  dir 1 
% LFP_M3_SK009_17 dir 1
NConds=size(DCM.xU.X,2);
for iC = 1:NConds
    
    % extrinsic (driving) connections
    %----------------------------------------------------------------------           
    CM = DCM.Ep.A{3}+DCM.Ep.B{iC};
    name=DCM.xU.name(iC);

    figure; 
    % imagesc(exp(DCM.Ep.A{3}));
    imagesc(exp(CM));

    xticks(1:2)
    xticklabels({'S','K'});
    yticks(1:2)
    yticklabels({'S','K'});
    title(name)
    xlabel('from')
    ylabel('to');
    colormap gray;

end


plot_spm_dcm_csd_results(DCM,'Coupling (A)',figure);
% plot_spm_dcm_csd_results(DCM,'spectral data',figure);  
plot_spm_dcm_csd_results(DCM,'Coupling (B)',figure);
plot_spm_dcm_csd_results(DCM,'trial-specific effects',figure);

end