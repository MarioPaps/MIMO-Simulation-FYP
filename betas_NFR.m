%output beta for MAI users for given NFR level
function [beta_complex] = betas_NFR(MAI_powers,M)

    beta_mag= sqrt(MAI_powers);
    beta_complex= beta_mag.*exp(1i*2*pi*rand(size(beta_mag)));

end