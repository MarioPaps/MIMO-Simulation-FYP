function[RMSE]= findRMSE(theor,est)
            
    diff= (theor-est).^2;
    avg= sum(diff)/ length(diff);
    RMSE=sqrt(avg);
end