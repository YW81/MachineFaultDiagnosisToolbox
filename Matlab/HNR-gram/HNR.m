function [HNR1] = HNR (y, fs, noOfFrames)

    %find the maximum lag M
    M = fs*2; %25 is the default minimum fundatmental frequency
    % IS IT? FOR WHAT SAMPLE RATE? HOW DO YOU DETERMIMNE THIS?
    for i = 1:(noOfFrames - 1)
       NA = xcorr(y(:,i), y(:,i), M, 'coeff');
       NA = NA( ceil(length(NA)/2) : end );
       % find first zero crossing
       sample1 = NA(1);

       for lag = 2:length(NA)
           sample2 = NA(lag);
           if(( sample1 > 0 ) && ( sample2 < 0 ))
               zero_position = lag;
               break;
           elseif(( sample1 == 0 ) || ( sample2 == 0))
               zero_position = lag;
               break;
           else
               sample1 = sample2;
           end
       end
       
       NA = NA( zero_position : end );
       [max_value max_position] = max(NA);

%        [max_value max_position] = max(abs(NA));

       HR(i) = max_value;
 
%        HNR1(i)= 10*log10(HR(i)/(1-HR(i)));
    
      HNR1(i)= (HR(i)/(1-HR(i)));
      
%       HNR1(i) =  HR(i);
       
%         HNR1(i)= HR(i);
    end
    
    
    
    
    
     
    
    
    
