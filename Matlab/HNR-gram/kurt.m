function K = kurt(x,opt)%计算信号x的峭度

global fr;
global N;
global snr;
global fs;


xn = length(x);
ratio = N/xn;
newfs = fs/ratio;



if strcmp(opt,'kurt2')%判断opt是否为kurt2
   if all(x == 0), K = 0; E = 0;return;end 
   

   y = abs(hilbert(real(x)))-mean(abs(hilbert(real(x))));
%    y = abs(x)-mean(abs(x));


   x = x - mean(x);
   
%    E = mean(abs(x).^2);
%    if E < eps, K = 0; return;end
%    K = mean(abs(x).^4)/E^2;
% 
% 
%    if all(isreal(x))
%       K = K - 3;								% real signal
%    else
%       K = K - 2;
%    end

    K = HNR(y,round(newfs),2);


elseif strcmp(opt,'kurt1')
   if all(x == 0), K = 0; E = 0;return;end
   
      y = abs(hilbert(real(x)))-mean(abs(hilbert(real(x))));
   
   x = x - mean(x);
   E = mean(abs(x));
   if E < eps, K = 0; return;end
    K = mean(abs(x).^2)/E^2;

   if all(isreal(x))
      K = K - 1.57;							% real signal
   else
      K = K - 1.27;
   end

end




% ------------------------------------------------------------------------
