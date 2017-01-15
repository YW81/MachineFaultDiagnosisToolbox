function [Bw,fc] = HNR_gram(x,nlevel,Fs,plot)
% 输入：
%  x——原始信号
%  nlevel——分解层数
%  Fs——采样频率
%  plot——是否作图（为1时画图）
% 
% 输出：
% Bw ——最有滤波频带的带宽
% fc ——最有滤波频带的中心频率
% -------------------
% Origin from J. Antoni
% Xiaoqiang @2017.1
% -------------------

opt1 =1;
opt2 = 1;

N = length(x);
N2 = log2(N) - 7;%最大允许分解层
if nlevel > N2%判断输入分解层是否大于最大允许分解层
   error('Please enter a smaller number of decomposition levels');
end


% Fast computation of the kurtogram
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if opt1 == 1
   % 1) Filterbank-based kurtogram
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Analytic generating filters
   N = 16;			fc = .4;					% a short filter is just good enough!
   h = fir1(N,fc).*exp(2i*pi*(0:N)*.125);%二分段低通滤波器
   n = 2:N+1;
   g = h(1+mod(1-n,N)).*(-1).^(1-n);%二分段高通滤波器
   % 
   N = fix(3/2*N);
   h1 = fir1(N,2/3*fc).*exp(2i*pi*(0:N)*.25/3);%三分段第一段滤波器
   h2 = h1.*exp(2i*pi*(0:N)/6);%三分段第二段滤波器
   h3 = h1.*exp(2i*pi*(0:N)/3);%三分段第三段滤波器  
   % 
   if opt2 == 1
      Kwav = K_wpQ(x,h,g,h1,h2,h3,nlevel,'kurt2');				% kurtosis of the complex envelope
   else
      Kwav = K_wpQ(x,h,g,h1,h2,h3,nlevel,'kurt1');				% variance of the envelope magnitude
   end
   Kwav = Kwav.*(Kwav>0);												% keep positive values only!
   
else 
   % 2) STFT-based kurtogram
   %%%%%%%%%%%%%%%%%%%%%%%%%
   Nfft = 2.^[3:nlevel+2];				% level 1 of wav_kurt roughly corresponds to a 4-sample hanning window with stft_kurt
   %											  or a 8-sample flattop	第一层分8段	
   temp = [3*Nfft(1)/2 3*Nfft(1:end-2);Nfft(2:end)];%第一行为三分层段数，第二行为二分层段数
   Nfft = [Nfft(1) temp(:)'];%用矩阵计算而不是循环赋值的方法，所有分层段数
   if opt2 == 1
      Kstft = Kf_fft(x,Nfft,1,'kurt2');							% kurtosis of the complex envelope
      Kx = kurt(x,'kurt2');%未分层前的信号峭度
   else
      Kstft = Kf_fft(x,Nfft,1,'kurt1');							% variance of the envelope magnitude
      Kx = kurt(x,'kurt1');
   end
   Kstft = [Kx*ones(1,size(Kstft,2));Kstft];%所有峭度汇总
   Kstft = Kstft.*(Kstft>0);%只取峭度大于零的部分，峭度小于零的部分视为零		% keep positive values only!
   
end

% Graphical results
%%%%%%%%%%%%%%%%%%%

if plot ==1
    figure
    if opt1 == 1%分叉树滤波结果
        Level_w = 1:nlevel;	Level_w = [Level_w;Level_w+log2(3)-1];	Level_w = Level_w(:); Level_w = [0 Level_w(1:2*nlevel-1)'];%图形纵坐标
        freq_w = Fs*((0:3*2^nlevel-1)/(3*2^(nlevel+1)) + 1/(3*2^(2+nlevel)));%图形横坐标Fs/2*((0:3*2^nlevel)+0.5)/(3*2^nlevel)
        imagesc(freq_w,1:2*nlevel,Kwav),colorbar,[I,J,M] = max_IJ(Kwav);%绘图并求最大值
        xlabel('Frequency [Hz]'),set(gca,'ytick',1:2*nlevel,'yticklabel',round(Level_w*10)/10),ylabel('Level k')%标注坐标
        fi = (J-1)/3/2^(nlevel+1);   fi = fi + 2^(-2-Level_w(I));%中心频率
        if opt2 == 1
            %title(['fb-kurt.2 - K_{max}=',num2str(round(10*M)/10),' @ level ',num2str(fix(10*Level_w(I))/10),', Bw= ',num2str(Fs*2^-(Level_w(I)+1)),'Hz, f_c=',num2str(Fs*fi),'Hz'])
                    title(['HNR_{max}=',num2str(round(10*M)/10),' @ level ',num2str(fix(10*Level_w(I))/10),', Bw= ',num2str(Fs*2^-(Level_w(I)+1)),'Hz, f_c=',num2str(Fs*fi),'Hz'])
            %      title(['Periodicity_{max}=',num2str(round(10*M)/10),' @ level ',num2str(fix(10*Level_w(I))/10),', Bw= ',num2str(Fs*2^-(Level_w(I)+1)),'Hz, f_c=',num2str(Fs*fi),'Hz'])
            Bw= Fs*2^-(Level_w(I)+1);
            fc=Fs*fi;
        else
            title(['fb-kurt.1 - K_{max}=',num2str(round(10*M)/10),' @ level ',num2str(fix(10*Level_w(I))/10),', Bw= ',num2str(Fs*2^-(Level_w(I)+1)),'Hz, f_c=',num2str(Fs*fi),'Hz'])
        end
    else%短时傅利叶变换结果
        LNw_stft = [0 log2(Nfft)];%纵坐标
        freq_stft = Fs*((0:Nfft(end)/2-1)/Nfft(end) + 1/Nfft(end)/2);%横坐标
        %freq_stft = Fs*(0:Nfft(end)/2-1)/Nfft(end);
        imagesc(freq_stft,1:2*nlevel,Kstft),colorbar,[I,J,M] = max_IJ(Kstft);
        fi = (J-1)/Nfft(end);
        xlabel('frequency [Hz]'),set(gca,'ytick',1:2*nlevel,'yticklabel',round(LNw_stft*10)/10),ylabel('level: log2(Nw)')
        if opt2 == 1
            title(['stft-kurt.2 - K_{max}=',num2str(round(10*M)/10),' @ Nw=2^{',num2str(fix(10*LNw_stft(I))/10),'}, fc=',num2str(Fs*fi),'Hz'])
        else
            title(['stft-kurt.1 - K_{max}=',num2str(round(10*M)/10),' @ Nw=2^{',num2str(fix(10*LNw_stft(I))/10),'}, fc=',num2str(Fs*fi),'Hz'])
        end
    end
    
end


% Signal filtering
%%%%%%%%%%%%%%%%%%
%c = [];
% test = input('Do you want to filter out transient signals from the kurtogram (yes = 1 ; no = 0): ');
% while test == 1
%    fi = input(['	Enter the optimal carrier frequency (btw 0 and ',num2str(Fs/2),') where to filter the signal: ']);
%    fi = fi/Fs;
%    if opt1 == 1
%       lev = input(['	Enter the optimal level (btw 0 and ',num2str(nlevel),') where to filter the signal: ']);
%       if opt2 == 1
%          [c,Bw,fc] = Find_wav_kurt(x,h,g,h1,h2,h3,nlevel,lev,fi,'kurt2',Fs);
%       else
%          [c,Bw,fc] = Find_wav_kurt(x,h,g,h1,h2,h3,nlevel,lev,fi,'kurt1',Fs);
%       end
%    else
%       lev = input(['	Enter the optimal level (btw 0 and ',num2str(nlevel+2),') where to filter the signal: ']);
%       if opt2 == 1
%          [c,Nw,fc] = Find_stft_kurt(x,nlevel,lev,fi,'kurt2',Fs);
%       else
%          [c,Nw,fc] = Find_stft_kurt(x,nlevel,lev,fi,'kurt1',Fs);
%       end
%    end
%    test = input('Do you want to keep on filtering out transients (yes = 1 ; no = 0): ');%选择是否退出循环
% end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%