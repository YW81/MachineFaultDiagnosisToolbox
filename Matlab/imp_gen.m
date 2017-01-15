% 轴承冲击信号模拟器
% fs-----信号的采样频率，Hz
% phs----信号的瞬时相位，度
% m------二阶系统的质量，千克
% c------二阶系统的阻尼，N/(m/s)
% k------二阶系统的刚度，N/m
% ang----故障的角度间隔
% gain---响应信号的幅值增益
% fault_type --故障类型，字符串型，分为：'none'，'inner'，'outer','rolling','rand_impact'


function imp = imp_gen(fs,phs,m,c,k,ang,gain,fault_type)

% 如果没有故障，直接置零
if strcmp(fault_type,'none')
    imp = 0*phs;
    return
end

N = length(phs);
    
% 计算冲击的固有频率   
fc = (k/m)^0.5/2/pi;
fd = fc*(1-(c/2/m/fc/2/pi)^2)^0.5;
% disp(['-----------------------------------------']);
% disp([fault_type, '故障固有频率']);
% disp(['无阻尼固有频率： ',num2str(fc),' Hz' ]);
% disp(['有阻尼固有频率： ',int2str(fd),' Hz' ]);
% disp(['-----------------------------------------']);

% 计算冲击信号的长度(根据信号的衰减率确定，可参考振动力学)
pn = round(20*m/c*fs);
temp = gain*impulse(tf(1,[m,c,k]),(0:pn-1)/fs);
imp = zeros(size(phs));


if strcmp(fault_type,'rand_impact')

    % 设置随机冲击发生的次数
    N_impact = 3;
    
    % 随机选取冲击发生的序号
    index = round(rand(1,N_impact)*(N-pn-1)); 
    
    % 由低到高排序
    index = sort(index);
    
    for i =1:N_impact
       % 在序号的位置加入冲击 
       imp(index(i):index(i)+pn-1) =  imp(index(i):index(i)+pn-1)+(0.5+rand(1)/2)*gain*temp;
        
    end
        
    return
end


% 按角度间隔在时域信号中添加冲击
for i = 1:N-pn
    if round(phs(i+1)/ang)>round(phs(i)/ang)
        imp(i:i+pn-1) =  imp(i:i+pn-1)+temp;
    end
end


% 对冲击序列进行幅值调试，模拟传递路径变化和载荷变化对内圈故障的影响
if strcmp(fault_type,'inner')
    am = 0.5*(1-1*cos(phs/180*pi));
elseif strcmp(fault_type,'outer')
    am = ones(size(imp));
elseif strcmp(fault_type,'rolling')
    am = ones(size(imp));
end
imp = imp.*am;

end