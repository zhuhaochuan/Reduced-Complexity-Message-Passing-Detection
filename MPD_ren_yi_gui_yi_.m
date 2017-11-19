%MIMO MPD 
clear all;
close all;
%做归一化

%输入参数的设置 
K =16;%用户天线数
N =128;%基站天线数
t=2;%选择调制阶数，按照下面的M中的顺序
s = 7;%算法迭代的次数
SNRdBs =6;%信噪比SNR的范围
MAX_nobit=100;%最大误码数

%基本参数
M_tiaozhi=[4,16,64,256];%调制阶数
M_cs=[2,4,6,8];%不同调制阶数对应不同的参数
M=M_tiaozhi(t);%调制阶数
can_s=M_cs(t);%调制阶数对应的参数
Es_=[2,10,42,170];%每个发送符号的能量
Es=Es_(t);%不同调制阶数下符号的能量
E_guiyi=sqrt(Es);
syn_per_ant = 360 ;%每个天线发送的数据长度
Num_pack =1 ;%发送的数据包个数
bit_len = syn_per_ant * K * Num_pack ;% 总的发射数据量=每个天线发射的数据包长度 * 天线数 * 发射的包数
%不同调制阶数下对于符号空间sym_的选择
if t==1
    sym_=[-1,+1]/E_guiyi;
elseif t==2
        sym_=[-3:2:3]/E_guiyi;
    elseif t==3
        sym_=[-7:2:7]/E_guiyi;
else 
       sym_=[-15:2:15]/E_guiyi;
end


fprintf('\tSNRdBs\t\tBER\n');%打印出每次信噪比下的误码率


tic%计算整个过程的时间
for i_SNR = 1:length(SNRdBs)  %对于信噪比变化的大循环
  
    SNRdB = SNRdBs(i_SNR) ;%信噪比
        rand('seed',1);
	    randn('seed',1);  
    nobit = 0 ;%误码数的初始化
    
    x_guji =zeros(2*K,1);
    
    x_bit_guji = zeros(K,1);
    
    count=0;%计数器的初始化
     
   while(1)
       
       count=count+1;%计数器加1
       
       T_bits = randi([0 1],bit_len, 1);%发送的bit数据
   
       receive_bit =zeros(length(bit_len),1);%接受bit数据的初始化
            
     %输入数据并检测的循环
      for tx_time = 1 : (syn_per_ant * Num_pack/can_s) %每次进入8K个bit数据，这是传输次数的循环
          
            tx_bit = T_bits( (( tx_time - 1 ) *can_s* K + 1 ): tx_time *can_s * K );%每次传输的8K个bit数据
          
            mod = modem.qammod('M',M,'InputType','Bit','SymbolOrder','gray');
       
	        x_1 = modulate(mod,tx_bit);%QAM调制之后的符号数据（复数）
            
            x_ = x_1./E_guiyi;%能量归一化
            
            N0 = (K/2)*1*10^(-SNRdB/10);%噪声方差Es=1 归一化后
            
             sym_shi_shu=zeros(2*K,1);
             
            %构建信道矩阵并转化成实数形式
            H = zeros(2*N,2*K);
             hseed = sqrt(0.5)* (randn(N,K) + 1j * randn(N,K));
             hseed1 = real(hseed);
	         hseed2 = imag(hseed);
	        for si1 = 1 : N
	    	for sj1 = 1 : K
	    				hi1 = 2 * si1 - 1;
	    				hj1 = 2 * sj1 - 1;
	    				H(hi1,hj1) = hseed1(si1,sj1);
	    				H(hi1+1, hj1+1) = hseed1(si1,sj1);
            end
	        end
	    		    for si2 = 1 : N
	    			for sj2 = 1 : K
	    				hi2 = 2 * si2;
	    				hj2 = 2 * sj2 - 1;
	    				hi3 = 2 * si2 - 1;
	    				hj3 = 2 * sj2;
	    				H(hi2,hj2) = hseed2(si2,sj2);
	    				H(hi3, hj3) = -hseed2(si2,sj2);
	    			end
                    end  
                    
                    %将调制之后的复数符号转化成实数形式，
                    sseed = x_;
	    		    sseed1 = real(sseed);
	    		    sseed2 = imag(sseed);
	    		   for si = 1 : K
	    			    sym_shi_shu(2*si-1) = sseed1(si);
	    			    sym_shi_shu(2*si) = sseed2(si);
                   end
                   x=sym_shi_shu;
                
               %噪声的建立
             noise = zeros(2*N,1);
             nseed = sqrt(N0)*(randn(N,1) + 1j * randn(N,1));
             nseed1 = real(nseed);
	    	 nseed2 = imag(nseed);
    	    	    	   for ni = 1 : N
                               noise(2*ni-1) = nseed1(ni);
                               noise(2*ni) = nseed2(ni);
                           end
                    
        %   接收到的数据 
             y = H * x + noise; 
              
        %信道硬化  近似 J,Z调用子函数
        [J,Z,N0v] = ESTIMATE(N0,N,K,y,1,H);%归一化后Es=1
      
        %MPD迭代算法求出 x的LLR估计, s 为迭代次数，调用子函数
         [L,pp,dd,bb,LL]= MPD_ren_yi_gui_yi(K,J,Z,N0v,s,t,1,E_guiyi);
         %[L,pp,dd,bb,LL]= MPD_16qam_G_J(K,J,Z,N0v,s,t,1,E_guiyi);
         
        %算法输出是LLR 用LLR判决 
        %输出判决过程
        for n=1:2*K
            L_=max(L(n,:));
            index=find(L(n,:)==max(L(n,:)));
            if L_>0
                x_guji(n,1)=sym_(index);
            else
                x_guji(n,1)=sym_(1);
            end
        end
        
        
        %将实数符号化回复数形式
        for n_=1:K 
            x_bit_guji(n_) = x_guji(2*n_-1,1) + x_guji(2*n_,1)*1j;
        end
         x_bit_guji= x_bit_guji* E_guiyi;
       %解调过程
       demod = modem.qamdemod('M',M,'OutputType','Bit','DecisionType','hard decision','SymbolOrder','gray','NoiseVariance',N0);
        x_bit_1= demodulate(demod,x_bit_guji);%输出的解调的数据比特这代码和N0没有关系 
        x_bit_=x_bit_1(:)';
        receive_bit(can_s*(tx_time-1)*K+1:can_s*tx_time*K,1) =x_bit_;
      end    %到这里一次传输倒到接受估计结束，再次循环输入数据直到数据全部输入
      
      %计算误码数
      for i_comp=1:bit_len
           if T_bits(i_comp)~=receive_bit(i_comp) 
                nobit=nobit+1;   
           end
      end
           if nobit>MAX_nobit
              break;%如果误码数大于设定的最大误码数即误码率大于所设定的最大误码数对应的误码率则跳出whlie循环进行下一个SNR下BER的计算
           end 
           %if count>=100
              %break;
           %end
   end%结束while循环
   
   
   %误码率的计算
    BER(i_SNR) = nobit / (bit_len*count);
  
    fprintf('\t%d\t\t%e\t\n',SNRdB,BER(i_SNR));%打印出每次误码率的结果
             
end
toc%计算算法执行时间

 %   制作相应表格显示误比特率和信噪比的关系
    semilogy(SNRdBs(1:length(BER)),BER,'b<-','LineWidth',2);
    xlabel('Average SNR in dB');%横坐标
    ylabel('Uncoded BER');%纵坐标
    title('N=128,K=16,16QAM');   %标题
    legend('MPD');%标注
    grid on%加网格线