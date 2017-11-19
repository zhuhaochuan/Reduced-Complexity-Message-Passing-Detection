%   uncoded information are transmitted on the channel
%   16QAM modulation
%   rayleigh fade channel
%   AWGN, i.i.d., ~CN(0,N0)
%   H i.i,d, ~CN(0,1)
%   SNR = U * Ex/N0
%   2017/04/19
    clear;
    close all;

%   configuration
    B = 128;
    U = 8;

    %   QAM64 configuration 
    M = 16;
    Ex = 1;
    E = sqrt(10);
	    MAXnobit = 50;

    k = log2 (M);
    SNRdBs = 0:14;			% although 16dB promised the performance of exact-float stuff, we should take 17dB into consideration
    data_per_ant = 300;
    raw_bit_len = data_per_ant * U;	% one antenna transmit 16*300=4800 bits per iteration

    K=4;				% the design that we planned is iteration 4 times

    EB = zeros(length(SNRdBs),1);
    BER = zeros(length(SNRdBs),1);	
    %   start transmission
    for i_SNR = 1:length(SNRdBs)
            SNRdB = SNRdBs(i_SNR);
	    rand('seed',1);
	    randn('seed',1);    
	    m = 0;
	    while (1)
		    m = m+1;
	    	%   constant transmit original data generation	    
	    	    b = randi([0 1], raw_bit_len, 1);  
	    	%   modulation
	    	    mod = modem.qammod('M',M,'InputType','Bit','SymbolOrder','gray');
	    	    xs = modulate(mod,b);
		    x = xs./E;
	    	    s = zeros(U,1);
		    rs = zeros(size(x));
		    [i,j] = size(x);	    
	    	%   transmission energy
	    	    N0 = U * Ex * 10^(-SNRdB/10);                  
	    	    for n = 1 : (i/U)
	    		    s = x((n-1)*U+1 : n*U);
	    		    H = sqrt(0.5)*(randn(B,U) + 1j * randn(B,U));
			    noise = sqrt(N0/2)*(randn(B,1) + 1j * randn(B,1));
    	                    y = H * s + noise;
			    rs ((n-1)*U+1 : n*U,1) = MMSE_detection(H,y,N0/Ex,U);
	            end
		%   demodulation
		    r = rs.*E;
	            demod = modem.qamdemod('M',M,'OutputType','Bit','DecisionType','hard decision','SymbolOrder','gray','NoiseVariance',N0);
	            d = demodulate(demod,r);
		%   calculate BER value
		    [e(m),ratio] = biterr(d,b);
		%   fprintf ('wrong-number = %d\n',e(m));
		    EB(i_SNR) = EB(i_SNR) + e(m);
		    if (EB(i_SNR)>MAXnobit) 
			    break;
		    end
	    end
	    BER(i_SNR) = EB(i_SNR)/(raw_bit_len * m);
            fprintf('SNR = %d\t\t number = %d\t\t BER = %e\t\t m = %d\n',SNRdB,EB(i_SNR),BER(i_SNR),m);
            semilogy(SNRdBs,BER,'rx-','LineWidth',2); 
    end


