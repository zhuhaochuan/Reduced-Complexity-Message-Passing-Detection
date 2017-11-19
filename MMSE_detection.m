%   MMSE Detection
%   输入：   H:    channel matrix
%           y:     receive bits
%           sigma:   子信道的噪声方差
%           M:    发射天线数
%   输出：   x：    估计的输入数据
function x=MMSE_detection(H,r,sigma,M)
% function x=MMSE_detection(H,r,M)
    W = (H' * H) + sigma.*eye(M);
    b = H'*r;
    x = inv(W)*b;
   
