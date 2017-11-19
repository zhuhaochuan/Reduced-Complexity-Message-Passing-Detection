function [out1,out2,out3,out4,out5] = MPD_zhuge_genxin(K,J,Z,N0v,s,t,Es,E_guiyi)
%行对应2K个数据 列对应16个符号的概率
D=0.33;
if t==1
    sym_=[-1,+1]/E_guiyi;
    p_cs=0.5;
    cs=2;
elseif t==2
        sym_=[-3:2:3]/E_guiyi;
        p_cs=0.25;
        cs=4;
    elseif t==3
        sym_=[-7:2:7]/E_guiyi;
        p_cs=0.125;
        cs=8;
        else 
       sym_=[-15:2:15]/E_guiyi;
       p_cs=0.0625;
       cs=16;
end

pro=(p_cs)*ones(2*K,cs);%初始化 每一个符号的概率为1、16

L=zeros(2*K,cs);

for t=1:s% 迭代的循环 
    for i_=1:2*K%数据xi的循环
         for j_=1:2*K%求均值方差的通项的循环
              a(j_) =J(i_,j_) *(sym_*pro(j_,:)');%求xj得均值不包含xi的所有xi
              b(j_)=(J(i_,j_).^2)*((sym_.^2)*pro(j_,:)'-(abs(sym_*pro(j_,:)'))^2);    
         end
         bb(:,i_,t)=b;
       %求和
         c(i_)= sum(a(:))-a(i_);
         d(i_)= sum(b(:)) - b(i_) + N0v;
        %到这里一次迭代中所要的消息已经更新完毕 
         %求每个符号的对数似然比
         for n_=1:cs
              L(i_,n_)=(J(i_,i_)*(sym_(n_)-sym_(1)))*(2*Z(i_)-2*c(i_)-J(i_,i_)*sym_(n_)-J(i_,i_)*sym_(1))/(2*d(i_));
         end 
         for i_c=1:cs
              if L(i_,i_c)>709
                 L(i_,:)=L(i_,:).*0.5;
              end
         end
          for n=1:cs
             pro_(i_,n)=exp(L(i_,n))/(sum(exp(L(i_,:))));%更新的概率
             pro(i_,n)=(1-D).* pro_(i_,n)+D.*pro(i_,n);%加上阻尼因子
         end 
    end 
      %防止方差太小 导致llr的值大于709 会使e指数值超过MATLAB计算上限     
    LL(:,:,t)=L;
  %这部分只是计算出概率的部分
dd(:,t)=d;
pp(:,:,t) =pro;  
end

out1=L;
out2=pp;
out3=dd;
out4=bb;
out5=LL;


