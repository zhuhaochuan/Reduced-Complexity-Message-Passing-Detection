
function [out1,out2,out3,out4,out5] = MPD_256qam_G_J(K,J,Z,N0v,s,t,Es,E_guiyi)
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
youxiao=zeros(2*K,8);%取8个中的4个
index=zeros(2*K,8);%对应取的概率的位置
index_=zeros(2*K,16);
L=zeros(2*K,cs);

for t=1% 迭代的循环 
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
    end   
   for i_c=1:2*K
       for j_c=1:cs
           if L(i_c,j_c)>709
              L(i_c,:)=L(i_c,:).*0.5;
           end
       end
   end  
    LL(:,:,t)=L;
  %这部分只是计算出概率的部分
    for k=1:2*K
       for n=1:cs
             pro_(k,n)=exp(L(k,n))/(sum(exp(L(k,:))));%更新的概率
             pro(k,n)=(1-D).* pro_(k,n)+D.*pro(k,n);%加上阻尼因子
       end 
    end
dd(:,t)=d;
pp(:,:,t) =pro;  
end


%第一次迭代正常全部更新 第二次第三次取两个点更新
for t=2:6% 迭代的循环 
     for k_i=1:2*K
        [paixu(k_i,:),index_(k_i,:)]=sort(pro(k_i,:));%sort从小到大排序
        youxiao(k_i,:)=paixu(k_i,9:16);
        index(k_i,:)=index_(k_i,9:16);
    end
    for i_=1:2*K%数据xi的循环
         for j_=1:2*K%求均值方差的通项的循环
              a(j_) =J(i_,j_) *(sym_(index(j_,:))*pro(j_,index(j_,:))');%求xj得均值不包含xi的所有xi
              b(j_)=(J(i_,j_).^2)*((sym_(index(j_,:)).^2)*pro(j_,index(j_,:))'-(abs(sym_(index(j_,:)))*pro(j_,index(j_,:))')^2); 
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
    end     
   for i_c=1:2*K
       for j_c=1:cs
           if L(i_c,j_c)>709
              L(i_c,:)=L(i_c,:).*0.5;
           end
       end
   end    
    LL(:,:,t)=L;
  %这部分只是计算出概率的部分
    for k=1:2*K
       for n=1:cs
             pro_(k,n)=exp(L(k,n))/(sum(exp(L(k,:))));%更新的概率
             pro(k,n)=(1-D).* pro_(k,n)+D.*pro(k,n);%加上阻尼因子
       end 
    end
dd(:,t)=d;
pp(:,:,t) =pro;  
end
%第四次迭代开始不再更新概率点的位置 而是固定概率点继续更新 这样带来的误差就是如果后期发生概率值得跳变将会导致错误的结果
for t=7:s% 迭代的循环 
    for i_=1:2*K%数据xi的循环
        for j_=1:2*K%求均值方差的通项的循环
              a(j_) =J(i_,j_) *(sym_(index(j_,:))*pro(j_,index(j_,:))');%求xj得均值不包含xi的所有xi
              b(j_)=(J(i_,j_).^2)*((sym_(index(j_,:)).^2)*pro(j_,index(j_,:))'-(abs(sym_(index(j_,:)))*pro(j_,index(j_,:))')^2);         
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
    end     
   for i_c=1:2*K
       for j_c=1:cs
           if L(i_c,j_c)>709
              L(i_c,:)=L(i_c,:).*0.5;
           end
       end
   end  
    LL(:,:,t)=L;
  %这部分只是计算出概率的部分
  %现在这部分只是更新所固定的两种符号概率 其他的都不再计算 这样是又进一步减少了存储和计算复杂度
    for k=1:2*K
       for n=index(k,1):index(k,2)
             pro_(k,n)=exp(L(k,n))/(sum(exp(L(k,:))));%更新的概率
             pro(k,n)=(1-D).* pro_(k,n)+D.*pro(k,n);%加上阻尼因子
       end 
    end
dd(:,t)=d;
pp(:,:,t) =pro;  
end
out1=L;
out2=pp;
out3=dd;
out4=bb;
out5=LL;


