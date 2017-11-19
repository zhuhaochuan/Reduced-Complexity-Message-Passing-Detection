
function [out1,out2,out3,out4,out5,out6,out7,out8] = dynamic_64(K,J,Z,N0v,s,t,Es,E_guiyi)
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
        sym_=[-7:2:7,0]/E_guiyi;
        p_cs=0.125;
        cs=8;
        else 
       sym_=[-15:2:15]/E_guiyi;
       p_cs=0.0625;
       cs=16;
end

pro=(p_cs)*ones(2*K,cs+1);%初始化 每一个符号的概率为1、16
pro(:,cs+1)=0;
L=zeros(2*K,cs);
paixu=zeros(2*K,8);%对于64qam而言

index=zeros(2*K,4);%对应取的概率的位置
index_=zeros(2*K,8);
L=zeros(2*K,cs);
dy = zeros(2*K,4,13);
for t=1% 迭代的循环 
    for i_=1:2*K%数据xi的循环
         for j_=1:2*K%求均值方差的通项的循环
              a(j_) =J(i_,j_) *(sym_*pro(j_,:)');%求xj得均值不包含xi的所有xi
              b(j_)=(J(i_,j_).^2)*((sym_.^2)*pro(j_,:)'-(abs(sym_*pro(j_,:)'))^2);
              %b(j_)=(J(i_,j_).^2)*(Es/2);%方差不更新
         end
         
         bb(:,i_,t)=b;
         
       %求和
      
         c(i_)= sum(a(:))-a(i_);
         d(i_)= sum(b(:)) - b(i_) + N0v;
         
      
        %到这里一次迭代中所要的消息已经更新完毕 
         %求每个符号的对数似然比
         for n_=1:cs
            % L_(i_,n_)=(2*J(i_,i_)*(Z(i_)-c(i_))*(sym_(n_)-sym_(1))+(J(i_,i_)^2)*(sym_(1)^2-sym_(n_)^2))/(2*d(i_));
              L(i_,n_)=(J(i_,i_)*(sym_(n_)-sym_(1)))*(2*Z(i_)-2*c(i_)-J(i_,i_)*sym_(n_)-J(i_,i_)*sym_(1))/(2*d(i_));
             % L(i_,n_)=(1-D)* L_(i_,n_)+D*L(i_,n_);
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
for t=2:s% 迭代的循环 
    %这是对上次迭代结果的概率排序的确对每一个符号排序现在考虑增加条件判断
     for k_i=1:2*K%动态的决定要取几个
        [paixu(k_i,:),index_(k_i,:)]=sort(pro(k_i,1:cs));%sort从小到大排序
        
        %这么写出现的都是最后一个点的情况 肯定不对的所以还得改
        %这里的动态标准是 出现0.5-0.9之间的就取2个 出现0.9以上的就缩小到1个其余4个
        %其余情况都是取8个这样动态的选择会进一步降低计算的复杂度
        for j_k=1:cs%逐个判断只要出现在范围内的就跳出当前循环执行下次循环否则会出现上面提到的问题
            if (pro(k_i,j_k)>0.5)&&(pro(k_i,j_k)<0.9)%大于0.5取2个
                index(k_i,1:2)=index_(k_i,7:8);
                index(k_i,3:4)=9;
                break;
           
            elseif  pro(k_i,j_k)>=0.9%大于0.9取1个
                
                index(k_i,1)=index_(k_i,8);
                index(k_i,2:4)=9;
                break;
            else%其他情况取4个
                 index(k_i,:)=index_(k_i,5:8);  
                 
            end
        end
               
     end
    in(:,:,t)=index;

    
 
    
    for i_=1:2*K%数据xi的循环
         for j_=1:2*K%求均值方差的通项的循环
              a(j_) =J(i_,j_) *(sym_(index(j_,:))*pro(j_,index(j_,:))');%求xj得均值不包含xi的所有xi
              b(j_)=(J(i_,j_).^2)*((sym_(index(j_,:)).^2)*pro(j_,index(j_,:))'-(abs(sym_(index(j_,:)))*pro(j_,index(j_,:))')^2);
              %b(j_)=(J(i_,j_).^2)*(Es/2);%方差不更新
         end
         
         bb(:,i_,t)=b;
         
       %求和
      
         c(i_)= sum(a(:))-a(i_);
         d(i_)= sum(b(:)) - b(i_) + N0v;
         
      
        %到这里一次迭代中所要的消息已经更新完毕 
         %求每个符号的对数似然比
         for n_=1:cs
            % L_(i_,n_)=(2*J(i_,i_)*(Z(i_)-c(i_))*(sym_(n_)-sym_(1))+(J(i_,i_)^2)*(sym_(1)^2-sym_(n_)^2))/(2*d(i_));
              L(i_,n_)=(J(i_,i_)*(sym_(n_)-sym_(1)))*(2*Z(i_)-2*c(i_)-J(i_,i_)*sym_(n_)-J(i_,i_)*sym_(1))/(2*d(i_));
             % L(i_,n_)=(1-D)* L_(i_,n_)+D*L(i_,n_);
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


    %动态更新的复杂度分析：
    dy = in;
     dy(dy~=9)= 1;
     dy(dy==9)= 0;
     num=sum(sum(sum(dy)));
    
out1=L;
out2=pp;
out3=dd;
out4=bb;
out5=LL;
out6=in;
out7=num;
out8=dy;


