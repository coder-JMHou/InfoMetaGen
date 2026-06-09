%-------------------------------------------------
%创建衍射运算矩阵A
%X(a) = %%x(n)exp(-jkr)  r表示运算点和各个矩阵单元的距离
%输入      矩阵单元阵列个数和矩阵的长宽
%输入      频率，计算的距离
%num_l     矩阵单元X方向个数
%num_w     矩阵单元Y方向个数
%p    矩阵单元长度，假设长和宽一样
%lambda    波长
%d_dis     平面的距离
%mode_select  mode_select = 0 代表分母为1，mode_select = 1，分母带上距离的约束
%dis 运算方向  0 --》Y 1--》 X
%设计假设矩阵单元从上到下，从左到右进行排列成一竖
%-------------------------------------------------

function [metrix,metrix_inv] = Creat_metrix(num_x,num_y,p,lambda,distance,mode_select,dis)
	%计算出各个单元中心对应的坐标，假设以原点为中心,从上到下，从左到右
	X_lable_max=(num_x-1) * p/2;
	Y_lable_max=(num_y-1) * p/2;
    
    X_lable= - X_lable_max: p : X_lable_max;
    Y_lable= Y_lable_max:-p:-Y_lable_max;
	index = 1;
	total_unit_num = num_x * num_y;
	lable=zeros(total_unit_num,2); % 用于存放坐标值

	for i=1:num_x
		for j=1:num_y
            if dis == 0    %先从上到下，再从左到右
                lable(index,:) = [X_lable(i),Y_lable(j)];%第m个单元坐标 X_lable = lable(m,1)  Y_lable = lable(m,2)，从上到下从左到右存储坐标
            % elseif dis == 1  %X方向变化快
            %     lable(index,:) = [X_lable(j),Y_lable(i)];%第m个单元坐标 X_lable = lable(m,1)  Y_lable = lable(m,2)              
            end
            index = index +1;
		end
	end 
	
	%计算对应点的映射矩阵
	%假设为第m个对应点
	%与第一个矩形单元对应的距离为
	%dis_r(m,1) = sqrt(d_dis^2+ (lable(m,1) - lable(1,1)）^2 + (lable(m,2) - lable(1,2))^2)
	%与第k个矩形单元对应距离 dis_r(m,k) = sqrt(d_dis^2+ (lable(m,1) - lable(k,1)）^2 + (lable(m,2) - lable(k,2))^2)
	
	%计算距离矩阵
	dis_r=zeros(total_unit_num,total_unit_num);%用于存放坐标值
	%phi=zeros(total_unit_num,total_unit_num);%OAM
    
   for m=1:total_unit_num %变化的是映射单元，从第一个算起，变化的慢
        for k=1:total_unit_num %变化的是矩阵单元，从第一个算起，变化的快		
			dis_r(m,k) = sqrt(distance^2+ (lable(m,1) - lable(k,1))^2 + (lable(m,2) - lable(k,2))^2);
        end
   end 
    
	k = 2*pi/lambda; 
    
	if mode_select == 0
        metrix = exp((-1i*k).*dis_r);
        metrix_inv = exp((1i*k).*dis_r);
        %metrix = exp((-1i*k).*dis_r+li*OAM_mode*phi);
        %metrix_inv = exp((1i*k).*dis_r-li*OAM_mode*phi);
    else
        metrix = exp((-1i*k).*dis_r)./dis_r;
        metrix_inv = exp((1i*k).*dis_r)./dis_r./dis_r;
	end
end