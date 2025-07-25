function res = enlarge(image,x,y)


%% 设置参数区——所有用到的参数都在这里设置===================================================
%========================================================================================
img = image;%加载要操作图像的地址，在当前文件夹只写文件名
pt = [x, y];%框的左上角起始坐标点
wSize = [50,50];%以pt为原点，截取终点的坐标
lineSize=3;%框线的宽度
flag=2;  %flag=1: 有缺口的框  %flag=2: 无缺口的框
color=[255 0 0];%框的颜色(RGB值)
magnification=4;%小图放大倍数
lineSize_1=3;%小图框线的宽度
%========================================================================================
%========================================================================================
 
%% 大图上画框
des = drawRect(img,pt,wSize,lineSize,color,flag);
 
%% 截取出小图
img_1 =imcrop(img,[pt(1),pt(2),wSize(1),wSize(2)]);  %切割图像，起始坐标点（x1,y1）截取到终止坐标点(x2,y2)
%% 小图放大后画框
img_1=imresize(img_1,magnification);
pt_1=[lineSize_1,lineSize_1];[M2,N2,p]=size(img_1);wSize_1=[N2-lineSize_1,M2-lineSize_1];
des_1 = drawRect(img_1,pt_1,wSize_1,lineSize_1,color,flag);[M2,N2,p]=size(des_1);
% figure; imshow(des_1);
[M1,N1,c]=size(img);%大框的尺寸
if M2>M1||N2>N1
        disp('小框的大小已经超过大框了！请减小小图的放大倍数');
        return;
end
%% 两图重叠
% %小图放到大图的左下角
% for i=1:M2
%     for j=1:N2
%         des((M1-M2+i),j,:)=des_1(i,j,:);
%     end
% end

% %右下角
for i=1:M2
    for j=1:N2
        des((M1-M2+i), (N1-N2+j), :) = des_1(i,j,:);
    end
end

%右上角
% for i=1:M2
%     for j=1:N2
%         des(i, (N1-N2+j), :) = des_1(i,j,:);
%     end
% end

res = des;

end



%===================================主函数结束=============================================================
 
%%
function [ dest ] = drawRect( src, pt, wSize,  lineSize, color,flag )
%简介：
% %将图像画上有颜色的框图，如果输入是灰度图，先转换为彩色图像，再画框图
% 图像矩阵
% 行向量方向  是  y
% 列向量方向  是  x
%----------------------------------------------------------------------
%输入：
% src：        原始图像，可以为灰度图，可为彩色图
% pt：         左上角坐标   [x1, y1]
% wSize：   框的大小      [wx, wy]
% lineSize： 线的宽度
% color：     线的颜色      [r,  g,  b] 
%----------------------------------------------------------------------
%输出：
% dest：           画好了的图像
%----------------------------------------------------------------------
 
%flag=1: 有缺口的框
%flag=2: 无缺口的框
 
 
%判断输入参数个数
if nargin < 6
   
end
 
if nargin < 5
    lineSize = 1;
end
 
if nargin < 4
    disp('输入参数不够 !!!');
    return;
end
 
%判断框的边界问题
[yA, xA, z] = size(src);
x1 = pt(1);
y1 = pt(2);
wx = wSize(1);
wy = wSize(2);
if  x1>xA || ...
        y1>yA||...
        (x1+wx)>xA||...
        (y1+wy)>yA
 
    disp('画的框将超过图像 !!!');
    return;
end
 
%如果是单通道的灰度图，转成3通道的图像
if 1==z
    dest(:, : ,1) = src;
    dest(:, : ,2) = src;
    dest(:, : ,3) = src;
else
    dest = src;
end
 
%开始画框图
for c = 1 : 3                 %3个通道，r，g，b分别画
    for dl = 1 : lineSize   %线的宽度，线条是向外面扩展的
        d = dl - 1;
        if  1==flag %有缺口的框
            dest(  y1-d ,            x1:(x1+wx) ,  c  ) =  color(c); %上方线条
            dest(  y1+wy+d ,     x1:(x1+wx) , c  ) =  color(c); %下方线条
            dest(  y1:(y1+wy) ,   x1-d ,           c  ) =  color(c); %左方线条
            dest(  y1:(y1+wy) ,   x1+wx+d ,    c  ) =  color(c); %左方线条
        elseif 2==flag %无缺口的框
            dest(  y1-d ,            (x1-d):(x1+wx+d) ,  c  ) =  color(c); %上方线条
            dest(  y1+wy+d ,    (x1-d):(x1+wx+d) ,  c  ) =  color(c); %下方线条
            dest(  (y1-d):(y1+wy+d) ,   x1-d ,           c  ) =  color(c); %左方线条
            dest(  (y1-d):(y1+wy+d) ,   x1+wx+d ,    c  ) =  color(c); %左方线条
        end
    end    
end %主循环尾
 
 
end %函数尾