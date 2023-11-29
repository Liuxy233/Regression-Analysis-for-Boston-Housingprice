install.packages("mlbench") # 安装"mlbench"包。只需在计算机上安装一次即可反复调用，不必每次重新安装
install.packages("zoo")
library(mlbench) # 调用"mlbench"包
library(lmtest) # 调用"lmtest"包
library (MASS) # 调用"MASS"包
data("BostonHousing") # 载入"BostonHousing"数据
? BostonHousing # 查看BonstonHousing数据基本信息
# 在分析中，请以房价中位数"medv"为应变量，以其它十三个变量为自变量。其中"chas"为定性变量，不在课程范围内，
# 可删去，只考虑其它十二个自变量。
data1=BostonHousing
data1=data1[,-4]# 删去"chas"这一列
names(data1)[13]=c("y")
#names(data1)[1:13]=c("x1","x2","x3","x4","x5","x6","x7","x8","x9","x10","x11","x12","y")# 给自变量和因变量改名
data0=data.frame(scale(data1))
n = nrow(data0) # 清点样本量
p = ncol(data0) - 1 # 清点自变量个数（减去数据中包含的因变量一列）

lm1 = lm(y ~ ., data = data0) # 拟合线性模型
summary(lm1) # 输出拟合结果

e = resid(lm1) # 计算残差
sigma = sqrt(sum(e ^ 2) / (n - p - 1)) # 估计随机误差标准差sigma

X1 = as.matrix(data0[, 1 : 12]) # 抽取数据中的自变量矩阵
X = cbind(1, X1) # 在自变量矩阵右侧加一行常数1
XX = t(X) %*% X # 计算X'X
covBeta = sigma ^ 2 * solve(XX) # 计算最小二乘估计betahat的协方差矩阵
covBeta # 输出OLSE的协方差矩阵

r = matrix(nrow = p + 1, ncol = p + 1) # 定义一个零矩阵用于储存相关系数
for(i in 1 : (p + 1))
{
  for(j in 1 : (p + 1)) r[i, j] = covBeta[i, j] / sqrt(covBeta[i, i] * covBeta[j, j]) # 计算相关系数并储存在矩阵r中
}
r # 输出OLSE的相关系数矩阵

###-----变量选择——逐步回归-----###
summary(step(lm1, direction = "both"))
#indus和age变量没有被选择，而且此时的adjusted R-square比原模型高，所以删去这两个变量
data0=data0[,c(1,2,4,5,7,8,9,10,11,12,13)]
lm2=lm(y ~ .,data=data0)
lm2
#(appendix)后退法可得出同样结果
#Anova(lm1, type = "III") # 利用Anova函数对每一个自变量做F检验，输出结果

#lm_dropage = lm(y ~ . - age, data = data0) # 后退法，先去除Anova结果中最不显著的x6，用y对剩余自变量拟合线性模型
#summary(lm_dropage) # 输出拟合结果

#lm_droptwo = lm(y ~ .- age-indus, data = data0) # 再去除Anova结果中最不显著的x3，用y对剩余自变量拟合线性模型
#summary(lm_droptwo) # 输出拟合结果
# 用后退法发现删去indus和age后其他自变量都显著

######回归诊断######
###-----计算各样本点的Cook距离与异常点F检验p值-----###

ck=cooks.distance(lm2) # 计算各样本点Cook距离
which(ck>1) # 无大于1的值
#说明样本中没有强影响点

r = rstandard(lm2) # 计算各样本点SRE
p1 = ncol(data0)+1 # 自变量个数加一
n1 = nrow(data0) # 清点样本点个数
F = (n1 - p1 - 1) * r ^ 2 / (n1 - p1 -r ^ 2) # 计算异常点F检验统计量
p.value = 1 - pf(F, 1, n1 - p1 - 1) # 计算异常点F检验p值
del=which(p.value < 0.05) 
del # 在0.05显著性水平下,找出异常点
data0=data0[-del,] # 删去异常值点
lm3=lm(y ~ ., data=data0)#重新拟合模型
summary(lm3)

#-----检验异方差性detect heteroscedasticity-----#

abse = abs(resid(lm3)) # 计算残差绝对值
for(i in 1:(ncol(data0)-1))
{out=cor.test(data0[,i], abse, method = "spearman",exact=FALSE) # 计算残差与xi的相关系数
if(out$p.value<0.05)
  print(names(data0)[i])
}
# 结果表明残差绝对值与自变量存在显著相关关系,即数据存在异方差性
ymin=min(data0$y)
a=-2*ymin
data0$y=data0$y+a
bc = boxcox(y ~ ., data = data0, lambda = seq(-2, 2, 0.01))
# 计算不同lambda值对应BoxCox变换的似然函数，lambda取值区间为[-2, 2]，步长为0.01
# 输出的图像展示对数似然函数随lambda增长的变化
lambda = bc $ x[which.max(bc $ y)] # 选取使似然函数达到最大值的lambda值
lambda # 输出lambda值
y_bc = (data0 $ y ^ lambda - 1) / lambda # 计算变换后的y值

lm1_bc = lm(y_bc ~ ., data = data0) # 使用变换后的y值拟合线性模型
summary(lm1_bc) # 输出拟合结果
abse_bc = abs(resid(lm1_bc)) # 计算残差绝对值
for(i in 1:(ncol(data0)-1))
{out=cor.test(data0[,i], abse_bc, method = "spearman",exact=FALSE) # 计算残差与xi的相关系数
if(out$p.value<0.05)
  print(names(data0)[i])
}
# 借助boxcox变换，数据异方差性并未消除。
data0$y=data0$y-a# 把data0中的y值复原

s = seq(-0.5, 0.5, 0.01) # 产生数列(-0.49,-0.48,...,0.5)作为权重函数幂指数备选
result1 = vector(length = length(s), mode = "list") # 产生一个与s维度相同的空向量以储存后续结果
result2 = vector(length = length(s), mode = "list") # 产生一个与s维度相同的空向量以储存后续结果

p=(ncol(data0)-1) #自变量个数
corr2 = vector(length = p, mode = "list") # 产生一个长度为p的向量，用来存储rho的绝对值
for(i in 1:p)
{corr2[i]=abs(cor.test(data0[,i], abse, method = "spearman",exact=FALSE)$estimate) # 把rho的绝对值存储在向量corr2中
}
xk=which.max(corr2) # corr2的最大的元素对应等级相关系数最大的自变量
xk # 找出等级相关系数最大的自变量是第四个，也就是rm

for(j in 1 : length(s))
{
  
  w = data0 $ rm ^ (-s[j]) # 计算权重
  lm1_wlse = lm(y ~ ., weights = w, data0) # 用WLSE拟合线性模型
  result1[[j]] = logLik(lm1_wlse) # 储存拟合模型对应的对数似然函数值
  result2[[j]] = summary(lm1_wlse) # 储存模型拟合结果
  
}
result1 # 输出对数似然函数值
s[which.max(result1)]
result2[which.max(result1)] # 输出对应对数函数最大值的模型 
# 此时Ra square=0.8908,增大了
# 说明对于这个数据集加权最小二乘估计的拟合效果好于普通最小二乘估计


###-----检验自相关性-----###
dwtest(lm1, alternative = "two.sided") # 进行DW检验
# 检验表明数据存在显著的一阶自相关性。
e=resid(lm1)
N = length(e) # 清点样本点个数
plot(e[1 : (N - 1)], e[2 : N]) # 画出e_{i-1}与e_i关系图
# 图中显示e_{i-1}与e_i存在明显正相关。
rhohat = ( sum(e[1 : (N - 1)] ^ 2) * sum(e[2 : N] ^ 2) ) ^ (- 1 / 2) * sum(e[1 : (N - 1)] * e[2 : N]) 
#估计一阶自相关系数
rhohat
# 一阶自相关系数为0.4886222，表明数据存在明显的一阶自相关性。

aft=data0
for(i in 1 : (p + 1))
{
  j=N
  while(j>1)
  {
    aft[j,i]=aft[j,i]-rhohat*aft[j-1,i];
    j=j-1;
  }
}
aft=aft[-1,]
 # 用迭代法处理数据,以消除自相关性
lm_iterate = lm(y ~ ., data = aft) # 用迭代后数据拟合线性模型
summary(lm_iterate) # 输出拟合结果
dwtest(lm_iterate, alternative = "two.sided") # 进行DW检验
# 检验表明迭代后数据的一阶自相关性在显著性水平0.01下不明显。

newe = resid(lm_iterate) # 计算迭代后模型残差
m = length(newe) # 清点迭代后样本点个数
plot(newe[1 : (m - 1)], newe[2 : m]) # 画出newe_{i-1}与newe_i关系图
# 图中显示newe_{i-1}与newe_i不存在明显关联。
rhohat = ( sum(newe[1 : (m - 1)] ^ 2) * sum(newe[2 : m] ^ 2) ) ^ (- 1 / 2) * sum(newe[1 : (m - 1)] * newe[2 : m]) 
#估计一阶自相关系数
rhohat
# 一阶自相关系数为-0.1167036，绝对值较原始数据大幅降低。

###-----判定线性模型中是否存在多重共线性-----###
lm3=lm(y ~ ., data = aft) # 用处理后的数据拟合线性模型
summary(lm3)
vif(lm3) # 计算方差扩大因子
# 方差扩大因子都小于10
XX = cor(aft[1:m, -11]) # 计算设计矩阵的相关矩阵
kappa(XX, exact = T) # 计算条件数
# 条件数为25.78646，小于100，表明模型中多重共线性程度很小









