#### 微网A中用户要考虑的问题：用电效益最大+成本最小：

$$
\min_{D_A} -U_A(D_A)+\lambda_A D_A
$$

这里$D_A$指A微网中用户的用电功率/需求功率；$U_A(D_A)$为效益/效用函数；用电价格$\lambda_A$乘$D_A$表示用电成本。

当然可以考虑微网内的馈线约束来限制$D_A$​,可以考虑微网A中有多个用户(用户可以是EV电车，可以是单个居民，也可是商场，医院)。目前只实现简单约束(设置一个用电上界$D_{A,\max}$)：
$$
0\leq D_A\leq D_{A,\max}
$$
其他微网用户同上。



### 微网A的虚拟代理商需要考虑的问题：购电成本最小+售电收益最大+损耗最小：

$$
\min_{G_A,P_{A2B},P_{A2C}} {\color{red}C_{grid}(G_A)}-{\color{blue}\lambda_A(G_A+PV_A-T_{A2B}-T_{A2C}-B_A)}-[\mu_1,\mu_2] [T_{A2B};T_{A2C}]+\alpha T_{A2B}^2+\alpha T_{A2C}^2\\
$$

购电成本包括：从公用电网购电的成本为$\color{red}C_{grid}(G_A)$,其中$G_A$表示相应的购电功率；从B微网购电的成本:$\color{red}\mu_1T_{A2B}(T_{A2B}<0)$;以及从C微网购电的成本$\color{red}\mu_2T_{A2C}(T_{A2C}<0)$.需要指出的是$\mu_1$表示AB微网之间的买卖电价（是非负的），$T_{A2B}$表示微网A传输给微网B的功率（可正可负，$T_{A2B}$为负值时表示B传输给A，因此当$T_{A2B}$为负值时$\mu_1T_{A2B}$表示成本，反之表示收益）

售电收益包括：微网A的虚拟代理商售电给A微网用户的收益${\color{blue}\lambda_A(G_A+PV_A-T_{A2B}-T_{A2C}-B_A)}$,其中$\lambda_A$表示用户购电电价，A微网用户的总购电功率为代理商从公用电网购电功率$G_A$+光伏发电功率$PV_A$-传递给B微网的功率$T_{A2B}$-传递给C微网的功率$T_{A2C}$-电池充电功率$B_A$。另一部分收益是传输给B,C微网的收益：$\color{blue}\mu_1T_{A2B}(T_{A2B}>0)$;$\color{blue}\mu_2T_{A2C}(T_{A2C}>0)$.

微网之间买卖电的损耗与功率的平方成正比:$\alpha T_{A2B}^2+\alpha T_{A2C}^2$.这里$\alpha$为可调参数。添加这一项从算法优化的角度来讲，是为了保证目标函数严格凸，这是对偶上升方法快速收敛的前提。

决策变量为$G_A,P_{A2B},P_{A2C}$,其约束条件相对简单：从公共电网购电量不大于最大功率需求，微网间的传递功率之和不应高于当前微网光伏发电量(如果出现不合理的解可以尝试进一步约束决策变量，但目前来看，已经够了)

其他微网虚拟代理商同上；

### AB微网之间买卖电价调整策略：购买功率大于出售功率则价格上调；购买功率小于出售功率则价格下调

$$
\mu_1=\mu_1+\bar{T}_{A2B}-T_{A2B}
$$

$\bar{T}_{A2B}$表示B微网计划从A微网购买的功率即购买功率，$T_{A2B}$表示A微网计划售卖B微网的功率即出售功率，因此电价调整测量可以设置为$\mu_1=\mu_1+\gamma(\bar{T}_{A2B}-T_{A2B})$,$\gamma=1$为调整幅度。



## 虚拟代理商更新电价的策略：需求大则电价涨，需求少则电价跌。

因此微网A代理商更新电价策略为：
$$
\lambda_A=\lambda_A+D_A-(G_A+PV_A-T_{A2B}-T_{A2C}-B_A)
$$
其中$D_A$为微网A总需求功率，$G_A+PV_A-T_{A2B}-T_{A2C}-B_A$为微网A虚拟代理商总供给功率。



将上述四个过程整合可得如下算法：

#### 下层更新需求(Demand)功率：

微网A:
$$
\min_{D_A} -U_A(D_A)+\lambda_A D_A
$$
微网B:
$$
\min_{D_B} -U_A(D_B)+\lambda_B D_B
$$

微网C:
$$
\min_{D_C} -U_C(D_C)+\lambda_C D_C
$$

#### 上层更新供给(Supply)功率：

微网A:
$$
\min_{G_A,P_{A2B},P_{A2C}} C_{grid}(G_A)-\lambda_A(G_A+PV_A-T_{A2B}-T_{A2C}-B_A)-\mu_1T_{A2B}-\mu_2T_{A2C}+\alpha T_{A2B}^2+\alpha T_{A2C}^2\\
s.t. G_A\in[0,D_{max}];\\
T_{A2B}\in[-PV_{B},PV_{A}];\\
T_{A2C}\in[-PV_{C},PV_{A}];\\
T_{A2B}+T_{A2C}\leq PV_A
$$
微网B:
$$
\min_{G_B,\bar{T}_{A2B},T_{B2C}} C_{grid}(G_B)-\lambda_B(G_B+PV_B+\bar{T}_{A2B}-\bar{T}_{B2C}-B_B)+\mu_1\bar{T}_{A2B}-\mu_3T_{B2C}+\alpha \bar{T}_{A2B}^2+\alpha T_{B2C}^2\\
s.t. G_B\geq 0;\\
\bar{T}_{A2B}\in[-PV_{B},PV_{A}];\\
T_{B2C}\in[-PV_{C},PV_{B}];\\ 
-\bar{T}_{A2B}+T_{B2C}\leq PV_B
$$
微网C:
$$
\min_{G_C,\bar{T}_{A2C},\bar{T}_{B2C}} C_{grid}(G_C)-\lambda_C(G_C+PV_C+\bar{T}_{A2C}+\bar{T}_{B2C}-B_C)+\mu_2\bar{T}_{A2C}+\mu_3\bar{T}_{B2C}+\alpha \bar{T}_{A2C}^2+\alpha \bar{T}_{B2C}^2\\
s.t. G_C\geq 0;\\
\bar{T}_{A2C}\in[-PV_{C},PV_{A}];\\ 
\bar{T}_{B2C}\in[-PV_{C},PV_{B}];\\ 
-\bar{T}_{A2C}-\bar{T}_{B2C}\leq PV_C
$$
更新对偶变量$\mu$:
$$
\begin{bmatrix}
\mu_1 \\
\mu_2 \\
\mu_3
\end{bmatrix}=\begin{bmatrix}
\mu_1 \\
\mu_2 \\
\mu_3
\end{bmatrix}
+\begin{bmatrix}
\bar{T}_{A2B}-T_{A2B} \\
\bar{T}_{A2C}-T_{A2C} \\
\bar{T}_{B2C}-T_{B2C}
\end{bmatrix}
$$


#### 虚拟代理商更新电价：

微网A代理商：
$$
\lambda_A=\lambda_A+D_A-(G_A+PV_A-T_{A2B}-T_{A2C}-B_A)
$$
微网B代理商：
$$
\lambda_B=\lambda_B+D_B-(G_B+PV_B+T_{A2B}-T_{B2C}-B_B)\\
$$

微网C代理商：
$$
\lambda_C=\lambda_C+D_C-(G_C+PV_C+T_{A2C}+T_{B2C}-B_C)\\
$$

也可以从另一个角度(全局角度)思考问题：将上述优化问题加在一起得
$$
\min_{D_A,D_B,D_C,G_A,G_B,G_C,T_{A2B},T_{A2C},T_{B2C}} \sum_{\beta\in\{A,B,C\}}\Big(-U_\beta(D_\beta) +C_{grid}(G_\beta)\Big) +2\alpha {T}_{B2C}^2+2\alpha {T}_{A2B}^2+2\alpha {T}_{A2C}^2\\ 
s.t.
D_A=(G_A+PV_A-T_{A2B}-T_{A2C}-B_A)\\
D_B=(G_B+PV_B+T_{A2B}-T_{B2C}-B_B)\\
D_C=(G_C+PV_C+T_{A2C}+T_{B2C}-B_C)\\
$$
注意到电池充电功率$B_A,B_B,B_C$在上述算法中是已经确定的。事实上我们通过一个简单的判断函数来确定电池充电功率$B_A,B_B,B_C$:
第一步：令$B_A=B_B=B_C=0$，并执行上述算法得到电价$\lambda_A,\lambda_B,\lambda_C$;

第二步：$B_A=f(\lambda_A),B_B=f(\lambda_B),B_C=f(\lambda_C)$,$f$为分段函数，当电价高于阈值时返回放电功率$P_{dis}$，低于设定阈值时返回充电功率$P_{char}$.

第三步：重新执行上述算法



## 数值实验参数：

$$
U_A(D_A)=\min(10*\sqrt{D_A},100)\\
C_{grid}(G_A)=10*G_A^2+G_A\\
PV_A=0.5;P_{A,Max}=40;\\
\\
U_B(D_B)=\min(1*\sqrt{D_B},0.5)\\
C_{grid}(G_B)=1*G_B^2+G_B\\ 
PV_A=1;P_{A,Max}=40;\\
\\U_C(D_C)=\min(1*\sqrt{D_C},0.5)\\
C_{grid}(G_C)=1*G_C^2+G_C\\ 
PV_C=1;P_{C,Max}=40;\\
\\
$$

