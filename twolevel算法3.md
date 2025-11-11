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
\min_{G_A,P_{A2B},P_{A2C}} {\color{red}C_{grid}(G_A)}-\lambda_A(G_A+PV_A-T_{A2B}-T_{A2C})+[\mu_1,\mu_2] [T_{A2B};T_{A2C}]+\alpha T_{A2B}^2+\alpha T_{A2C}^2\\
s.t. G_A\in[0,D_{max}];\\
T_{A2B}\in[-PV_{B},PV_{A}];\\
T_{A2C}\in[-PV_{C},PV_{A}];\\
T_{A2B}+T_{A2C}\leq PV_A
$$
微网B:
$$
\min_{G_B,\bar{T}_{A2B},T_{B2C}} C_{grid}(G_B)-\lambda_B(G_B+PV_B+\bar{T}_{A2B}-\bar{T}_{B2C})+[\mu_1,\mu_3] [-\bar{T}_{A2B},T_{B2C}]+\alpha \bar{T}_{A2B}^2+\alpha T_{B2C}^2\\
s.t. G_B\geq 0;\\
\bar{T}_{A2B}\in[-PV_{B},PV_{A}];\\
T_{B2C}\in[-PV_{C},PV_{B}];\\ 
-\bar{T}_{A2B}+T_{B2C}\leq PV_B
$$
微网C:
$$
\min_{G_C,\bar{T}_{A2C},\bar{T}_{B2C}} C_{grid}(G_C)-\lambda_C(G_C+PV_C+\bar{T}_{A2C}+\bar{T}_{B2C})+[\mu_2,\mu_3] [-\bar{T}_{A2C},-\bar{T}_{B2C}]+\alpha \bar{T}_{A2C}^2+\alpha \bar{T}_{B2C}^2\\
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
T_{A2B}-\bar{T}_{A2B} \\
T_{A2C}-\bar{T}_{A2C} \\
T_{B2C}-\bar{T}_{B2C}
\end{bmatrix}
$$


#### 虚拟代理商更新电价：

微网A代理商：
$$
\lambda_A=\lambda_A+D_A-(G_A+PV_A-T_{A2B}-T_{A2C})
$$
微网B代理商：
$$
\lambda_B=\lambda_B+D_B-(G_B+PV_B+T_{A2B}-T_{B2C})\\
$$

微网C代理商：
$$
\lambda_C=\lambda_C+D_C-(G_C+PV_C+T_{A2C}+T_{B2C})\\
$$


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

C_{tran}(x)=0.5*x\\
$$

数值计算结果：
$$
[D_A,D_B,D_C]=[2.6049, 0.1486, 0.1486]\\
[G_A,G_B,G_C]=[0.1049, 0.1486, 0.1486]\\
[\lambda_A,\lambda_B,\lambda_C]= [3.0980 1.2972 1.2972]\\
P_{A2B}=-1.0000\\
P_{A2C}=-1.0000\\
P_{B2C}=0\\
\mu=[-2.5000
   -2.5000
   -2.0000]
$$


```
 			PowerDemand: [2.6049 0.1486 0.1486]
            PowerSupply: [2.6049 0.1486 0.1486]
                 Lambda: [3.0980 1.2972 1.2972]
    number_of_Microgrid: 3
                     MA: [1×1 Micro_Community_A]
                     MB: [1×1 Micro_Community_B]
                     MC: [1×1 Micro_Community_C]
              TotalCost: -16.3544
              PowerGrid: [0.1049 0.1486 0.1486]
                     mu: [3×1 double]
           Powertrana2b: [-1.0000 -1.0000]
           Powertrana2c: [-1.0000 -1.0000]
           Powertranb2c: [2.3373e-08 -2.3383e-08]

```

