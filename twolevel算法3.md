#### 下层更新需求(Demand)功率：

$$
U(x,x_0)=\begin{cases}
\displaystyle
\beta\,\ln\!\big(1+\sigma\,x\big),  0 \le x \le x_0, \\[10pt]
\displaystyle
\beta\,\ln\!\big(1+\sigma\,x_0\big),   x \ge x_0.
\end{cases}
$$

微网A[居民(m)+EV(n)]:
$$
\min_{x^A_{i,t}}  \sum_i^{m+n}\left[-U(x^A_{i,t},D^A_{i,t})+\lambda^A x^A_{i,t}\right]\\
$$
微网B[CBD+EV(n)]:
$$
\min_{x^B_{j,t}}  \sum_j^{1+n}\left[-U(x^B_{j,t},D^B_{j,t})+\lambda^B x^B_{j,t}\right]
$$

微网C(Hospital):
$$
\min_{x^C_{t}}   \left[-U(x^C_{t},D^C_{t})+\lambda^C x^C_{t}\right]
$$

#### 上层更新供给(Supply)功率：

$$
C_{grid}(x)= a_0x^2+b_0x+c_0;\\ 
C_{EBSS}(x,SOC)=   (1-SOC)/SOC*exp(3)*0.1*(exp(x)-1).*(x <= 0) + (0.01*x.^2 + (1-SOC)/SOC*exp(3)*0.1*x).*(x > 0); 
$$

微网A:
$$
\min_{G_A,P_{A2B},P_{A2C},B_A} 
& C_{grid}(G_A)-\lambda^A(G_A+PV_A-T_{A2B}-T_{A2C}-B_A)\\
&+C_{EBSS}(B_A,SOC)-\mu_1T_{A2B}-\mu_2T_{A2C}+\alpha T_{A2B}^2+\alpha T_{A2C}^2\\
s.t.& G_A\in[0,D_{max}]; 
T_{A2B}\in[-PV_{B},PV_{A}]; 
T_{A2C}\in[-PV_{C},PV_{A}]; 
T_{A2B}+T_{A2C}\leq PV_A\\
EBSS:&SOC=SOC+\frac{B_A}{E_A}
$$
微网B:
$$
\min_{G_B,\bar{T}_{A2B},T_{B2C},B_B}& C_{grid}(G_B)-\lambda_B(G_B+PV_B+\bar{T}_{A2B}-\bar{T}_{B2C}-B_B)\\
&+C_{EBSS}(B_B,SOC)+\mu_1\bar{T}_{A2B}-\mu_3T_{B2C}+\alpha \bar{T}_{A2B}^2+\alpha T_{B2C}^2\\
s.t.& G_B\geq 0; 
\bar{T}_{A2B}\in[-PV_{B},PV_{A}]; 
T_{B2C}\in[-PV_{C},PV_{B}]; 
-\bar{T}_{A2B}+T_{B2C}\leq PV_B\\
EBSS:&SOC=SOC+\frac{B_B}{E_B}
$$
微网C:
$$
\min_{G_C,\bar{T}_{A2C},\bar{T}_{B2C},B_C}& C_{grid}(G_C)-\lambda_C(G_C+PV_C+\bar{T}_{A2C}+\bar{T}_{B2C}-B_C)\\
&+C_{EBSS}(B_C,SOC)+\mu_2\bar{T}_{A2C}+\mu_3\bar{T}_{B2C}+\alpha \bar{T}_{A2C}^2+\alpha \bar{T}_{B2C}^2\\
s.t.& G_C\geq 0; 
\bar{T}_{A2C}\in[-PV_{C},PV_{A}]; 
\bar{T}_{B2C}\in[-PV_{C},PV_{B}]; 
-\bar{T}_{A2C}-\bar{T}_{B2C}\leq PV_C\\
EBSS:&SOC=SOC+\frac{B_C}{E_C}
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
\lambda_A=\lambda_A+\sum_i^{m+n} x^A_{i,t} -(G_A+PV_A-T_{A2B}-T_{A2C}-B_A)
$$
微网B代理商：
$$
\lambda_B=\lambda_B+\sum_j^{1+n} x^B_{j,t}-(G_B+PV_B+T_{A2B}-T_{B2C}-B_B)\\
$$

微网C代理商：
$$
\lambda_C=\lambda_C+ x^C_{i,t}-(G_C+PV_C+T_{A2C}+T_{B2C}-B_C)\\
$$

 



$$
\min_{x^A_{i,t},G_A,B_A}&  \sum_i^{n}\left[-U(x^A_{i,t},D^A_{i,t})\right] +C_{grid,EBSS}(G_A,B_A,SOC)\\
s.t.&\sum_i^{m+n} x^A_{i,t}= G_A+PV_A+B_A
$$

$$
x_i=argmin -U(x_i)+\lambda x_i +\frac{\rho}{2}\|\sum_{j\neq i}^{n} x_{j,t}^{old}- (G_A+PV_A+B_A)\|^2\\
(G_A,B_A)=argmin C_{grid,EBSS}(G_A,B_A,SOC)-\lambda (G_A+PV_A+B_A) +\frac{\rho}{2}\|\sum_i^{m+n} x^{old}_{i,t}- (G_A+PV_A+B_A)\|^2\\
\lambda=\lambda+\rho(\sum_{j\neq i}^{n} x_{j,t}^{old}- (G_A+PV_A+B_A))
$$


$$

$$

