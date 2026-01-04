$$
\min \sum_i^3(-U_i(x_i)+C(g_i)+E(B_i))+\alpha T_{12}^2+\alpha T_{13}^2+\alpha T_{23}^2\\
st.x_1=g_1+PV_1+B_1-T_{12}-T_{13}\\
   x_2=g_2+PV_2+B_2+T_{12}-T_{23}\\
   x_3=g_3+PV_3+B_3+T_{13}+T_{23}\\
$$

$$
\min  (-U_1(x_1)+C(g_1)+E(B_1))+\alpha/2 T_{12}^2+\alpha/2 T_{13}^2 \\
st.x_1=g_1+PV_1+B_1-T_{12}-T_{13}\\ 

\min  (-U_2(x_2)+C(g_2)+E(B_2))+\alpha/2 T_{12}^2+\alpha/2 T_{23}^2 \\
st.x_2=g_2+PV_2+B_2+T_{12}-T_{23}\\ 

\min  (-U_3(x_3)+C(g_3)+E(B_3))+\alpha/2 T_{13}^2+\alpha/2 T_{23}^2 \\
st.x_3=g_3+PV_3+B_3+T_{13}+T_{33}\\ 
$$


$$
L(x_i,g_i,B_i,T_{ij})=\sum_i^3(-U_i(x_i)+C(g_i)+E(B_i))+\alpha T_{12}^2+\alpha T_{13}^2+\alpha T_{23}^2\\
+\lambda_1(-x_1+g_1+PV_1+B_1-T_{12}-T_{13})+\frac{\rho}{2}\|-x_1+g_1+PV_1+B_1-T_{12}-T_{13}\|^2\\
+\lambda_2(-x_2+g_2+PV_2+B_2+T_{12}-T_{23})+\frac{\rho}{2}\|-x_2+g_2+PV_2+B_2+T_{12}-T_{23}\|^2\\
+\lambda_3(-x_3+g_3+PV_3+B_3+T_{13}+T_{23})+\frac{\rho}{2}\|-x_3+g_3+PV_3+B_3+T_{13}+T_{23}\|^2\\
$$

$$
x_i=argmin\quad\sum_i^3(-U_i(x_i)+C(g_i)+E(B_i))+\alpha T_{12}^2+\alpha T_{13}^2+\alpha T_{23}^2\\
+\lambda_1(-x_1+g_1+PV_1+B_1-T_{12}-T_{13})+\frac{\rho}{2}\|-x_1+g_1+PV_1+B_1-T_{12}-T_{13}\|^2\\
+\lambda_2(-x_2+g_2+PV_2+B_2+T_{12}-T_{23})+\frac{\rho}{2}\|-x_2+g_2+PV_2+B_2+T_{12}-T_{23}\|^2\\
+\lambda_3(-x_3+g_3+PV_3+B_3+T_{13}+T_{23})+\frac{\rho}{2}\|-x_3+g_3+PV_3+B_3+T_{13}+T_{23}\|^2\\
$$


$$
\min -U(-T_{12}-T_{13})+\mu_{12}T_{12}+\mu_{13}T_{13} \\
\min -U(-T_{12}-T_{13})+\alpha T_{12}^2+\alpha T_{13}^2+\alpha T_{23}^2\\ 
\min -U(-T_{12}-T_{13})+\alpha T_{12}^2+\alpha T_{13}^2+\alpha T_{23}^2\\ 
$$