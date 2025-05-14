import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
# 定义常量
fs = 290e3;
Lr = 12e-6;
Cr = 48e-9;
Z0 = np.sqrt(Lr/Cr);  # Use np.sqrt instead of sqrt
fr = 1/(2*np.pi*np.sqrt(Lr*Cr));  # Use np.sqrt instead of sqrt
td = 50e-9;
Dd = td*fs;
r = fr/fs;
omega_r = 2*np.pi*fr;
Ln = 100;
Vin = 800;
Vo = 1000;
ILm0 = np.pi*r*Vo/(2*Ln*Z0);5
Ioss = 0.06; 



# 定义 g(Dy1, Dy2, Dp, fn=1.2) 函数（用户提供的公式）
def g_function(vars):
    Dy1, Dy2, Dp = vars

    # 设定 Vin, Vo, Z0, r 的值（假设为 1，用户需提供具体值）
    Vin = 800
    Vo = 660
    Z0 = np.sqrt(12e-6/48e-9)
    r = 1 /1.2

    # 计算 g 的值
    g = Z0/Vin*((1/2/np.pi/r)*(-(2*Vin**2*np.sin(Dy1*np.pi*r) + 2*Vo**2*np.sin(Dy2*np.pi*r) + 2*Vin**2*np.sin(np.pi*r*(Dy1 - 1)) + 2*Vo**2*np.sin(np.pi*r*(Dy2 - 1)) + 2*Vin**2*np.sin(np.pi*r) + 2*Vo**2*np.sin(np.pi*r) + 2*Vin*Vo*np.sin((np.pi*r*(2*Dp - Dy1 - Dy2 + 2))/2) - 2*Vin*Vo*np.sin((np.pi*r*(Dy1 - 2*Dp + Dy2))/2) + 2*Vin*Vo*np.sin((np.pi*r*(2*Dp + Dy1 + Dy2 - 2))/2) + 2*Vin*Vo*np.sin((np.pi*r*(2*Dp + Dy1 + Dy2 - 4))/2) + 2*Vin*Vo*np.sin((np.pi*r*(2*Dp + Dy1 - Dy2))/2) + 2*Vin*Vo*np.sin((np.pi*r*(2*Dp - Dy1 + Dy2))/2) + 2*Vin*Vo*np.sin((np.pi*r*(2*Dp + Dy1 - Dy2 - 2))/2) + 2*Vin*Vo*np.sin((np.pi*r*(2*Dp - Dy1 + Dy2 - 2))/2) - 2*Vin**2*np.pi*r - 2*Vo**2*np.pi*r + 2*Vin**2*np.pi*r*np.cos(Dy1*np.pi*r) + 2*Vo**2*np.pi*r*np.cos(Dy2*np.pi*r) - 2*Vin*Vo*np.pi*r*np.cos((np.pi*r*(Dy1 - 2*Dp + Dy2))/2) + 4*Vin*Vo*np.pi*r*np.cos((np.pi*r*(2*Dp + Dy1 + Dy2 - 2))/2) + 2*Vin*Vo*np.pi*r*np.cos((np.pi*r*(2*Dp + Dy1 + Dy2 - 4))/2) - 2*Dy1*Vin**2*np.pi*r*np.cos(Dy1*np.pi*r) - 2*Dy2*Vo**2*np.pi*r*np.cos(Dy2*np.pi*r) + 2*Vin*Vo*np.pi*r*np.cos((np.pi*r*(2*Dp + Dy1 - Dy2))/2) + 2*Vin*Vo*np.pi*r*np.cos((np.pi*r*(2*Dp - Dy1 + Dy2))/2) - 2*Dy1*Vin**2*np.pi*r*np.cos(np.pi*r*(Dy1 - 1)) - 2*Dy2*Vo**2*np.pi*r*np.cos(np.pi*r*(Dy2 - 1)) - 2*Dp*Vin*Vo*np.pi*r*np.cos((np.pi*r*(2*Dp - Dy1 - Dy2 + 2))/2) + Dy1*Vin*Vo*np.pi*r*np.cos((np.pi*r*(2*Dp - Dy1 - Dy2 + 2))/2) + Dy2*Vin*Vo*np.pi*r*np.cos((np.pi*r*(2*Dp - Dy1 - Dy2 + 2))/2) - 2*Dp*Vin*Vo*np.pi*r*np.cos((np.pi*r*(Dy1 - 2*Dp + Dy2))/2) - 2*Dp*Vin*Vo*np.pi*r*np.cos((np.pi*r*(2*Dp + Dy1 + Dy2 - 2))/2) - 2*Dp*Vin*Vo*np.pi*r*np.cos((np.pi*r*(2*Dp + Dy1 + Dy2 - 4))/2) + Dy1*Vin*Vo*np.pi*r*np.cos((np.pi*r*(Dy1 - 2*Dp + Dy2))/2) + Dy2*Vin*Vo*np.pi*r*np.cos((np.pi*r*(Dy1 - 2*Dp + Dy2))/2) - Dy1*Vin*Vo*np.pi*r*np.cos((np.pi*r*(2*Dp + Dy1 + Dy2 - 2))/2) - Dy2*Vin*Vo*np.pi*r*np.cos((np.pi*r*(2*Dp + Dy1 + Dy2 - 2))/2) - Dy1*Vin*Vo*np.pi*r*np.cos((np.pi*r*(2*Dp + Dy1 + Dy2 - 4))/2) - Dy2*Vin*Vo*np.pi*r*np.cos((np.pi*r*(2*Dp + Dy1 + Dy2 - 4))/2) - 2*Dp*Vin*Vo*np.pi*r*np.cos((np.pi*r*(2*Dp + Dy1 - Dy2))/2) - 2*Dp*Vin*Vo*np.pi*r*np.cos((np.pi*r*(2*Dp - Dy1 + Dy2))/2) - 2*Dp*Vin*Vo*np.pi*r*np.cos((np.pi*r*(2*Dp + Dy1 - Dy2 - 2))/2) - 2*Dp*Vin*Vo*np.pi*r*np.cos((np.pi*r*(2*Dp - Dy1 + Dy2 - 2))/2) - Dy1*Vin*Vo*np.pi*r*np.cos((np.pi*r*(2*Dp + Dy1 - Dy2))/2) + Dy1*Vin*Vo*np.pi*r*np.cos((np.pi*r*(2*Dp - Dy1 + Dy2))/2) + Dy2*Vin*Vo*np.pi*r*np.cos((np.pi*r*(2*Dp + Dy1 - Dy2))/2) - Dy2*Vin*Vo*np.pi*r*np.cos((np.pi*r*(2*Dp - Dy1 + Dy2))/2) - Dy1*Vin*Vo*np.pi*r*np.cos((np.pi*r*(2*Dp + Dy1 - Dy2 - 2))/2) + Dy1*Vin*Vo*np.pi*r*np.cos((np.pi*r*(2*Dp - Dy1 + Dy2 - 2))/2) + Dy2*Vin*Vo*np.pi*r*np.cos((np.pi*r*(2*Dp + Dy1 - Dy2 - 2))/2) - Dy2*Vin*Vo*np.pi*r*np.cos((np.pi*r*(2*Dp - Dy1 + Dy2 - 2))/2))/(4*Z0**2*(np.cos(np.pi*r) + 1))));
    return g

def f_function(vars):
    Dy1, Dy2, Dp = vars

    # 设定 Vin, Vo, Z0, r 的值（假设为 1，用户需提供具体值）
    Vin = 800
    Vo = 660
    Z0 = np.sqrt(12e-6/48e-9)
    r = 1 /1.2

    # 计算 g 的值
    f = ((np.cos(np.pi*r*(Dp-0.5))*np.cos(np.pi*r*(1-Dy1)/2)*np.cos(np.pi*r*(1-Dy2)/2))/(np.pi*r*np.cos(np.pi*r/2))-1/(np.pi*r));
    return f

# 参考注释，别删了
# IL0_RegionE = -(((np.cos(np.pi*r) + 1)*(Vo*np.sin((r*np.pi*(2*Dp - Dy1 + Dy2 - 2))/2) + Vo*np.sin(np.pi*r) + Vin*np.sin(np.pi*Dy1*r) - Vo*np.sin((r*np.pi*(Dy1 - 2*Dp + Dy2))/2)))/Z0 + (np.sin(np.pi*r)*(Vin - Vo + Vo*np.cos((r*np.pi*(2*Dp - Dy1 + Dy2 - 2))/2) - Vo*np.cos(np.pi*r) - Vin*np.cos(np.pi*Dy1*r) + Vo*np.cos((r*np.pi*(Dy1 - 2*Dp + Dy2))/2)))/Z0)/(2*(np.cos(np.pi*r) + 1));
# IL1_RegionE = -(Vo*np.sin((r*np.pi*(2*Dp + Dy1 - Dy2 - 2))/2) + Vin*np.sin(np.pi*r) + Vo*np.sin((r*np.pi*(2*Dp + Dy1 + Dy2 - 4))/2) + Vin*np.sin(np.pi*r*(Dy1 - 1)) + Vo*np.cos((r*np.pi*(2*Dp + Dy1 + Dy2 - 4))/2)*np.sin(np.pi*r) + Vo*np.sin((r*np.pi*(2*Dp + Dy1 + Dy2 - 4))/2)*np.cos(np.pi*r) + Vin*np.cos(np.pi*r*(Dy1 - 1))*np.sin(np.pi*r) + Vin*np.sin(np.pi*r*(Dy1 - 1))*np.cos(np.pi*r) + Vo*np.cos((r*np.pi*(2*Dp + Dy1 - Dy2 - 2))/2)*np.sin(np.pi*r) + Vo*np.sin((r*np.pi*(2*Dp + Dy1 - Dy2 - 2))/2)*np.cos(np.pi*r))/(2*Z0*(np.cos(np.pi*r) + 1));
# IL2_RegionE = (((np.cos(np.pi*r) + 1)*(Vin*np.sin((r*np.pi*(2*Dp - Dy1 + Dy2))/2) - Vin*np.sin(np.pi*r) + Vo*np.sin(np.pi*Dy2*r) + Vin*np.sin((r*np.pi*(2*Dp + Dy1 + Dy2 - 2))/2)))/Z0 + (np.sin(np.pi*r)*(Vin + Vo - Vin*np.cos((r*np.pi*(2*Dp - Dy1 + Dy2))/2) + Vin*np.cos(np.pi*r) - Vo*np.cos(np.pi*Dy2*r) - Vin*np.cos((r*np.pi*(2*Dp + Dy1 + Dy2 - 2))/2)))/Z0)/(2*(np.cos(np.pi*r) + 1));
# IL3_RegionE = (Vin*np.sin((r*np.pi*(2*Dp + Dy1 - Dy2))/2) + Vo*np.sin(np.pi*r) + Vin*np.sin((r*np.pi*(2*Dp - Dy1 - Dy2 + 2))/2) + Vo*np.sin(np.pi*r*(Dy2 - 1)) + Vin*np.cos(np.pi*r)*np.sin((r*np.pi*(2*Dp - Dy1 - Dy2 + 2))/2) - Vin*np.sin(np.pi*r)*np.cos((r*np.pi*(2*Dp - Dy1 - Dy2 + 2))/2) + Vo*np.cos(np.pi*r*(Dy2 - 1))*np.sin(np.pi*r) + Vo*np.sin(np.pi*r*(Dy2 - 1))*np.cos(np.pi*r) - Vin*np.cos((r*np.pi*(2*Dp + Dy1 - Dy2))/2)*np.sin(np.pi*r) + Vin*np.sin((r*np.pi*(2*Dp + Dy1 - Dy2))/2)*np.cos(np.pi*r))/(2*Z0*(np.cos(np.pi*r) + 1));


# 定义约束函数（ZVS 约束）
def constraint1(vars):
    Dy1, Dy2, Dp = vars
    I_zvs_p12_RegionE = -(Vo*np.sin((r*np.pi*(2*Dp + Dy1 - Dy2 - 2))/2) + Vin*np.sin(np.pi*r) + Vo*np.sin((r*np.pi*(2*Dp + Dy1 + Dy2 - 4))/2) + Vin*np.sin(np.pi*r*(Dy1 - 1)) + Vo*np.cos((r*np.pi*(2*Dp + Dy1 + Dy2 - 4))/2)*np.sin(np.pi*r) + Vo*np.sin((r*np.pi*(2*Dp + Dy1 + Dy2 - 4))/2)*np.cos(np.pi*r) + Vin*np.cos(np.pi*r*(Dy1 - 1))*np.sin(np.pi*r) + Vin*np.sin(np.pi*r*(Dy1 - 1))*np.cos(np.pi*r) + Vo*np.cos((r*np.pi*(2*Dp + Dy1 - Dy2 - 2))/2)*np.sin(np.pi*r) + Vo*np.sin((r*np.pi*(2*Dp + Dy1 - Dy2 - 2))/2)*np.cos(np.pi*r))/(2*Z0*(np.cos(np.pi*r) + 1))
    return I_zvs_p12_RegionE # I_zvs_p12_RegionE should be less than 0

def constraint2(vars):
    Dy1, Dy2, Dp = vars
    I_zvs_s12_RegionE = (Vin*np.sin((r*np.pi*(2*Dp + Dy1 - Dy2))/2) + Vo*np.sin(np.pi*r) + Vin*np.sin((r*np.pi*(2*Dp - Dy1 - Dy2 + 2))/2) + Vo*np.sin(np.pi*r*(Dy2 - 1)) + Vin*np.cos(np.pi*r)*np.sin((r*np.pi*(2*Dp - Dy1 - Dy2 + 2))/2) - Vin*np.sin(np.pi*r)*np.cos((r*np.pi*(2*Dp - Dy1 - Dy2 + 2))/2) + Vo*np.cos(np.pi*r*(Dy2 - 1))*np.sin(np.pi*r) + Vo*np.sin(np.pi*r*(Dy2 - 1))*np.cos(np.pi*r) - Vin*np.cos((r*np.pi*(2*Dp + Dy1 - Dy2))/2)*np.sin(np.pi*r) + Vin*np.sin((r*np.pi*(2*Dp + Dy1 - Dy2))/2)*np.cos(np.pi*r))/(2*Z0*(np.cos(np.pi*r) + 1))
    return I_zvs_s12_RegionE - 0  # I_zvs_s12_RegionE should be greater than 0

def constraint3(vars):
    Dy1, Dy2, Dp = vars
    I_zvs_p34_RegionE = (((np.cos(np.pi*r) + 1)*(Vo*np.sin((r*np.pi*(2*Dp - Dy1 + Dy2 - 2))/2) + Vo*np.sin(np.pi*r) + Vin*np.sin(np.pi*Dy1*r) - Vo*np.sin((r*np.pi*(Dy1 - 2*Dp + Dy2))/2)))/Z0 + (np.sin(np.pi*r)*(Vin - Vo + Vo*np.cos((r*np.pi*(2*Dp - Dy1 + Dy2 - 2))/2) - Vo*np.cos(np.pi*r) - Vin*np.cos(np.pi*Dy1*r) + Vo*np.cos((r*np.pi*(Dy1 - 2*Dp + Dy2))/2)))/Z0)/(2*(np.cos(np.pi*r) + 1))
    return I_zvs_p34_RegionE - 0  # I_zvs_p34_RegionE should be greater than 0

def constraint4(vars):
    Dy1, Dy2, Dp = vars
    I_zvs_s34_RegionE = -(((np.cos(np.pi*r) + 1)*(Vin*np.sin((r*np.pi*(2*Dp - Dy1 + Dy2))/2) - Vin*np.sin(np.pi*r) + Vo*np.sin(np.pi*Dy2*r) + Vin*np.sin((r*np.pi*(2*Dp + Dy1 + Dy2 - 2))/2)))/Z0 + (np.sin(np.pi*r)*(Vin + Vo - Vin*np.cos((r*np.pi*(2*Dp - Dy1 + Dy2))/2) + Vin*np.cos(np.pi*r) - Vo*np.cos(np.pi*Dy2*r) - Vin*np.cos((r*np.pi*(2*Dp + Dy1 + Dy2 - 2))/2)))/Z0)/(2*(np.cos(np.pi*r) + 1))
    return I_zvs_s34_RegionE  # I_zvs_s34_RegionE should be less than 0

# 新增约束函数，确保 PoN 在 0.2 到 1.2 的范围内
def constraint5(vars, PoN_target):
    PoN = f_function(vars)
    return PoN - PoN_target  # PoN should be equal to PoN_target

# 进行优化
# 变量边界
bounds = [(0,1), (0,1), (0,0.5)]
initial_guess = [0.5, 0.5, 0.25]

# 存储优化结果
optimal_points = []

# 对每一个 PoN 进行优化求解
for PoN_target in np.linspace(0.2, 1.2, 80):
    constraints = [{'type': 'ineq', 'fun': constraint1},
                   {'type': 'ineq', 'fun': constraint2},
                   {'type': 'ineq', 'fun': constraint3},
                   {'type': 'ineq', 'fun': constraint4},
                   {'type': 'eq', 'fun': lambda vars: constraint5(vars, PoN_target)}]
    result = minimize(g_function, initial_guess, bounds=bounds, constraints=constraints)
    optimal_points.append(result.x)

# 可视化最优点
fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111, projection='3d')

# 绘制优化后的点
for point in optimal_points:
    ax.scatter(point[0], point[1], point[2], color='red', marker='*', s=10)

# 将点连起来，不要首尾相连
optimal_points = np.array(optimal_points)
ax.plot(optimal_points[:, 0], optimal_points[:, 1], optimal_points[:, 2], color='blue')

# 轴标签
ax.set_xlim([0, 1])
ax.set_ylim([0, 1])
ax.set_zlim([0, 0.5])
ax.set_xlabel('Dy1')
ax.set_ylabel('Dy2')
ax.set_zlabel('Dp')
ax.set_title('Optimized g(Dy1, Dy2, Dp) for different PoN values')

plt.show()

# 返回优化结果
optimal_points
