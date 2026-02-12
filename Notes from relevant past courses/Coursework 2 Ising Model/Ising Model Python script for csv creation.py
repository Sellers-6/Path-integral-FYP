# Making csv for MFT 3D

import pandas as pd
import numpy as np

T_values = np.linspace(2.5, 5.5, 301)

def m(T):
    return (-(T - 4.51) / 4.51) ** 0.5
m_values = m(T_values)

def m2(T):
    return (-(T - 4.51) / 4.51) ** 0.10
m2_values = m2(T_values)

def chi(T):
    return (abs(T - 4.51) / 4.51) ** - 1
chi_values = chi(T_values)

data = pd.DataFrame({"Temperature": T_values, "avg_mag" : m_values, "avg_mag2" : m2_values, "var_mag" : chi_values})
data.to_csv("function_values.csv", index=False)

# Concatenating the csvs from different seeds (1)

data1 = pd.read_csv("C:\\Users\\gs010\\Desktop\\Bath\\Maths and Physics\\Year 3 CW\\Computational Physics 1B - PH30056\\CW2\\Ising model - Copy\\IMQ1.csv")
data2 = pd.read_csv("C:\\Users\\gs010\\Desktop\\Bath\\Maths and Physics\\Year 3 CW\\Computational Physics 1B - PH30056\\CW2\\Ising model - Copy (2)\\IMQ1.csv")
data3 = pd.read_csv("C:\\Users\\gs010\\Desktop\\Bath\\Maths and Physics\\Year 3 CW\\Computational Physics 1B - PH30056\\CW2\\Ising model - Copy (3)\\IMQ1.csv")
data4 = pd.read_csv("C:\\Users\\gs010\\Desktop\\Bath\\Maths and Physics\\Year 3 CW\\Computational Physics 1B - PH30056\\CW2\\Ising model - Copy (4)\\IMQ1.csv")
data5 = pd.read_csv("C:\\Users\\gs010\\Desktop\\Bath\\Maths and Physics\\Year 3 CW\\Computational Physics 1B - PH30056\\CW2\\Ising model - Copy (5)\\IMQ1.csv")
data6 = pd.read_csv("C:\\Users\\gs010\\Desktop\\Bath\\Maths and Physics\\Year 3 CW\\Computational Physics 1B - PH30056\\CW2\\Ising model - Copy (6)\\IMQ1.csv")

data = pd.concat([data1, data2, data3, data4, data5, data6])

data.to_csv('IMQ1c.csv', index=False)

# Concatenating the csvs from different seeds (3)

data1 = pd.read_csv("C:\\Users\\gs010\\Desktop\\Bath\\Maths and Physics\\Year 3 CW\\Computational Physics 1B - PH30056\\CW2\\Ising model - Copy\\IMQ3.csv")
data2 = pd.read_csv("C:\\Users\\gs010\\Desktop\\Bath\\Maths and Physics\\Year 3 CW\\Computational Physics 1B - PH30056\\CW2\\Ising model - Copy (2)\\IMQ3.csv")
data3 = pd.read_csv("C:\\Users\\gs010\\Desktop\\Bath\\Maths and Physics\\Year 3 CW\\Computational Physics 1B - PH30056\\CW2\\Ising model - Copy (3)\\IMQ3.csv")
data4 = pd.read_csv("C:\\Users\\gs010\\Desktop\\Bath\\Maths and Physics\\Year 3 CW\\Computational Physics 1B - PH30056\\CW2\\Ising model - Copy (4)\\IMQ3.csv")
data5 = pd.read_csv("C:\\Users\\gs010\\Desktop\\Bath\\Maths and Physics\\Year 3 CW\\Computational Physics 1B - PH30056\\CW2\\Ising model - Copy (5)\\IMQ3.csv")
data6 = pd.read_csv("C:\\Users\\gs010\\Desktop\\Bath\\Maths and Physics\\Year 3 CW\\Computational Physics 1B - PH30056\\CW2\\Ising model - Copy (6)\\IMQ3.csv")

data = pd.concat([data1, data2, data3, data4, data5, data6])

data.to_csv('IMQ3c.csv', index=False)

# Concatenating the csvs from different seeds (3x)

datax01 = pd.read_csv("C:\\Users\\gs010\\Desktop\\Bath\\Maths and Physics\\Year 3 CW\\Computational Physics 1B - PH30056\\CW2\\Ising model - Copy\\IMQ3.csv")
datax02 = pd.read_csv("C:\\Users\\gs010\\Desktop\\Bath\\Maths and Physics\\Year 3 CW\\Computational Physics 1B - PH30056\\CW2\\Ising model - Copy (2)\\IMQ3.csv")
datax03 = pd.read_csv("C:\\Users\\gs010\\Desktop\\Bath\\Maths and Physics\\Year 3 CW\\Computational Physics 1B - PH30056\\CW2\\Ising model - Copy (3)\\IMQ3.csv")
datax04 = pd.read_csv("C:\\Users\\gs010\\Desktop\\Bath\\Maths and Physics\\Year 3 CW\\Computational Physics 1B - PH30056\\CW2\\Ising model - Copy (4)\\IMQ3.csv")
datax05 = pd.read_csv("C:\\Users\\gs010\\Desktop\\Bath\\Maths and Physics\\Year 3 CW\\Computational Physics 1B - PH30056\\CW2\\Ising model - Copy (5)\\IMQ3.csv")
datax06 = pd.read_csv("C:\\Users\\gs010\\Desktop\\Bath\\Maths and Physics\\Year 3 CW\\Computational Physics 1B - PH30056\\CW2\\Ising model - Copy (6)\\IMQ3.csv")

datax11 = pd.read_csv("C:\\Users\\gs010\\Desktop\\Bath\\Maths and Physics\\Year 3 CW\\Computational Physics 1B - PH30056\\CW2\\Ising model - Copy\\IMQ310.csv")
datax12 = pd.read_csv("C:\\Users\\gs010\\Desktop\\Bath\\Maths and Physics\\Year 3 CW\\Computational Physics 1B - PH30056\\CW2\\Ising model - Copy (2)\\IMQ310.csv")
datax13 = pd.read_csv("C:\\Users\\gs010\\Desktop\\Bath\\Maths and Physics\\Year 3 CW\\Computational Physics 1B - PH30056\\CW2\\Ising model - Copy (3)\\IMQ310.csv")
datax14 = pd.read_csv("C:\\Users\\gs010\\Desktop\\Bath\\Maths and Physics\\Year 3 CW\\Computational Physics 1B - PH30056\\CW2\\Ising model - Copy (4)\\IMQ310.csv")
datax15 = pd.read_csv("C:\\Users\\gs010\\Desktop\\Bath\\Maths and Physics\\Year 3 CW\\Computational Physics 1B - PH30056\\CW2\\Ising model - Copy (5)\\IMQ310.csv")
datax16 = pd.read_csv("C:\\Users\\gs010\\Desktop\\Bath\\Maths and Physics\\Year 3 CW\\Computational Physics 1B - PH30056\\CW2\\Ising model - Copy (6)\\IMQ310.csv")

datax21 = pd.read_csv("C:\\Users\\gs010\\Desktop\\Bath\\Maths and Physics\\Year 3 CW\\Computational Physics 1B - PH30056\\CW2\\Ising model - Copy\\IMQ350.csv")
datax22 = pd.read_csv("C:\\Users\\gs010\\Desktop\\Bath\\Maths and Physics\\Year 3 CW\\Computational Physics 1B - PH30056\\CW2\\Ising model - Copy (2)\\IMQ350.csv")
datax23 = pd.read_csv("C:\\Users\\gs010\\Desktop\\Bath\\Maths and Physics\\Year 3 CW\\Computational Physics 1B - PH30056\\CW2\\Ising model - Copy (3)\\IMQ350.csv")
datax24 = pd.read_csv("C:\\Users\\gs010\\Desktop\\Bath\\Maths and Physics\\Year 3 CW\\Computational Physics 1B - PH30056\\CW2\\Ising model - Copy (4)\\IMQ350.csv")
datax25 = pd.read_csv("C:\\Users\\gs010\\Desktop\\Bath\\Maths and Physics\\Year 3 CW\\Computational Physics 1B - PH30056\\CW2\\Ising model - Copy (5)\\IMQ350.csv")
datax26 = pd.read_csv("C:\\Users\\gs010\\Desktop\\Bath\\Maths and Physics\\Year 3 CW\\Computational Physics 1B - PH30056\\CW2\\Ising model - Copy (6)\\IMQ350.csv")

datax31 = pd.read_csv("C:\\Users\\gs010\\Desktop\\Bath\\Maths and Physics\\Year 3 CW\\Computational Physics 1B - PH30056\\CW2\\Ising model - Copy\\IMQ3100.csv")
datax32 = pd.read_csv("C:\\Users\\gs010\\Desktop\\Bath\\Maths and Physics\\Year 3 CW\\Computational Physics 1B - PH30056\\CW2\\Ising model - Copy (2)\\IMQ3100.csv")
datax33 = pd.read_csv("C:\\Users\\gs010\\Desktop\\Bath\\Maths and Physics\\Year 3 CW\\Computational Physics 1B - PH30056\\CW2\\Ising model - Copy (3)\\IMQ3100.csv")
datax34 = pd.read_csv("C:\\Users\\gs010\\Desktop\\Bath\\Maths and Physics\\Year 3 CW\\Computational Physics 1B - PH30056\\CW2\\Ising model - Copy (4)\\IMQ3100.csv")
datax35 = pd.read_csv("C:\\Users\\gs010\\Desktop\\Bath\\Maths and Physics\\Year 3 CW\\Computational Physics 1B - PH30056\\CW2\\Ising model - Copy (5)\\IMQ3100.csv")
datax36 = pd.read_csv("C:\\Users\\gs010\\Desktop\\Bath\\Maths and Physics\\Year 3 CW\\Computational Physics 1B - PH30056\\CW2\\Ising model - Copy (6)\\IMQ3100.csv")

datax1 = pd.concat([datax11, datax12, datax13, datax14, datax15, datax16])
datax1.to_csv('IMQ310.csv', index=False)

datax2 = pd.concat([datax21, datax22, datax23, datax24, datax25, datax26])
datax2.to_csv('IMQ350.csv', index=False)

datax3 = pd.concat([datax31, datax32, datax33, datax34, datax35, datax36])
datax3.to_csv('IMQ3100.csv', index=False)

# Concatenating the csvs from different seeds (4)

data1 = pd.read_csv("C:\\Users\\gs010\\Desktop\\Bath\\Maths and Physics\\Year 3 CW\\Computational Physics 1B - PH30056\\CW2\\Ising model - Copy\\IMQ4.csv")
data2 = pd.read_csv("C:\\Users\\gs010\\Desktop\\Bath\\Maths and Physics\\Year 3 CW\\Computational Physics 1B - PH30056\\CW2\\Ising model - Copy (2)\\IMQ4.csv")
data3 = pd.read_csv("C:\\Users\\gs010\\Desktop\\Bath\\Maths and Physics\\Year 3 CW\\Computational Physics 1B - PH30056\\CW2\\Ising model - Copy (3)\\IMQ4.csv")
data4 = pd.read_csv("C:\\Users\\gs010\\Desktop\\Bath\\Maths and Physics\\Year 3 CW\\Computational Physics 1B - PH30056\\CW2\\Ising model - Copy (4)\\IMQ4.csv")
data5 = pd.read_csv("C:\\Users\\gs010\\Desktop\\Bath\\Maths and Physics\\Year 3 CW\\Computational Physics 1B - PH30056\\CW2\\Ising model - Copy (5)\\IMQ4.csv")
data6 = pd.read_csv("C:\\Users\\gs010\\Desktop\\Bath\\Maths and Physics\\Year 3 CW\\Computational Physics 1B - PH30056\\CW2\\Ising model - Copy (6)\\IMQ4.csv")

data = pd.concat([data1, data2, data3, data4, data5, data6])

data.to_csv('IMQ4c.csv', index=False)

# Concatenating the csvs from different seeds (4Tc)

data1 = pd.read_csv("C:\\Users\\gs010\\Desktop\\Bath\\Maths and Physics\\Year 3 CW\\Computational Physics 1B - PH30056\\CW2\\Ising model - Copy\\IMQ4Tc.csv")
data2 = pd.read_csv("C:\\Users\\gs010\\Desktop\\Bath\\Maths and Physics\\Year 3 CW\\Computational Physics 1B - PH30056\\CW2\\Ising model - Copy (2)\\IMQ4Tc.csv")
data3 = pd.read_csv("C:\\Users\\gs010\\Desktop\\Bath\\Maths and Physics\\Year 3 CW\\Computational Physics 1B - PH30056\\CW2\\Ising model - Copy (3)\\IMQ4Tc.csv")
data4 = pd.read_csv("C:\\Users\\gs010\\Desktop\\Bath\\Maths and Physics\\Year 3 CW\\Computational Physics 1B - PH30056\\CW2\\Ising model - Copy (4)\\IMQ4Tc.csv")
data5 = pd.read_csv("C:\\Users\\gs010\\Desktop\\Bath\\Maths and Physics\\Year 3 CW\\Computational Physics 1B - PH30056\\CW2\\Ising model - Copy (5)\\IMQ4Tc.csv")
data6 = pd.read_csv("C:\\Users\\gs010\\Desktop\\Bath\\Maths and Physics\\Year 3 CW\\Computational Physics 1B - PH30056\\CW2\\Ising model - Copy (6)\\IMQ4Tc.csv")

data = pd.concat([data1, data2, data3, data4, data5, data6])

data.to_csv('IMQ4Tc.csv', index=False)

# Concatenating the csvs from different seeds (5)

data1 = pd.read_csv("C:\\Users\\gs010\\Desktop\\Bath\\Maths and Physics\\Year 3 CW\\Computational Physics 1B - PH30056\\CW2\\Ising model - Copy\\IMQ5.csv")
data2 = pd.read_csv("C:\\Users\\gs010\\Desktop\\Bath\\Maths and Physics\\Year 3 CW\\Computational Physics 1B - PH30056\\CW2\\Ising model - Copy (2)\\IMQ5.csv")
data3 = pd.read_csv("C:\\Users\\gs010\\Desktop\\Bath\\Maths and Physics\\Year 3 CW\\Computational Physics 1B - PH30056\\CW2\\Ising model - Copy (3)\\IMQ5.csv")
data4 = pd.read_csv("C:\\Users\\gs010\\Desktop\\Bath\\Maths and Physics\\Year 3 CW\\Computational Physics 1B - PH30056\\CW2\\Ising model - Copy (4)\\IMQ5.csv")
data5 = pd.read_csv("C:\\Users\\gs010\\Desktop\\Bath\\Maths and Physics\\Year 3 CW\\Computational Physics 1B - PH30056\\CW2\\Ising model - Copy (5)\\IMQ5.csv")
data6 = pd.read_csv("C:\\Users\\gs010\\Desktop\\Bath\\Maths and Physics\\Year 3 CW\\Computational Physics 1B - PH30056\\CW2\\Ising model - Copy (6)\\IMQ5.csv")

data = pd.concat([data1, data2, data3, data4, data5, data6])

data.to_csv('IMQ5c.csv', index=False)

# Concatenating the csvs from different seeds (7)

data1 = pd.read_csv("C:\\Users\\gs010\\Desktop\\Bath\\Maths and Physics\\Year 3 CW\\Computational Physics 1B - PH30056\\CW2\\Ising model - Copy\\IMQ7.csv")
data2 = pd.read_csv("C:\\Users\\gs010\\Desktop\\Bath\\Maths and Physics\\Year 3 CW\\Computational Physics 1B - PH30056\\CW2\\Ising model - Copy (2)\\IMQ7.csv")
data3 = pd.read_csv("C:\\Users\\gs010\\Desktop\\Bath\\Maths and Physics\\Year 3 CW\\Computational Physics 1B - PH30056\\CW2\\Ising model - Copy (3)\\IMQ7.csv")
data4 = pd.read_csv("C:\\Users\\gs010\\Desktop\\Bath\\Maths and Physics\\Year 3 CW\\Computational Physics 1B - PH30056\\CW2\\Ising model - Copy (4)\\IMQ7.csv")
data5 = pd.read_csv("C:\\Users\\gs010\\Desktop\\Bath\\Maths and Physics\\Year 3 CW\\Computational Physics 1B - PH30056\\CW2\\Ising model - Copy (5)\\IMQ7.csv")
data6 = pd.read_csv("C:\\Users\\gs010\\Desktop\\Bath\\Maths and Physics\\Year 3 CW\\Computational Physics 1B - PH30056\\CW2\\Ising model - Copy (6)\\IMQ7.csv")

data = pd.concat([data1, data2, data3, data4, data5, data6])

data.to_csv('IMQ7c.csv', index=False)

# Concatenating the csvs from different seeds (8)

data1 = pd.read_csv("C:\\Users\\gs010\\Desktop\\Bath\\Maths and Physics\\Year 3 CW\\Computational Physics 1B - PH30056\\CW2\\Ising model - Copy\\IMQ8.csv")
data2 = pd.read_csv("C:\\Users\\gs010\\Desktop\\Bath\\Maths and Physics\\Year 3 CW\\Computational Physics 1B - PH30056\\CW2\\Ising model - Copy (2)\\IMQ8.csv")
data3 = pd.read_csv("C:\\Users\\gs010\\Desktop\\Bath\\Maths and Physics\\Year 3 CW\\Computational Physics 1B - PH30056\\CW2\\Ising model - Copy (3)\\IMQ8.csv")
data4 = pd.read_csv("C:\\Users\\gs010\\Desktop\\Bath\\Maths and Physics\\Year 3 CW\\Computational Physics 1B - PH30056\\CW2\\Ising model - Copy (4)\\IMQ8.csv")
data5 = pd.read_csv("C:\\Users\\gs010\\Desktop\\Bath\\Maths and Physics\\Year 3 CW\\Computational Physics 1B - PH30056\\CW2\\Ising model - Copy (5)\\IMQ8.csv")
data6 = pd.read_csv("C:\\Users\\gs010\\Desktop\\Bath\\Maths and Physics\\Year 3 CW\\Computational Physics 1B - PH30056\\CW2\\Ising model - Copy (6)\\IMQ8.csv")

data = pd.concat([data1, data2, data3, data4, data5, data6])

data.to_csv('IMQ8c.csv', index=False)

# Concatenating the csvs from different seeds (8Tc)

data1 = pd.read_csv("C:\\Users\\gs010\\Desktop\\Bath\\Maths and Physics\\Year 3 CW\\Computational Physics 1B - PH30056\\CW2\\Ising model - Copy\\IMQ8Tc.csv")
data2 = pd.read_csv("C:\\Users\\gs010\\Desktop\\Bath\\Maths and Physics\\Year 3 CW\\Computational Physics 1B - PH30056\\CW2\\Ising model - Copy (2)\\IMQ8Tc.csv")
data3 = pd.read_csv("C:\\Users\\gs010\\Desktop\\Bath\\Maths and Physics\\Year 3 CW\\Computational Physics 1B - PH30056\\CW2\\Ising model - Copy (3)\\IMQ8Tc.csv")
data4 = pd.read_csv("C:\\Users\\gs010\\Desktop\\Bath\\Maths and Physics\\Year 3 CW\\Computational Physics 1B - PH30056\\CW2\\Ising model - Copy (4)\\IMQ8Tc.csv")
data5 = pd.read_csv("C:\\Users\\gs010\\Desktop\\Bath\\Maths and Physics\\Year 3 CW\\Computational Physics 1B - PH30056\\CW2\\Ising model - Copy (5)\\IMQ8Tc.csv")
data6 = pd.read_csv("C:\\Users\\gs010\\Desktop\\Bath\\Maths and Physics\\Year 3 CW\\Computational Physics 1B - PH30056\\CW2\\Ising model - Copy (6)\\IMQ8Tc.csv")

data = pd.concat([data1, data2, data3, data4, data5, data6])

data.to_csv('IMQ8Tc.csv', index=False)
