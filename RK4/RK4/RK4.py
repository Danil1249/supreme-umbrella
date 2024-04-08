from turtle import title
from matplotlib import pyplot as plt
import math
import numpy as np

def RK4(z, q, vv, mm, nn, hh):
  t = [0] * (z + 1)
  V = [0] * (z + 1)
  m = [0] * (z + 1)
  n = [0] * (z + 1)
  h = [0] * (z + 1)
  
  
  Box = [[], [], [], [], []]
  t[0] = 0
  Box[0].append(t[0])
  V[0] = vv
  Box[1].append(V[0])
  m[0] = mm
  Box[2].append(m[0])
  n[0] = nn
  Box[3].append(n[0])
  h[0] = hh
  Box[4].append(h[0])
  for i in range(1, z + 1):
      
      mk1 = q * func_m(m[i-1], V[i-1])
      nk1 = q * func_n(n[i-1], V[i-1])
      hk1 = q * func_h(h[i-1], V[i-1])
      Vk1 = q * func_V(m[i-1], n[i-1], h[i-1], V[i-1])
      
      
      mk2 = q * func_m(m[i-1] + mk1/2, V[i-1])
      nk2 = q * func_n(n[i-1] + nk1/2, V[i-1])
      hk2 = q * func_h(h[i-1] + hk1/2, V[i-1])
      Vk2 = q * func_V(m[i-1], n[i-1], h[i-1], V[i-1] + Vk1/2)
      
      
      mk3 = q * func_m(m[i-1] + mk2/2, V[i-1])
      nk3 = q * func_n(n[i-1] + nk2/2, V[i-1])
      hk3 = q * func_h(h[i-1] + hk2/2, V[i-1])
      Vk3 = q * func_V(m[i-1], n[i-1], h[i-1], V[i-1] + Vk2/2)
      
      
      mk4 = q * func_m(m[i-1] + mk3, V[i-1])
      nk4 = q * func_n(n[i-1] + nk3, V[i-1])
      hk4 = q * func_h(h[i-1] + hk3, V[i-1])
      Vk4 = q * func_V(m[i-1], n[i-1], h[i-1], V[i-1] + Vk3)
        
      t[i] = t[i-1] + q
      Box[0].append(t[i])
      
      m[i] = m[i-1] + (mk1 + 2*mk2 + 2*mk3 + mk4) / 6
      n[i] = n[i-1] + (nk1 + 2*nk2 + 2*nk3 + nk4) / 6
      h[i] = h[i-1] + (hk1 + 2*hk2 + 2*hk3 + hk4) / 6
      V[i] = V[i-1] + (Vk1 + 2*Vk2 + 2*Vk3 + Vk4) / 6

      Box[1].append(V[i])
      Box[2].append(m[i])
      Box[3].append(n[i])
      Box[4].append(h[i])

  return Box
    
def func_m(m, V):
  return 1000 * (alpha_m(V) * (1 - m) - beta_m(V) * m)

def func_n(n, V):
  return 1000 * (alpha_n(V) * (1 - n) - beta_n(V) * n)

def func_h(h, V):
  return 1000 * (alpha_h(V) * (1 - h) - beta_h(V) * h)
   
def func_V(m, n, h, V):
  return 1000 * ((gNa * m**3 * h * (ENa - V) + gK * n * (EK - V) + gLeak * (ELeak - V) + iapp)) / C

def alpha_m(V):
  return 0.182 * (V + 35) / (1 - np.exp(-1*(V + 35) / 9))

def beta_m(V):
  return -0.124 * (V + 35) / (1 - np.exp((V + 35)/9))

def alpha_n(V):
  return 0.02 * (V - 25) / (1 - np.exp(-1 * (V - 25)/9))

def beta_n(V):
  return -0.002 * (V - 25) / (1 - np.exp((V - 25)/9))

def alpha_h(V):
  return 0.25 * np.exp(-1 * (V + 90)/12)

def beta_h(V):
  return 0.25 * np.exp((V + 62) / 6) / np.exp((V + 90)/12)

def Maximum():
  Cage = [[], []]
  for i in range(0, len(Box[1]) - 1):
    if (Box[1][i] > Box[1][i + 1]) and (Box[1][i] > Box[1][i - 1]):
      Cage[0].append(Box[0][i])
      Cage[1].append(Box[1][i])
  return Cage

def Period():
  Container = [[], []]
  for i in range(0, len(Cage[0]) - 1):
    Container[0].append(Cage[0][i])
    Container[1].append(Cage[0][i + 1] - Cage[0][i])
  return Container

def Frequency():
  Container = [[], []]
  for i in range(0, len(Cage[0]) - 1):
    Container[0].append(Cage[0][i])
    Container[1].append(1 / (Cage[0][i + 1] - Cage[0][i]))
  return Container

def Middle_Frequency():
  sum = 0
  for i in Container[1]:
    sum += i
  if len(Container[1]) != 0:
    ssum = sum / len(Container[1])
    return 1 / ssum
  else:
    return 0.0

C = 1
gNa = 40
gK = 35
gLeak = 0.3
ENa = 55
EK = -77
ELeak = -65
iapp = 1.0361


z = 1000000
q = 0.00001

file = open('initial_conditions.txt', 'r')
data = []

while True:
  line = file.readline()  
  if not line:
    break
  data.append(float(line))

file.close()

vv = data[0]
mm = data[1]
nn = data[2]
hh = data[3]

# # data1 = []
# # file6 = open('Iapp.txt', 'r')

# # while True:
# #   line = file6.readline()  
# #   d = ''
# #   if not line:
# #     break
# #   for i in line:
# #     if i == ' ' or i == '\n' or i == '\r':
# #       data1.append(float(d))
# #       d = ''
# #     else:
# #       d += i
      
# # data2 = [[], []]

# # for i in range(len(data1)):
# #   if i % 2 == 0:
# #     data2[0].append(data1[i])
# #   if i % 2 != 0:
# #     data2[1].append(data1[i])
  
# # file6.close()

# # data1 = []
# # file6 = open('Iapp2.txt', 'r')

# # while True:
# #   line = file6.readline()  
# #   d = ''
# #   if not line:
# #     break
# #   for i in line:
# #     if i == ' ' or i == '\n' or i == '\r':
# #       data1.append(float(d))
# #       d = ''
# #     else:
# #       d += i
      
# # data3 = [[], []]

# # for i in range(len(data1)):
# #   if i % 2 == 0:
# #     data3[0].append(data1[i])
# #   if i % 2 != 0:
# #     data3[1].append(data1[i])
  
# # file6.close()

Box = RK4(z, q, vv, mm, nn, hh)
# Box[0] = Box[0][int(0.25*len(Box[0])):len(Box[0])]
# Box[1] = Box[1][int(0.25*len(Box[1])):len(Box[1])]
# Box[2] = Box[2][int(0.25*len(Box[2])):len(Box[2])]
# Box[3] = Box[3][int(0.25*len(Box[3])):len(Box[3])]
# Box[4] = Box[4][int(0.25*len(Box[4])):len(Box[4])]
Cage = Maximum()
Container = Period()
Storage = Frequency()
middle = Middle_Frequency()



# # file = open('V(t).txt', 'w')
# # file1 = open('Maximum.txt', 'w')
# # file2 = open('Period.txt', 'w')
# # file3 = open('Frequency.txt', 'w')
file4 = open('initial_conditions.txt', 'w')
file5 = open('Iapp.txt', 'a')

# # for i in range(len(Box[1])):
# #   if i < len(Box[1]) - 1:
# #     g = str(Box[0][i])
# #     p = str(Box[1][i])
# #     file.write(g + ' ' + p + '\n')
# #   else:
# #     g = str(Box[0][i])
# #     p = str(Box[1][i])
# #     file.write(g + ' ' + p)
    
# # for i in range(len(Cage[1])):
# #   if i < len(Cage[1]) - 1:
# #     g = str(Cage[0][i])
# #     p = str(Cage[1][i])
# #     file1.write(g + ' ' + p + '\n')
# #   else:
# #     g = str(Cage[0][i])
# #     p = str(Cage[1][i])
# #     file1.write(g + ' ' + p)
    
# # for i in range(len(Container[1])):
# #   if i < len(Container[1]) - 1:
# #     g = str(Container[0][i])
# #     p = str(Container[1][i])
# #     file2.write(g + ' ' + p + '\n')
# #   else:
# #     g = str(Container[0][i])
# #     p = str(Container[1][i])
# #     file2.write(g + ' ' + p)
    
# # for i in range(len(Storage[1])):
# #   if i < len(Storage[1]) - 1:
# #     g = str(Storage[0][i])
# #     p = str(Storage[1][i])
# #     file3.write(g + ' ' + p + '\n')
# #   else:
# #     g = str(Storage[0][i])
# #     p = str(Storage[1][i])
# #     file3.write(g + ' ' + p)

p = str(Box[1][-1])
j = str(Box[2][-1])
l = str(Box[3][-1])
w = str(Box[4][-1])
file4.write(p + '\n')
file4.write(j + '\n')
file4.write(l + '\n')
file4.write(w)

file5.write(str(iapp) + ' ' + str(middle) + '\n')

# plt.style.use('ggplot')
# plt.figure(figsize=(9,6)) 
# plt.subplots_adjust(hspace=0.4)

# plt.subplot(3, 1, 1, title="Решение")
# # plt.xlabel("t")
# # plt.ylabel("V(t)")
# plt.plot(Box[0], Box[1], "red")

# # plt.savefig('Graphic.png')


# # ax = plt.axes()
# # plt.xlabel("t")
# # plt.ylabel("V(t)")
# # plt.plot(Box[0], Box[1], "green")
# # for i in range(len(Cage[0])):
# #   ax.plot(Cage[0][i], Cage[1][i], '-ro', label='marker only')
# # plt.title("Точки максимума функции V(t)")
# # plt.savefig('Maximum.png')

# plt.subplot(3, 1, 2, title="Периоды колебаний", ylim=(0, 0.1))

# # plt.xlabel("t")
# # plt.ylabel("T(t)")
# plt.plot(Container[0], Container[1], "blue")

# # plt.savefig('Period.png')

# plt.subplot(3, 1, 3, title="Частоты колебаний", ylim=(0, 100))
# # plt.xlabel("t")
# # plt.ylabel("Frequency(t)")
# plt.plot(Storage[0], Storage[1], "blue")
# plt.savefig('3Graphics21.png')
# plt.show()

plt.plot(Box[4], Box[1])
plt.title("Iapp = 1.0369")
plt.xlabel('h')
plt.ylabel('V, mV')
plt.show()

# # plt.ylabel("Частоты")
# # plt.xlabel("Iapp")
# # plt.xticks(np.arange(0.0, 2.1, 0.2))
# # plt.plot(data2[0], data2[1], "red")
# # plt.plot(data3[0], data3[1], "blue")
# # plt.savefig('Iapp.png')
# # plt.show()


# file.close()
# file1.close()
# file2.close()
# file3.close()
# # plt.close()
file4.close()
file5.close()

