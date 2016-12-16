import re
import sys
import matplotlib.pyplot as plt


nMatch_list = []
sim_list = []
nMisMatch_list = []
nMatch_list1 = []
sim_list1 = []
nMisMatch_list1 = []
nMatch_list2 = []
sim_list2 = []
nMisMatch_list2 = []


with open(sys.argv[1]) as f:
    for line in f:
        a = line
        a = re.sub(r'\s+', '', a)
        a = re.split(':', a)
        if a[0] == 'nMatch':
            nMatch_list.append(a[1])
        if a[0] == 'nIns':
            nMisMatch_list.append(a[1])
        if a[0] == '%sim':
            sim_list.append(a[1])


with open(sys.argv[2]) as f:
    for line in f:
        a = line
        a = re.sub(r'\s+', '', a)
        a = re.split(':', a)
        if a[0] == 'nMatch':
            nMatch_list1.append(a[1])
        if a[0] == 'nIns':
            nMisMatch_list1.append(a[1])
        if a[0] == '%sim':
            sim_list1.append(a[1])

with open(sys.argv[3]) as f:
    for line in f:
        a = line
        a = re.sub(r'\s+', '', a)
        a = re.split(':', a)
        if a[0] == 'nMatch':
            nMatch_list2.append(a[1])
        if a[0] == 'nIns':
            nMisMatch_list2.append(a[1])
        if a[0] == '%sim':
            sim_list2.append(a[1])

e = []
x = []
y = []
z = []
x1 = []
y1 = []
z1 = []
x2 = []
y2 = []
z2 = []

for i in range(0,len(nMisMatch_list)):
    if int(nMisMatch_list[i])!=0 :
        x.append(round(int(nMatch_list[i])/int(nMisMatch_list[i]),3))
        y.append(round(float(sim_list[i]),3))
        z.append(int(nMatch_list[i]))

for i in range(0, len(nMisMatch_list1)):
    if int(nMisMatch_list1[i]) != 0:
        x1.append(round(int(nMatch_list1[i]) / int(nMisMatch_list1[i]), 3))
        y1.append(round(float(sim_list1[i]), 3))
        z1.append(int(nMatch_list1[i]))

for i in range(0, len(nMisMatch_list2)):
    if int(nMisMatch_list2[i]) != 0:
        x2.append(round(int(nMatch_list2[i]) / int(nMisMatch_list2[i]), 3))
        y2.append(round(float(sim_list2[i]), 3))
        z2.append(int(nMatch_list2[i]))


t = sorted(zip(x,y))
t1 = sorted(zip(x1,y1))
t2 = sorted(zip(x2,y2))

x,y = zip(*t)
x1,y1 = zip(*t1)
x2,y2 = zip(*t2)
#ylim=(0.0, 1.01)
xlim = (0,20)
plt.figure()
plt.title("Comaprison of Baseline vs Star_Forward_Reverse ordering algo with Edge_Weight scoring fn")
plt.xlabel("Match/Insertion ratio")
plt.ylabel("Sim Score")
plt.grid()
# if ylim is not None:
#     plt.ylim(*ylim)
if xlim is not None:
    plt.xlim(*xlim)
plt.plot(x, y, '-', color="b",
         label="pbccs (Baseline)")
plt.plot(x1, y1, '-', color="r",
         label="MaxScore traversal algorithm")
plt.plot(x2, y2, '-', color="g",
         label="Opt_Random traversal algorithm")
plt.legend(loc="best")
plt.show()

# x,y = zip(*t)
# x1,y1 = zip(*t1)
# x2,y2 = zip(*t2)
# #ylim=(0.0, 1.01)
# xlim = (0,20)
# plt.figure()
# plt.title("Comaprison of Baseline vs Star_Forward_Reverse ordering algo with Edge_Weight scoring fn")
# plt.xlabel("Match/MisMatch ratio")
# plt.ylabel("Sim Score")
# plt.grid()
# # if ylim is not None:
# #     plt.ylim(*ylim)
# if xlim is not None:
#     plt.xlim(*xlim)
# plt.plot(x, y, '-', color="b",
#          label="pbccs (Baseline)")
# plt.plot(x1, y1, '-', color="r",
#          label="MaxScore traversal algorithm")
# plt.plot(x2, y2, '-', color="g",
#          label="Opt_Random traversal algorithm")
# plt.legend(loc="best")
# plt.show()

# t1 = sorted(zip(z,y))
# print(t)
# z,y = zip(*t1)
# #ylim=(0.0, 1.01)
# xlim = (0,1000)
# plt.figure()
# plt.title("Match  vs SimScore")
# plt.xlabel("Match length")
# plt.ylabel("Sim Score")
# plt.grid()
# # if ylim is not None:
# #     plt.ylim(*ylim)
# if xlim is not None:
#     plt.xlim(*xlim)
# plt.plot(z, y, 'o-', color="r",
#          label="Test config2")
# plt.legend(loc="best")
# plt.show()


