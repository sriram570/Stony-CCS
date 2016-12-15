import re
import sys
from prettytable import PrettyTable
import matplotlib.pyplot as plt


max_nMisMatch =  sys.maxsize
max_nIns = sys.maxsize
max_nDel = sys.maxsize
max_raw_score = sys.maxsize
max_nMatch = -1
max_sim = -1
counter_nMatch = 0
counter_nMisMatch = 0
counter_nIns = 0
counter_nDel = 0
counter_sim = 0
counter_raw_score = 0
total_nMatch  = 0
total_nMisMatch = 0
total_nIns = 0
total_nDel = 0
total_sim = 0
total_raw_score = 0
nMatch_list = []
sim_list = []
nDel_list = []
nIns_list = []
raw_score_list = []
nMisMatch_list = []
query_name_list = []
with open(sys.argv[1]) as f:
    for line in f:
        a = line
        a = re.sub(r'\s+', '', a)
        a = re.split(':', a)
        if a[0] == 'nMatch':
            counter_nMatch += 1
            max_nMatch = max(max_nMatch,int(a[1]))
            total_nMatch += int(a[1])
            nMatch_list.append(a[1])
        if a[0] == 'nMisMatch':
            counter_nMisMatch += 1
            max_nMisMatch = min(max_nMisMatch,int(a[1]))
            total_nMisMatch += int(a[1])
            nMisMatch_list.append(a[1])
        if a[0] == 'nIns':
            counter_nIns += 1
            max_nIns = min(max_nIns,int(a[1]))
            total_nIns += int(a[1])
            nIns_list.append(a[1])
        if a[0] == 'nDel':
            counter_nDel += 1
            max_nDel = min(max_nDel,int(a[1]))
            total_nDel += int(a[1])
            nDel_list.append(a[1])
        if a[0] == '%sim':
            counter_sim += 1
            max_sim = max(max_sim,float(a[1]))
            total_sim += float(a[1])
            sim_list.append(a[1])
        if a[0] == 'Rawscore':
            counter_raw_score += 1
            max_raw_score = min(max_raw_score,int(a[1]))
            total_raw_score += int(a[1])
            raw_score_list.append(a[1])
        if a[0] == 'Query':
            query_name_list.append(a[1])

e = []
x = []
y = []
z = []
for i in range(0,len(nMisMatch_list)):
    if int(nMisMatch_list[i])!=0 :
        x.append(round(int(nMatch_list[i])/int(nMisMatch_list[i]),3))
        y.append(round(float(sim_list[i]),3))
        z.append(int(nMatch_list[i]))
        t = (round(int(nMatch_list[i])/int(nMisMatch_list[i]),3),sim_list[i])
        e.append(t)
#print(e)
#print(x)
#print(y)


t = sorted(zip(x,y))
#print(t)
x,y = zip(*t)
#ylim=(0.0, 1.01)
xlim = (0,20)
plt.figure()
plt.title("Match/Mismatch ratio vs SimScore")
plt.xlabel("Match/Mismatch ratio")
plt.ylabel("Sim Score")
plt.grid()
# if ylim is not None:
#     plt.ylim(*ylim)
if xlim is not None:
    plt.xlim(*xlim)
plt.plot(x, y, 'o-', color="r",
         label="Test config")
plt.legend(loc="best")
plt.show()

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


