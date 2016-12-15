from __future__ import print_function
import re
import sys
from prettytable import PrettyTable
import os

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


t = PrettyTable(['Q#', 'nMatch','nMisMatch' ,'nIns' ,'nDel' ,'Sim#' ,'RawScore'])
for i in range(0,len(query_name_list)):
    t.add_row([query_name_list[i],nMatch_list[i],nMisMatch_list[i],nIns_list[i],nDel_list[i],sim_list[i],raw_score_list[i]])
sys.stdout = open(os.path.dirname(sys.argv[1])+"/scores.txt",'w')

print("Max_nMatch :", max_nMatch,"   Max_nMisMatch:", max_nMisMatch,"  Max_nIns:", max_nIns,"  Max_nDel:", max_nDel,"  Max_Sim:", max_sim,"  Max_RawScore:", max_raw_score)
print("Avg_nMatch :",round(total_nMatch/counter_nMatch,3),"  Avg_nMisMatch :",round(total_nMisMatch/counter_nMisMatch,3),"  Avg_nIns :",round(total_nIns/counter_nIns,3),"  Max_nDel :",round(total_nDel/counter_nDel,3),"  Avg_Sim :",round(total_sim/counter_sim,3),"  Avg_RawScore :",round(total_raw_score/counter_raw_score,3))
print(t)
