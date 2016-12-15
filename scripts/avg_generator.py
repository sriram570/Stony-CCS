from __future__ import print_function

import re
import sys

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

print("Max_nMatch :", max_nMatch,"   Max_nMisMatch:", max_nMisMatch,"  Max_nIns:", max_nIns,"  Max_nDel:", max_nDel,"  Max_Sim:", max_sim,"  Max_RawScore:", max_raw_score)

print("Avg_nMatch :",total_nMatch/counter_nMatch,"Avg_nMisMatch :",total_nMisMatch/counter_nMisMatch,"Avg_nIns :",total_nIns/counter_nIns,"Max_nDel :",total_nDel/counter_nDel,"Avg_Sim :",total_sim/counter_sim,"Avg_RawScore :",total_raw_score/counter_raw_score)

print("Q#\t\t\t\t\t\t\t nMatch\t nMisMatch\t nIns\t nDel\t Sim%\t RawScore")
for i in range(0,len(query_name_list)):
    print(query_name_list[i],'\t',nMatch_list[i],'\t',nMisMatch_list[i],'\t\t',nIns_list[i],'\t',nDel_list[i],'\t',sim_list[i],'\t',raw_score_list[i])
