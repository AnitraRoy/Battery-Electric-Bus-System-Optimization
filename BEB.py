#!/usr/bin/env python
# coding: utf-8

pip install gurobipy

from gurobipy import *
import gurobipy as gp
from itertools import product
import pandas as pd
import numpy as np

mod = Model(name="Deterministic Model")

f1 = pd.read_excel('or1.xlsx')
f2 = pd.read_excel('or2.xlsx')

D = f1['Depots'].unique()
D = np.delete(D, 1)

R = f1['Route no.'].unique()
c1 = f1['Origin terminal'].unique()
c2 = f1['Intermediate terminal'].unique()
c2 = np.delete(c2, -1)

T = np.unique(np.concatenate((c1, c2)))
T = np.delete(T, -6)

column_A = f1['Origin terminal']
column_B = f1['Intermediate terminal']
tr = np.column_stack((column_A, column_B))
TR = [tuple(row) for row in tr]

arr = f1.values

arr[9][2] = 'Badarpur Border'
arr[10][2] = 'Okhla Extn.'

RT = []
for i in T:
    x = []
    for j in arr:
        if j[2] == i or j[3] == i:
            x.append(j[1])
    RT.append(x)

for i in range(1, 6):
    arr[i][0] = 'Hari Nagar'
for i in range(7, 12):
    arr[i][0] = 'Sukhdev Vihar'
for i in range(13, 18):
    arr[i][0] = 'Subhash Place'

RD = []
for i in D:
    x = []
    for j in arr:
        if j[0] == i:
            x.append(j[1])
    RD.append(x)

dr = [i[0] for i in arr]
DR = range(len(dr))

BAT = [0, 1, 2]
SCS = [0, 1]
FCS = [0, 1, 2]
CAPb = [120, 250, 350]
CPsc = [20, 30]
CPfc = [150, 300, 450]
Cf1 = 2
Cf2 = 3
Cf3sc = [0.5, 0.7]
Cf4fc = [4.5, 8, 12]
Cb = [10, 20, 30]
Ce = 28
Eri = 2
ETECr = 1.5
SOCmax = 0.9
SOCmin = 0.2
NCTr = 5
ecdc = 0.05
CTc = 0

arr2 = f2.values

DHDr = f2['Deadhead(km)'].values
Hr = f2['Headway(min)'].values
Fr = f2['Trips per bus'].values
Lr = f2['Length(km)'].values
Orir = f1['Origin terminal'].values
Intr = f1['Intermediate terminal'].values
Dr = f1['Depots'].values
FSd = f2['Fleet size'].values

str = list(f2['Trips per bus'].values)

specific_pairs_e_arr = []
for i in range(0, len(R)):
    for j in range(1, str[i] + 1):
        specific_pairs_e_arr.append((R[i], Orir[i], j))
        specific_pairs_e_arr.append((R[i], Intr[i], j))
    specific_pairs_e_arr.append((R[i], Orir[i], 0))

specific_pairs_e_dep = []
for i in range(0, len(R)):
    for j in range(1, str[i] + 1):
        specific_pairs_e_dep.append((R[i], Orir[i], j))
        specific_pairs_e_dep.append((R[i], Intr[i], j))

specific_pairs_enr_t = []
for i in range(0, len(R)):
    specific_pairs_enr_t.append((R[i], Orir[i]))
    specific_pairs_enr_t.append((R[i], Intr[i]))

specific_pairs_enr_d = []
for i in range(0, len(R)):
    specific_pairs_enr_d.append((R[i], dr[i]))

x = mod.addVars(R, BAT, vtype=GRB.BINARY, name='x')
y = mod.addVars(D, SCS, vtype=GRB.BINARY, name='y')
z = mod.addVars(T, FCS, vtype=GRB.BINARY, name='z')

m = mod.addVars(D, SCS, vtype=GRB.INTEGER, lb=0, name='m')
n = mod.addVars(T, FCS, vtype=GRB.INTEGER, lb=0, name='n')

e_arr = mod.addVars(specific_pairs_e_arr, vtype=GRB.CONTINUOUS, lb=0, name='e_arr')
e_dep = mod.addVars(specific_pairs_e_dep, vtype=GRB.CONTINUOUS, lb=0, name='e_dep')
enr_t = mod.addVars(specific_pairs_enr_t, vtype=GRB.CONTINUOUS, lb=0, name='enr_t')
enr_d = mod.addVars(specific_pairs_enr_d, vtype=GRB.CONTINUOUS, lb=0, name='enr_d')

mod.addConstrs((gp.quicksum(x[r, b] for b in BAT) == 1 for r in R), name="constraint_2")

mod.addConstrs((e_arr[R[i], Orir[i], 0] == SOCmax * gp.quicksum(x[R[i], b] * CAPb[b] for b in BAT) - DHDr[i] * ETECr for i in range(len(R))), name="constraint_3")

mod.addConstrs((e_dep[R[i], Orir[i], 1] == e_arr[R[i], Orir[i], 0] for i in range(len(R))), name="constraint_4")

mod.addConstrs((e_arr[R[i], Intr[i], j] == e_dep[R[i], Orir[i], j] - (Eri * Lr[i]) / 2 for i in range(len(R)) for j in range(1, str[i] + 1)), name="constraint_5")

mod.addConstrs((e_dep[R[i], Intr[i], j] == e_arr[R[i], Intr[i], j] + enr_t[R[i], Intr[i]] for i in range(len(R)) for j in range(2, str[i] + 1)), name="constraint_6")

mod.addConstrs((e_arr[R[i], Orir[i], j] == e_dep[R[i], Intr[i], j] - (Eri * Lr[i]) / 2 for i in range(len(R)) for j in range(1, str[i] + 1)), name="constraint_7")

mod.addConstrs((e_dep[R[i], Orir[i], j] == e_arr[R[i], Orir[i], j - 1] + enr_t[R[i], Orir[i]] for i in range(len(R)) for j in range(2, str[i] + 1)), name="constraint_8")

mod.addConstrs((e_arr[R[i], Orir[i], str[i]] - DHDr[i] * ETECr + enr_d[R[i], dr[i]] == SOCmax * gp.quicksum(x[R[i], b] * CAPb[b] for b in BAT) for i in range(len(R))), name="constraint_9")

mod.addConstrs((SOCmin * gp.quicksum(x[i[0], b] * CAPb[b] for b in BAT) <= e_arr[i] for i in specific_pairs_e_arr), name="constraint_10")

mod.addConstrs((SOCmax * gp.quicksum(x[i[0], b] * CAPb[b] for b in BAT) >= e_arr[i] for i in specific_pairs_e_arr), name="constraint_11")

mod.addConstrs((SOCmin * gp.quicksum(x[i[0], b] * CAPb[b] for b in BAT) <= e_dep[i] for i in specific_pairs_e_dep), name="constraint_12")

mod.addConstrs((SOCmax * gp.quicksum(x[i[0], b] * CAPb[b] for b in BAT) >= e_dep[i] for i in specific_pairs_e_dep), name="constraint_13")

nR = {R[i]: f2[f2['Route no.'] == R[i]]['Fleet size'].values[0] for i in range(len(R))}

mod.addConstrs((gp.quicksum(nR[j] * (enr_d[j, D[i]] / NCTr) for j in RD[i]) <= gp.quicksum(CPsc[k] * m[D[i], k] for k in SCS) for i in range(len(D))), name="constraint_14")

mod.addConstrs((m[D[i], sc] <= 100 * y[D[i], sc] for i in range(len(D)) for sc in SCS), name="constraint_15")

mod.addConstrs((enr_t[r, T[t]] <= 100000 * gp.quicksum(z[T[t], fc] for fc in FCS) for t in range(len(T)) for r in RT[t]), name="constraint_16")

mod.addConstrs((n[T[t], fc] <= len(RT[t]) * z[T[t], fc] for fc in FCS for t in range(len(T))), name="constraint_17")

mod.addConstrs((gp.quicksum(z[t, fc] for fc in FCS) <= 1 for t in T), name="constraint_18")

xyz = {R[i]: f2[f2['Route no.'] == R[i]]['Trips per bus'].values[0] for i in range(len(R))}
no_of_buses = {R[i]: f2[f2['Route no.'] == R[i]]['Fleet size'].values[0] for i in range(len(R))}

mod.setObjective(
    gp.quicksum(Cb[b] * x[r, b] for b in BAT for r in R) +
    Cf1 * gp.quicksum(y[d, sc] for sc in SCS for d in D) +
    Cf2 * gp.quicksum(z[t, fc] for fc in FCS for t in T) +
    gp.quicksum(Cf3sc[sc] * m[d, sc] for sc in SCS for d in D) +
    gp.quicksum(Cf4fc[fc] * n[t, fc] for fc in FCS for t in T) +
    Ce * 365 * gp.quicksum(enr_d[r, D[d]] * no_of_buses[r] for d in range(len(D)) for r in RD[d]) +
    Ce * 365 * gp.quicksum(enr_t[r, T[t]] * xyz[r] * no_of_buses[r] for t in range(len(T)) for r in RT[t]),
    GRB.MINIMIZE
)

mod.optimize()

if mod.status == gp.GRB.Status.OPTIMAL:
    objective_value = mod.objVal
    print("Optimal objective value: $", objective_value)

    for v in mod.getVars():
        print(f"{v.varName}: {v.x}")

    print(f"Computational time: {mod.Runtime} seconds")
else:
    print("No optimal solution found.")
