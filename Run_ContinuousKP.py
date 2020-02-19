"""
author: Margarida Carvalho
License: MIT
PYTHON 2
"""


# Limitation in the number of cores
from gurobipy import *
setParam("Threads", 1)
setParam("OutputFlag", 0) # gurobi shows no output

# do not forget to order items
def PolynomialTimeSolver(P, F, L, c_F, c_L, n):
    best = sum(P)
    k = -1
    for c in range(1,n+1):
        # Solve CBK_c
        OPT_c, feasible, x_opt = Solve_CBK(P, F, L, c_F, c_L,n, c)
        if feasible:
            val = OPT_c
            if val < best:
                best = val
                x_best = x_opt.copy()
                k = c
    # no item is critical
    OPT_c, feasible, x_opt = Solve_CBK(P, F, L, c_F, c_L,n,n+1)
    if feasible:
        val = OPT_c
        if val < best:
            best = val
            x_best = x_opt.copy()
    return best, x_best

def Solve_CBK(P, F, L, c_F, c_L,n, c):
    CBK = Model("item c is critical")
    if c<=n:
        x = {i: CBK.addVar(vtype="C",ub=1,lb=0,obj=-P[i-1]+((P[c-1]*1.)/F[c-1])*F[i-1],name="x"+str(i)) for i in range(1,c)}
    else:
        x = {i: CBK.addVar(vtype="C",ub=1,lb=0,obj=-P[i-1],name="x"+str(i)) for i in range(1,c)}
    CBK.update()
    CBK.addConstr(quicksum(L[i-1]*x[i] for i in range(1,c)) <= c_L)
    CBK.update()
    CBK.addConstr(quicksum(F[i-1]*(1-x[i]) for i in range(1,c)) <= c_F)
    CBK.update()
    if c<=n:
        CBK.addConstr(F[c-1]+quicksum(F[i-1]*(1-x[i]) for i in range(1,c)) >= c_F)
        CBK.update()
    CBK.ModelSense = 1 # to minimize
    CBK.optimize()
    if CBK.status == 3: #infeasible
        return 0, False, []
    else:
        if c<=n:
            return CBK.ObjVal+((P[c-1]*1.)/F[c-1])*c_F+sum(P[i-1] - ((P[c-1]*1.)/F[c-1])*F[i-1] for i in range(1,c)), True, {i: (x[i].x if i < c else 0) for i in range(1,n+1)}
        else:
            return CBK.ObjVal+sum(P[i-1] for i in range(1,c)), True, {i: (x[i].x if i <c else 0) for i in range(1,n+1)}

def Order(P,F,L,n):
    ratio = [((P[i]*1.)/F[i],i) for i in range(n)]
    ratio.sort()
    ratio.reverse()
    ratio = [a[1] for a in ratio]
    return [P[i] for i in ratio], [F[i] for i in ratio], [L[i] for i in ratio]

#################################

if __name__ == "__main__":
    ################## Run a single instance #####################################
    n = 35
    P = [19, 87, 97, 22, 47, 22, 30, 5, 32, 54, 75, 70, 7, 76, 31, 61, 98, 27, 39, 91, 50, 10, 38, 25, 43, 37, 77, 18, 99, 9, 19, 97, 75, 26, 41]
    F = [1, 96, 67, 90, 13, 74, 22, 86, 23, 63, 89, 25, 100, 76, 23, 56, 5, 47, 45, 39, 47, 11, 100, 27, 60, 76, 75, 37, 1, 81, 51, 30, 40, 48, 68]
    L = [14, 85, 77, 26, 50, 45, 66, 79, 10, 3, 84, 44, 77, 1, 45, 73, 23, 95, 91, 4, 3, 55, 94, 39, 22, 43, 3, 23, 44, 50, 24, 24, 22, 46, 29]
    c_L = 152
    c_F = 162
    OPT_continuous, xOpt_continuous = PolynomialTimeSolver(P, F, L, c_F, c_L, n)
