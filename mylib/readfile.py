import numpy as np

def read_sol(sol_name):
    
    file = open(sol_name, "r")
 
    # fisrt line
    line = file.readline().split()
    model = int(line[1])
    n = int(line[2])
    t = float(line[3])

    x=np.zeros(n)
    sol = np.zeros((model, n))

    for i, line in enumerate(file):
        lline = line.split()
        x[i] = float(lline[0])
        sol[0,i] = float(lline[1])
        sol[1,i] = float(lline[2])
        if model>2:
            sol[2,i] = float(lline[3])
        
    return model, t, x, sol
