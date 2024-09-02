"""
Generating random trees for the paper "A Vector Representation for Phylogenetic Trees"
"""

import sys
import random
from datetime import datetime
random.seed(datetime.now().timestamp())

N = int(sys.argv[1]) # Number of taxa
M = int(sys.argv[1]) # Number of trees
OUT_FILE = sys.argv[3] # Output file

def create_tree():
    BRANCHES = {
        i: None
        for i in range(1,2*N)
    };

    LAST_BRANCH = 1;
    BRANCHES[1] = [0,1]

    for i in range(2,N+1):
        j = random.randint(1,2*(i-1)-1)
        X = BRANCHES[j][1]
        BRANCHES[j][1] = N-1+i
        LAST_BRANCH += 1
        BRANCHES[LAST_BRANCH] = [N-1+i, i]
        LAST_BRANCH += 1
        BRANCHES[LAST_BRANCH] = [N-1+i, X]

        CHILDREN = {
            f"t{i}": []
            for i in range(1,2*N)
        };
    
    for i,b in BRANCHES.items():
        BRANCHES[i] = [f"t{b[0]}",f"t{b[1]}",random.random()]
        if b[0] > 0:
            CHILDREN[BRANCHES[i][0]].append([BRANCHES[i][1], BRANCHES[i][2]])
            
    return BRANCHES,CHILDREN
                
NEWICK_STR = ""
def write_Newick(BRANCHES, CHILDREN):
    """
    Assumption: binary trees
    """
    global NEWICK_STR
    NEWICK_STR = ""
    def write_Newick_rec(X):
        global NEWICK_STR
        if CHILDREN[X] == []:
            NEWICK_STR += X
        else:
            NEWICK_STR += "("
            Y1 = CHILDREN[X][0]
            write_Newick_rec(Y1[0])        
            NEWICK_STR += f":{str(Y1[1])}"
            NEWICK_STR += ","
            Y2 = CHILDREN[X][1]
            write_Newick_rec(Y2[0])        
            NEWICK_STR += f":{str(Y2[1])}"
            NEWICK_STR += ")"
    write_Newick_rec(BRANCHES[1][1])
    return f"{NEWICK_STR};"   


with open(OUT_FILE,"w") as out_file:
    for i in range(M):
        BRANCHES,CHILDREN = create_tree()
        out_file.write(write_Newick(BRANCHES, CHILDREN))
        out_file.write("\n")




    
