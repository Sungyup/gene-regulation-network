'''
Created on 2015. 7. 14.

@author: Sungyup
'''

import networkx as nx
import matplotlib.pyplot as plt

def inferRegulation(g, gene1, gene2, f):

    paths = nx.all_simple_paths(g, gene1, gene2, 5)
    
    if nx.has_path(g, gene1, gene2) == False:
        return 
    if gene2 in g[gene1].keys():
        weight = 0
        for path in paths:
            for i in range(len(path) - 1):
                weight += g[path[i]][path[i + 1]]['weight']
                if (g[gene1][gene2]['weight'] != weight):
                    print ("cannot infer")
                    # gc.collect()
                    return
                weight = 0
        if g[gene1][gene2]['weight'] % 2 == 1:
            f.write("inferred down regulation : " + gene1 + ", " + gene2 + "\n")
            print("inferred down regulation : " + gene1 + ", " + gene2 + "\n")
        else:
            f.write("inferred up regulation : " + gene1 + ", " + gene2 + "\n")
            print("inferred up_regulation : " + gene1 + ", " + gene2 + "\n")
                             
    else:
        weight_old = 0
        weight_new = 0
        first = True
       
        for path in paths:    
            # print (path)
            for i in range(len(path) - 1):
                weight_new += g[path[i]][path[i + 1]]['weight']
                
            if first == True:
                weight_old = weight_new
                first = False
                continue
            
            if weight_new % 2 != weight_old % 2:
                print ("cannot infer")
                return
            else:
                weight_new = 0
        
        # print (weight_old)
        
        if weight_old % 2 == 1:
            print("inferred down regulation : " + gene1 + ", " + gene2 + "\n")
            f.write("inferred down regulation : " + gene1 + ", " + gene2 + "\n")
            # gc.collect()
            return 
        else:
            print("inferred up_regulation : " + gene1 + ", " + gene2 + "\n")
            f.write("inferred up_regulation : " + gene1 + ", " + gene2 + "\n")
            # gc.collect()
            return 

def inferRegulation_all(g):
    nodes = nx.nodes(g)
    
    f = open("inference result.txt", 'w')
    
    for gene1 in nodes:
        #gc.collect()
        for gene2 in nodes:
            if gene2 != gene1:
                inferRegulation(g, gene1, gene2, f)
def drawGraph(g):
    nx.draw(g)
    plt.savefig("gene_regulation.png") # save as png        

def main():
    f = open("inference_rules_canonical_reg.pl", 'r')
    
    lines = f.readlines()
    
    g = nx.DiGraph()
    
    for line in lines:
        index = line.find("(")
        gene = line[index + 1:]
        gene = gene.split(",")
        
        gene[1] = gene[1][1:gene[1].find(")") - 1]
    
        g.add_edge(gene[0], gene[1])    
        
        if line.find("positive_regulation") != -1:
            g[gene[0]][gene[1]]['weight'] = 1
        else:
            g[gene[0]][gene[1]]['weight'] = 0
    

    #inferRegulation_all(g)

if __name__ == "__main__":
    main()
