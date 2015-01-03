# -*- coding: utf-8 -*-
"""
Created on Thu Dec 18 10:50:22 2014

@author: rkelley
"""

from gurobipy import *
import numpy as np

time_cutoff = 10

def mycallback(model, where):
    if where == GRB.Callback.MIP:
        time = model.cbGet(GRB.Callback.RUNTIME)
        best = model.cbGet(GRB.Callback.MIP_OBJBST)
        if time > time_cutoff and best < GRB.INFINITY:
            model.terminate()

def build_model(c,Q,A,rhs,sense,lb,ub,vtype):

    # Maximizes:
    #
    #  f(x) = x.T Q x + c.T x
    #
    #  subject to: 
    #       Ax ? rhs, where ? can be <=, =, >=, <, > for each row
    #
    
    # c, numpy 1D matrix for the linear term

    # Q, numpy 2D matrix for the quadratic term
    # if Q is None, then it will build a linear model
    
    # A, numpy 2D matrix, constraint matrix
    # rhs, numpy 1D array, right hand side of constraint equation, e.g. Ax <= rhs
    # sense, sense of the inequalty, e.g. Ax <= rhs.  each entry cooresponds to a single row
    # lb, 1D np.array, lower bound of decision variables
    # ub, 1D np.array, upper bound of decision variables
    # vtype, list of varibles types (INTEGER, CONTINOUS, etc.)
    # solution, container for storing the solution                       

    # define the model
    model = Model()
    
    # configure for maximize (default is minimize)
    model.modelSense = -1
    
    num_rows = A.shape[0]
    num_cols = A.shape[1]

    # Add variables to model
    for j in range(num_cols):
      model.addVar(lb=lb[j], ub=ub[j], vtype=vtype[j])
    
    # update the model
    model.update()
    vars = model.getVars()
    
    # build constrainst
    # one constraint per row
    for i in range(num_rows):
        expr = LinExpr()
        
        # add up the variables from each column  
        # linear terms
        expr += quicksum(A[i,j]*vars[j] for j in range(num_cols))
                
        #add the constraint 
        model.addConstr(expr, sense[i], rhs[i])
    
    # Populate objective function
    if Q is not None :
        obj = QuadExpr()
    
#        # quadratic terms
#        obj += quicksum(Q[i,j]*vars[i]*vars[j]
#                        for i in range(num_cols) for j in range(num_cols))

        obj += quicksum(Q.flat[k]*vars[k / num_cols]*vars[k % num_cols] for k in range(num_cols*num_cols))    

#        # quadratic terms
#        obj += quicksum(Q.flat[k]*vars[k / num_cols]*vars[k % num_cols] for k in range(num_rows*num_cols))    

        # linear terms
        obj += quicksum(c[i,0]*vars[i] for i in range(num_cols))
        
    else :
        obj = LinExpr()
        
        # linear terms
        obj += quicksum(c[i,0]*vars[i] for i in range(num_cols))


    # lock in the objective function and update the model
    model.setObjective(obj)
    model.update()
    
    return(model)   
 
# %%    
def extract_solution(model) :

    decision_variables = model.getVars()
    Np = len(decision_variables)
    
    solution = np.matrix(np.zeros(shape = [Np,1] ))    
    for i in range(Np):
        solution[i,0] = decision_variables[i].x

    return solution          
    
    
# %% Build probability and weight data frames

def build_probability_matrices() :
    import os
    import pandas as pd

    print('\nBuilding Probability and Weight\'s matrices')    
    # change working directory
    os.chdir('/Users/rkelley/Development/rovi/gurobi/reach_study')
    print '\nchanged working directory to', os.getcwd()

    # Build probability matrix
    daypart_viewership_clean = pd.read_pickle('../data/daypart_viewership.pkl')
    
    ## build a probability and weight data frame
    
    Program_List = daypart_viewership_clean.tuple_id.unique()
    Persons_List = daypart_viewership_clean.person_id.unique()
    Persons_List.sort()
    
    # number of unique individuals
    Ni = Persons_List.shape[0]
    
    # number of unique programs
    Np = Program_List.shape[0]
    
    
    # initalize matrix to all zero
    prob_df = pd.DataFrame(np.zeros(shape=[Ni,Np]))
    prob_df.index = Persons_List[0:Ni]
    prob_df.index.name = 'Person_id'
    
    prob_df.columns = Program_List[0:Np]
    
    # the same for the weights
    weight_df = pd.DataFrame(np.zeros(shape=[Ni,Np]))
    weight_df.index = Persons_List[0:Ni]
    weight_df.index.name = 'Person_id'
    
    weight_df.columns = Program_List[0:Np]
    weight_df.head()
    
    
    for row in daypart_viewership_clean.index :
        person = daypart_viewership_clean.person_id[row]
        program = daypart_viewership_clean.tuple_id[row]
    
        weight = daypart_viewership_clean.weight[row]
        probability = daypart_viewership_clean.probability[row]
    
        weight_df.loc[person][program] = weight
        prob_df.loc[person][program] = probability
    
        if row % 50000 == 0 :
            print(row)
    
    print('\nsaving to ../data\n')
    
    prob_df.to_pickle('../data/Probabilities_dataframe.pkl')
    weight_df.to_pickle('../data/Weights_dataframe.pkl')
    
    return( [prob_df, weight_df] )

 

# %% load the matricies    
def load_probability_matrices() :
    
    import os
    import pandas as pd
    
    print('\nLoad Probability and Weight\'s matrices')    
    
    os.chdir('/Users/rkelley/Development/rovi/gurobi/reach_study')
    print '\nchanged working directory to', os.getcwd()

    prob_df   = pd.read_pickle('../data/Probabilities_dataframe.pkl')
    weight_df = pd.read_pickle('../data/Weights_dataframe.pkl')  
    
    return([prob_df, weight_df])