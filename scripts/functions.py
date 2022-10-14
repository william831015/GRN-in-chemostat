import numpy as np
import scipy.integrate
import scipy.optimize

def dilute(molecule_diluted,molecules_0,DR=0.2):  #input Object want to dilute and where the parameters is stored 
    
    molecules_0[molecule_diluted.idx] *= (1-DR)

def replenish(molecule_replenished, molecules_0, DR=0.2) : #input Object want to replenish and where the parameters is stored 

    molecules_0[molecule_replenished.idx] += molecule_replenished.lc * DR


def dilute_species(molecules_diluted,molecules_0,DR=0.2):  #dilute a list of molecules
    for molecule in (molecules_diluted):
        dilute(molecule,molecules_0,DR)

def replenish_species(molecules_replenished, molecules_0, DR=0.2) :   #replenish a list of molecules

   for molecule in (molecules_replenished):
        replenish(molecule,molecules_0,DR)



def run_model(model,t,parameters_list,molecules_0,dilute_list,replenish_list,result_all):  
    start_cycle,end_cycle = np.array(t)*4
    for n in range (start_cycle,end_cycle):
        #define time
        t_start= n*15
        t_end = (n+1)*15
        t= np.linspace(t_start,t_end,2)

        #solve equation and save result
        result = scipy.integrate.odeint(model, molecules_0, t, args=parameters_list)
        result_all = np.append(result_all,result[1])
        
        #update parameter
        molecules_0 = result.transpose()[:,-1]
        
        #dilution 
        ###diute out
        dilute_species((dilute_list),molecules_0)
        
        ###replenish 
        replenish_species((replenish_list),molecules_0)
    return result_all,molecules_0
    
