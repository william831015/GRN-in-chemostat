def dilute(molecule_diluted,molecules_0,DR=0.2):  #input Object want to dilute and where the parameters is stored 
    
    molecules_0[molecule_diluted.idx] *= (1-DR)

def replenish(molecule_replenished, molecules_0, DR=0.2) : #input Object want to replenish and where the parameters is stored 

    molecules_0[molecule_replenished.idx] += molecule_replenished.lc * DR


def dilute_out_species(molecules_diluted,molecules_0,DR=0.2):  #dilute a list of molecules
    for molecule in (molecules_diluted):
        dilute(molecule,molecules_0,DR)

def replenish_species(molecules, molecules_0, DR=0.2) :   #replenish a list of molecules

   for i in (molecules):
        replenish(molecules_0[i.idx],DR)


def run(model,t,molecules_0,molecules_list,dilution_list,replenish_list,result_all):  

    for n in range (Switch_cycle[0],Switch_cycle[1]):
        #define time
        t_start= n*15
        t_end = (n+1)*15
        t= np.linspace(t_start,t_end,2)

        #solve equation and save result
        result = scipy.integrate.odeint(Repressor_model, molecules_0, t, args=parameters_list)
        result_all = np.append(result_all,result[1])
        
        #update parameter
        molecules_0 = result.transpose()[:,-1]
        
        #dilution 
        ###replenish 
        molecules_0[R.idx] = R.ic*DR+(1-DR)*molecules_0[R.idx]
        

        ###diute out
        dilute_out_species((T7_RNA,GFP,GFP_RNA,Repressor,Repressor_RNA),molecules_0)


    return result_all, molecules_0