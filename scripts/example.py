import scripts.models
import scripts.species
import scripts.functions

#######Define species #######

#T7 RNAP
T7_DNA = DNA()
T7_RNA = RNA()
T7 = Protein()

#Repressor
Repressor_DNA = DNA()
Repressor_RNA = RNA()
Repressor = Protein()

#GFP
GFP_DNA = DNA()
GFP_RNA = RNA()
GFP = Protein()

#Resource
R      = Resource()

molecules_list = [T7_DNA,T7_RNA,T7,Repressor_DNA,Repressor_RNA,Repressor,GFP_DNA,GFP_RNA,GFP,R]

########define parameters #######




########define models #######

def resource_dependent_repressor_model(molecules,t,alpha,beta,K,KTL,KTX,lamb1,lamb2):
    R,dT,dG,mT,mG,pT,pG,TX,TL = molecules
    return np.array(
        [
            -alpha*TX*dT/(dT+dG+KTX)*pT*R/(R+K)-alpha*TX*dG/(dT+dG+KTX)*pT*R/(R+K)-beta*TL*mT/(mT+mG+KTL)*R/(R+K)-beta*TL*mG/(mT+mG+KTL)*R/(R+K),
            0,
            0,
            alpha*TX*pT*R/(R+K),
            alpha*TX*pT*R/(R+K),
            beta*TL*(R/(R+K)),
            beta*TL*R/(R+K),
            -lamb1*TX,
            -lamb1*TL
        ]
    )
    


########solve models #######