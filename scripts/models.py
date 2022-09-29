#ODEs of different models

#The simplest model of resource dependent model
def resource_dependent_model(molecules,t,alpha,beta,K,KTL,KTX,lamb1,lamb2):
    R,dT,dG,mT,mG,pT,pG,TX,TL = molecules
    return np.array(
        [
            -alpha*TX*dT*pT*R/(R+K)-alpha*TX*dG*pT*R/(R+K)-beta*TL*mT*R/(R+K)-beta*TL*mG*R/(R+K),
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
    
#resource dependent model that include saturation effect 
def resource_dependent_saturation_model(molecules,t,alpha,beta,K,KTL,KTX,lamb1,lamb2):
    R,dT,dG,mT,mG,pT,pG,TX,TL = molecules
    return np.array(
        [
            -alpha*TX*dT/(dT+dG+KTX)*pT*R/(R+K)-alpha*TX*dG/(dT+dG+KTX)*pT*R/(R+K)-beta*TL*mT/(mT+mG+KTL)*R/(R+K)-beta*TL*mG/(mT+mG+KTL)*R/(R+K),
            0,
            0,
            alpha*TX*dT/(dT+dG+KTX)*pT*R/(R+K),
            alpha*TX*dG/(dT+dG+KTX)*pT*R/(R+K),
            beta*TL*(mT/(mT+mG+KTL))*(R/(R+K)),
            beta*TL*mG/(mT+mG+KTL)*R/(R+K),
            -lamb1*TX,
            -lamb1*TL
        ]
    )
    