def dilution(molecule, dilution_ratio):

    molecule_diluted = molecule*(1-dilution_ratio)

    return molecule_diluted



def replenish(molecule_replenished, dilution_ratio, replenish_conc) : 

    molecule_replenished = replenish_conc * dilution_ratio

    return molecule_replenished;



