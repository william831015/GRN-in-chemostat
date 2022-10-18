#class of different species
# ic : initial concentration
# lc : loading concentration
# txr : transcription rate 
# tlr : translation rate 
# dr : degradation rate

class DNA:
    ic = 0
    lc = 0
    dr = 0
    result =[]
    pc_result = []
    nc_result = []


class RNA:
    txr = 1
    ic = 0
    lc = 0
    dr = 0
    result =[]
    pc_result = []
    nc_result = []


class Protein: 
    tlr = 1
    ic = 0
    lc = 0
    dr = 0
    result =[]
    pc_result = []
    nc_result = []

class Resource:
    ic = 100
    lc = 100
    result =[]
    pc_result = []
    nc_result = []
