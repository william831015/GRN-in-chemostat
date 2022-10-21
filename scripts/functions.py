import numpy as np
import scipy.integrate
import scipy.optimize
import bokeh.plotting
from bokeh.plotting import figure, output_file, show
import bokeh.io
from bokeh.models import Span

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
    

def plot_result(molecule):
    t = np.linspace(0, 15*(len(molecule)-1), len(molecule))
    p = bokeh.plotting.figure(
        plot_width=800,
        plot_height=400,
        x_axis_label="t",
        y_axis_type="linear",
    )

    colors = bokeh.palettes.d3["Category10"][3]

    # Populate glyphs
    p.line(
        t/60, molecule, line_width=2, color=colors[0]
    )
    vline1 = Span(location=4, dimension='height', line_color='black', line_width=1,line_dash='dashed')
    vline2 = Span(location=16, dimension='height', line_color='black', line_width=1,line_dash='dashed')
    p.add_layout(vline1)
    p.add_layout(vline2)
    show(p)


def plot_result_two_state(molecule):
    t = np.linspace(0, 15*(len(molecule)-1), len(molecule))
    p = bokeh.plotting.figure(
        plot_width=800,
        plot_height=400,
        x_axis_label="t",
        y_axis_type="linear",
    )

    colors = bokeh.palettes.d3["Category10"][3]

    # Populate glyphs
    p.line(
        t/60, molecule, line_width=2, color=colors[0]
    )
    vline1 = Span(location=4, dimension='height', line_color='black', line_width=1,line_dash='dashed')
    #vline2 = Span(location=16, dimension='height', line_color='black', line_width=1,line_dash='dashed')
    p.add_layout(vline1)
    #p.add_layout(vline2)
    show(p)