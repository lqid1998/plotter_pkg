from pymatgen.io.vasp.outputs import Eigenval
import warnings
import os

    
#DosPlotter.get_plot_ax = new_get_plot.get_plot_ax_func #Add this new method to DosPlotter


def band_gap(path): 
    gap = Eigenval(path+'/dos/' + 'EIGENVAL', occu_tol=1e-5).eigenvalue_band_properties[0]

    print("Band gap (eV) =" + str(gap))

# invers the axes of the graph and the xlim and ylim
# feel free to comment it out from the other functions

def invert_axes(plotter):
        
    gca = plotter.gca()
    
    for line in gca.lines:
        x = line.get_xdata()
        y = line.get_ydata()
        line.set_xdata(y)
        line.set_ydata(x)
    
    plotter.xlabel("Density of states")
    plotter.ylabel("Energies (eV)")
    
    xlim = plotter.ylim()
    ylim = plotter.xlim()
    
    plotter.xlim(xlim)
    plotter.ylim(ylim)
    #plotter.show()

    return plotter

def exchange_alpha_beta(plotter):
        
    gca = plotter.gca()
    
    for line in gca.lines:
        x = line.get_xdata()
        for i in range(len(x)):
            x[i]=-x[i]
        line.set_xdata(x)
    #plotter.show()

    return plotter

def creat_dir(path,dirname):
    if not os.path.exists('%s/%s' %(path,dirname)):
            os.mkdir('%s/%s' %(path,dirname))
    
def main():
    # a list of paths to the output folder of the calculations 
    # containing the vasprun.xml and EIGENVAL files
    # makes it easier to plot multiple files for comparison

    # ignores the warnings to make it nicer
    warnings.filterwarnings("ignore")


    
if __name__ == '__main__':
    main()

