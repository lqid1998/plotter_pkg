from pymatgen.io.vasp.outputs import Vasprun, Eigenval
from pymatgen.electronic_structure.plotter  import DosPlotter
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


from pymatgen.electronic_structure.cohp import Cohp
from pymatgen.io.lobster.outputs import Cohpcar,Doscar

def lob_d_orbs_dos(path, s, exchange = False):
    
    doscar= Doscar(path + 'DOSCAR.lobster',path + 'POSCAR')
    cdos= doscar.completedos
    structure=cdos.structure
    pdos = cdos.pdos

    new_orb_names = ["$d_{xy}$", "$d_{xz}$", "$d_{yz}$", "$d_{x^{2}-y^{2}}$", "$d_{z^{2}}$"]
    old_orb_names=['d_xy','d_xz','d_yz','d_x^2-y^2','d_z^2']
    orb_list = list(pdos[structure[s]].keys()) # see https://pymatgen.org/pymatgen.electronic_structure.core.html
    d_dos = {}
    title = "d orbitals dos %s" % ( structure[s].species)
    #testvar=Orbital(orb_indices[0])
    
    for i in range(len(orb_list)):
        if not ('p' in orb_list[i]):
            orb_in_d_dos_flag= True

            for j in range(len(new_orb_names)):
                if old_orb_names[j] in orb_list[i]:
                    orb_name_in_latex = orb_list[i].replace(old_orb_names[j],new_orb_names[j])
                    d_dos[orb_name_in_latex] = cdos.get_site_orbital_dos(structure[s], orb_list[i])
                    orb_in_d_dos_flag=False

            if orb_in_d_dos_flag :
                d_dos[orb_list[i]] = cdos.get_site_orbital_dos(structure[s], orb_list[i])
                orb_in_d_dos_flag=False
    
    plotter = DosPlotter(sigma=0.1)
    plotter.add_dos_dict(d_dos)
    
    new_plotter = plotter.get_plot(xlim=[-10, 10], ylim=[-2, 2])
    new_plotter.title(title, size=30)
    
    invert_plotter = invert_axes(new_plotter)
    if exchange:
        invert_plotter = exchange_alpha_beta(invert_plotter)
    creat_dir(path,'lob_dos_fig')
    invert_plotter.savefig('%s/lob_dos_fig/dorb_dos_%s%s.png' %(path,s+1,structure[s].species), format='png', pad_inches=1, bbox_inches='tight')


    return invert_plotter





    
def main():
    # a list of paths to the output folder of the calculations 
    # containing the vasprun.xml and EIGENVAL files
    # makes it easier to plot multiple files for comparison



    # ignores the warnings to make it nicer
    warnings.filterwarnings("ignore")


    
if __name__ == '__main__':
    main()

