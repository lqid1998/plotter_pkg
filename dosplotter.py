from pymatgen.electronic_structure.core import Orbital, OrbitalType
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.electronic_structure.plotter  import DosPlotter
import matplotlib.pyplot as plt
from pymatgen.util.plotting import pretty_plot
from plotter_pkg.plotter_main import creat_dir
from pymatgen.electronic_structure.core import Spin
import matplotlib.pyplot as plt
import numpy as np

def get_plot_ax_func(self, xlim=None, ylim=None,ax=None, color_list=None, color=None,element=None):
    """
    Plot the DOS in to an axes object of matplotlib. This method will not plot out the graph.
    
    I got this from the source code of .get_plot() method of the DosPlotter class in pymatgen.electronic_structure.plotter module
    I replace the plt. functions of matplotlib by the class method of axes. The main changes to the source code were shown in comments.

    Args:
        xlim: Specifies the x-axis limits. Set to None for automatic
            determination.
        ylim: Specifies the y-axis limits.
        ax: Specifies the target axes of matplotlib plot
    """

    if ax is None:
        ax=plt.gca()

    if color_list is None:
        ncolors = max(3, len(self._doses))
        ncolors = min(9, ncolors)
        import palettable

        # pylint: disable=E1101
        colors = palettable.colorbrewer.qualitative.Set1_9.mpl_colors
        color_list_flag=False
    else:
        color_list_flag=True

    if color is None:
        color_flag = False
    else:
        color_flag = True

    if element is not None:
        element_label = True
    else:
        element_label = False

    y = None
    alldensities = []
    allenergies = []

    
    #plt = pretty_plot(12, 8) This is moved to d_orbs_subdos.


    # Note that this complicated processing of energies is to allow for
    # stacked plots in matplotlib.
    for key, dos in self._doses.items():
        energies = dos["energies"]
        densities = dos["densities"]
        if not y:
            y = {
                Spin.up: np.zeros(energies.shape),
                Spin.down: np.zeros(energies.shape),
            }
        newdens = {}
        for spin in [Spin.up, Spin.down]:
            if spin in densities:
                if self.stack:
                    y[spin] += densities[spin]
                    newdens[spin] = y[spin].copy()
                else:

                    newdens[spin] = densities[spin]
        allenergies.append(energies)
        alldensities.append(newdens)

    keys = list(self._doses.keys())
    keys.reverse()
    alldensities.reverse()
    allenergies.reverse()
    allpts = []
    for i, key in enumerate(keys):
        x = []
        y = []
        for spin in [Spin.up, Spin.down]:
            if spin in alldensities[i]:
                densities = list(int(spin) * alldensities[i][spin])
                energies = list(allenergies[i])
                if spin == Spin.down:
                    energies.reverse()
                    densities.reverse()
                y.extend(energies)
                x.extend(densities)
        allpts.extend(list(zip(x, y)))      #some plt. were replaced by ax.
        
        if color_list_flag:
            color = color_list[key]
        elif color_flag is False:
            color=colors[i % ncolors]
        
        if element_label:
            label_text='%s-%s' %(element, str(key))
        else:
            label_text=str(key)


        
        if self.stack:
            ax.fill(x, y, color=color, label=label_text)        #use ax.fill instead of plt.fill
        else:
            ax.plot(x, y, color=color, label=label_text, linewidth=3)      #use ax.plot to plot the dos into the target axes, instead of plt.plot() directly 
        if not self.zero_at_efermi:
            ylim = ax.get_ylim()            #use ax.get_ylim() insteading of plt.ylim     
            ax.plot(                        #use ax.plot to plot the dos into the target axes, instead of plt.plot() directly
                [self._doses[key]["efermi"], self._doses[key]["efermi"]],
                ylim,
                color=colors[i % ncolors],
                linestyle="--",
                linewidth=2,
            )

    if xlim:
        ax.set_xlim(xlim)       #use ax.set_xlim() insteading of plt.xlim() 
    if ylim:
        ax.set_ylim(ylim)       #use ax.set_ylim() insteading of plt.ylim
    else:
        xlim = ax.get_xlim()    #use ax.get_xlim()
        relevanty = [p[1] for p in allpts if xlim[0] < p[0] < xlim[1]]
        ax.set_ylim((min(relevanty), max(relevanty)))

    if self.zero_at_efermi:
        ylim = ax.get_ylim()    #use ax.get_ylim()
        ax.plot([0, 0], ylim, "k--", linewidth=2)   #use ax.plot to plot the dos into the target ax

    ax.set_ylabel("Energies (eV)", fontsize=16)      #use ax.set_xlabel ax.set_ylabel instead of plt.xlabel plt.ylabel
    ax.set_xlabel("Density of states", fontsize=16)
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)
    ax.axhline(y=0, color="k", linestyle="--", linewidth=2)     #use ax.axhline instead of plt.axhline
    ax.legend(fontsize=30)         #use ax.legend() insteadof plt.legend

    #These lines can be replaced by ax.legend(fontsize=30)
    #leg = plt.gca().get_legend()
    #ltext = leg.get_texts()  # all the text.Text instance in the legend
    #plt.setp(ltext, fontsize=30)

    #This line is moved to d_orbs_subdos()
    #plt.tight_layout()

    #Modified from invert_axes(plotter) to exchange the x-axis and y-axis of the plot
    """
    for line in ax.lines:
        x = line.get_xdata()
        y = line.get_ydata()
        line.set_xdata(y)
        line.set_ydata(x)
    
    ax.set_ylabel("Energies (eV)")
    ax.set_xlabel("Density of states")
    
    xlim_invert = ylim
    ylim_invert = xlim
    
    ax.set_xlim(xlim_invert)
    ax.set_ylim(ylim_invert)
    """
    return ax
    
DosPlotter.get_plot_ax = get_plot_ax_func #Add this new method to DosPlotter


def set_color(orb_list):
    import palettable
    orb_index={"$d_{xy}$": 0, "$d_{xz}$": 1, "$d_{yz}$": 2, "$d_{x^{2}-y^{2}}$": 3, "$d_{z^{2}}$": 4\
        ,"$s$": 5,"$p_{x}$": 6,"$p_{y}$": 7,"$p_{z}$": 8,}
    colors = palettable.colorbrewer.qualitative.Set1_9.mpl_colors
    testvar=palettable.colorbrewer.qualitative.Dark2_8.mpl_colors
    color_list={}
    # pylint: disable=E1101
    for i in range(len(orb_list)):
        orb_name=orb_list[i]
        color_list[orb_name]=colors[orb_index[orb_name]]
    
    color_list['total']=(0.7,0.7,0.7)
    return color_list

def dos(path, s):
   
    v = Vasprun(path + 'vasprun.xml')
    cdos = v.complete_dos
    structure = cdos.structure.sites

    fig = pretty_plot(12, 8)
    title = "dos %s" % ( structure[s].species)
    fig.suptitle(title, size=30)
    #testvar=structure[0][s]
    spd_dos = cdos.get_site_spd_dos(structure[s])
    
    plotter = DosPlotter(sigma=0.1)
    plotter.add_dos_dict(spd_dos)

    ax=fig.subplot(111)
    plotter = DosPlotter(sigma=0.1)
    plotter.add_dos_dict(spd_dos)
    plotter.get_plot_ax(xlim=[-2, 2], ylim=[-10, 10], ax=ax)

    #This line is from DosPlotter.get_plot()ax.legend()                        
    plt.tight_layout()
    creat_dir(path,'dos_fig')
    plt.savefig('%s/dos_fig/dos_%s%s.png' %(path,s+1,structure[s].species), format='png', pad_inches=1, bbox_inches='tight')

    return plt
# MUST BE FROM A SINGLE POINT CALCULATION
# f is the calculation's folder position in the path array
# s is the atom you want to print (starts with 0)
# d is the d-orbital label

def d_orbs_dos(path, s):
   
    v = Vasprun(path + 'vasprun.xml')
    cdos = v.complete_dos
    structure = cdos.structure.sites

    orb_names = ["$d_{xy}$", "$d_{xz}$", "$d_{yz}$", "$d_{x^{2}-y^{2}}$", "$d_{z^{2}}$"]
    orb_indices = [4, 7, 5, 8, 6] # see https://pymatgen.org/pymatgen.electronic_structure.core.html
    d_dos = {}

    fig=pretty_plot(12,8)
    title = "d orbitals dos %s" % ( structure[s].species)
    fig.suptitle(title,size=30)

    color_list=set_color(orb_names)
    for i in range(5):
        d_dos[orb_names[i]] = cdos.get_site_orbital_dos(structure[s], Orbital(orb_indices[i]))
    
    ax=fig.subplot(111)
    plotter = DosPlotter(sigma=0.1)
    plotter.add_dos_dict(d_dos)
    
    plotter.get_plot_ax(xlim=[-2, 2], ylim=[-10, 10],ax=ax, color_list=color_list)
    
    #invert_plotter = invert_axes(new_plotter)
    plt.tight_layout()
    creat_dir(path,'dos_fig')
    plt.savefig('%s/dos_fig/dorb_dos_%s%s.png' %(path,s+1,structure[s].species), format='png', pad_inches=1, bbox_inches='tight')
    return plt

def p_orbs_dos(path, s):
   
    v = Vasprun(path + 'vasprun.xml')
    cdos = v.complete_dos
    structure = cdos.structure.sites
    
    orb_names = ["$p_{x}$", "$p_{y}$", "$p_{z}$"]
    orb_indices = [3,1,2] # see https://pymatgen.org/pymatgen.electronic_structure.core.html
    d_dos = {}
    fig=pretty_plot(12,8)
    title = "p orbitals dos %s" % ( structure[s].species)
    fig.suptitle(title,size=30)
    color_list=set_color(orb_names)
    for i in range(3):
        d_dos[orb_names[i]] = cdos.get_site_orbital_dos(structure[s], Orbital(orb_indices[i]))
    
    ax=fig.subplot(111)
    plotter = DosPlotter(sigma=0.1)
    plotter.add_dos_dict(d_dos)
    
    plotter.get_plot_ax(xlim=[-2, 2], ylim=[-10, 10],ax=ax,color_list=color_list)
    plt.tight_layout()

    creat_dir(path,'dos_fig')
    plt.savefig('%s/dos_fig/porb_dos_%s%s.png' %(path,s+1,structure[s].species), format='png', pad_inches=1, bbox_inches='tight')
    return plt




def d_orbs_pdos(path, s,d):
    
    v = Vasprun(path + 'vasprun.xml')
    cdos = v.complete_dos
    structure = cdos.structure.sites
    
    orb_names = ["$d_{xy}$", "$d_{xz}$", "$d_{yz}$", "$d_{x^{2}-y^{2}}$", "$d_{z^{2}}$"]
    orb_indices = [4, 7, 5, 8, 6] # see https://pymatgen.org/pymatgen.electronic_structure.core.html
    d_dos = {}
    fig=pretty_plot(12,8)
    title = "d orbitals dos (" + str(s) + ") %s %s " %(path, structure[s].species)
    fig.suptitle(title,size=30)
    color_list=set_color(orb_names)


    d_dos[orb_names[d]] = cdos.get_site_orbital_dos(structure[s], Orbital(orb_indices[d]))
    
    ax=fig.subplot(111)
    plotter = DosPlotter(sigma=0.05)
    plotter.add_dos_dict(d_dos)
    
    plotter.get_plot_ax(xlim=[-2, 2], ylim=[-10, 10],ax=ax,color_list=color_list)
    plt.tight_layout()
    
    creat_dir(path,'dos_fig')
    plt.savefig('%s/dos_fig/dorb_pdos_%s%s_%s.png' %(path,s+1,structure[s].species,d), format='png', pad_inches=1, bbox_inches='tight')
    return plt






def d_orbs_subdos(path, s):

    
    v = Vasprun(path + 'vasprun.xml')
    cdos = v.complete_dos
    structure = cdos.structure
    
    
    fig = pretty_plot(36, 16)      #This line is from DosPlotter.get.plot(), but tripled in the x-axis and doubled in the y-axis, in order to arraying the subplots as 2*3.

    title = "d orbitals dos %s %s" % (path , structure[s].species)
    fig.suptitle(title, size=30)

    orb_names = ["$d_{xy}$", "$d_{xz}$", "$d_{yz}$", "$d_{x^{2}-y^{2}}$", "$d_{z^{2}}$"]
    orb_indices = [4, 7, 5, 8, 6] # see https://pymatgen.org/pymatgen.electronic_structure.core.html
    d_dos = {}
    color_list=set_color(orb_names)
    #plot each d-orbital projected dos into subplot 1-5
    total_dos = cdos.get_site_dos(structure[s])
    total_plotter = DosPlotter(sigma=0.1,stack=True)
    total_plotter.add_dos("total", total_dos)

    for i in range(5):      
        ax = fig.subplot(2,3,(i+1))
        d_pdos=cdos.get_site_orbital_dos(structure[s], Orbital(orb_indices[i]))
        d_dos[orb_names[i]]= d_pdos
        total_plotter.get_plot_ax(xlim=[-2,2],ylim=[-10,10],ax=ax,color_list=color_list,element='Ru')
        plotter = DosPlotter(sigma=0.1)
        plotter.add_dos(orb_names[i],d_pdos)
        plotter.get_plot_ax(xlim=[-2, 2], ylim=[-10, 10], ax=ax, color_list=color_list,element='Ru')
    
    #plot all d-orbitals projected dos together into subplot 6
    ax=fig.subplot(2,3,6)
    plotter = DosPlotter(sigma=0.1)
    plotter.add_dos_dict(d_dos)
    plotter.get_plot_ax(xlim=[-2, 2], ylim=[-10, 10], ax=ax,color_list=color_list)


    #This line is from DosPlotter.get_plot()ax.legend()                        
    plt.tight_layout()
    plt.savefig('%s/dos_fig/dorb_subdos_%s%s.png' %(path,s+1,structure[s].species), format='png', pad_inches=1, bbox_inches='tight')
    return plt


def d_orbs_dos_with_total(path,s,d,element=None):
    v=Vasprun(path+"vasprun.xml")
    cdos = v.complete_dos
    structure = cdos.structure.sites

    orb_names = ["$d_{xy}$", "$d_{xz}$", "$d_{yz}$", "$d_{x^{2}-y^{2}}$", "$d_{z^{2}}$"]
    orb_indices = [4, 7, 5, 8, 6] # see https://pymatgen.org/pymatgen.electronic_structure.core.html
    d_dos = {}
    fig = pretty_plot(12,8)
    title = "dos %s" % ( structure[s].species)
    fig.suptitle(title, size=30)

    total_dos = cdos.get_site_dos(structure[s])
    spd_dos = cdos.get_site_spd_dos(structure[s])
    total_plotter = DosPlotter(sigma=0.1,stack=True)
    total_plotter.add_dos("total", total_dos)

    pdos_plotter = DosPlotter(sigma=0.1)
    d_dos[orb_names[d]] = cdos.get_site_orbital_dos(structure[s], Orbital(orb_indices[d]))
    pdos_plotter.add_dos_dict(d_dos)

    ax=fig.subplot(111)
    color_list=set_color(orb_names)
    total_plotter.get_plot_ax(xlim=[-2,2],ylim=[-10,10],ax=ax,color_list=color_list,element='Ru')
    pdos_plotter.get_plot_ax(xlim=[-2,2],ylim=[-10,10],ax=ax,color_list=color_list, element='Ru')
    
    plt.tight_layout()
    creat_dir(path,'dos_fig')
    plt.savefig('%s/dos_fig/dorb_pdos_with_total_%s%s_%s.png' %(path,s+1,structure[s].species,d), format='png', pad_inches=1, bbox_inches='tight')
    return plt

def element_dos_with_total(path,s1,s2):
    v=Vasprun(path+"vasprun.xml")
    cdos = v.complete_dos
    structure = cdos.structure.sites

    fig = pretty_plot(12,8)
    title = "dos %s" % ( structure[s1].species)
    fig.suptitle(title, size=30)

    total_dos = cdos.get_site_dos(structure[s2])

    total_plotter = DosPlotter(sigma=0.1,stack=True)
    total_plotter.add_dos("total", total_dos)

    sdos_plotter = DosPlotter(sigma=0.1)
    pdos_plotter = DosPlotter(sigma=0.1)
    spd_dos = cdos.get_site_spd_dos(structure[s1])

    #pdos_plotter.add_dos('N-s',spd_dos[OrbitalType(0)])
    sdos_plotter.add_dos('H-s',spd_dos[OrbitalType(0)])
    pdos_plotter.add_dos('N-p',spd_dos[OrbitalType(1)])
    ax=fig.subplot(111)
    total_plotter.get_plot_ax(xlim=[-2,2],ylim=[-10,10],ax=ax, color=(0.7,0.7,0.7),element='Ru')
    sdos_plotter.get_plot_ax(xlim=[-2,2],ylim=[-10,10],ax=ax,color=(0.6509803921568628, 0.33725490196078434, 0.1568627450980392))
    #pdos_plotter.get_plot_ax(xlim=[-2,2],ylim=[-10,10],ax=ax,color=(0.9058823529411765, 0.1607843137254902, 0.5411764705882353))
    
    plt.tight_layout()
    creat_dir(path,'dos_fig')
    plt.savefig('%s/dos_fig/%s%s_with_total.png' %(path,s1+1,structure[s1].species), format='png', pad_inches=1, bbox_inches='tight')
    return plt

