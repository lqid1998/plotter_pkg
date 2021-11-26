from pymatgen.electronic_structure.core import Orbital, OrbitalType
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.electronic_structure.plotter  import DosPlotter
import matplotlib.pyplot as plt
from pymatgen.io.lobster.outputs import Doscar
from pymatgen.util.plotting import pretty_plot
from plotter_pkg.plotter_main import creat_dir
from pymatgen.electronic_structure.core import Spin
import matplotlib.pyplot as plt
import numpy as np

def get_cdos_vasp(path):
    try:
        v = Vasprun(path + '/vasprun.xml')
        cdos = v.complete_dos
        structure = cdos.structure.sites
        return cdos, structure
    except Exception as err:
        print("an error occured when reading in vasprun.xml")
        print(err)

def get_cdos_lobster(path):
    try:
        doscar = Doscar(path + '/DOSCAR.lobster',path + '/POSCAR')
        cdos = doscar.completedos
        structure = cdos.structure.sites
        return cdos, structure
    except Exception as err:
        print("an error occured when reading in DOSCAR.lobster")
        print(err)

def get_plot_ax_func(self, xlim=None, ylim=None,ax=None, color_list=None, color=None,element=None, label=None):
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
        
        if label:
            label_text= '%s' %label
        elif element_label:
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

def full_color():
    import palettable
    orb_index = {"$d_{xy}$": 0, "$3d_{xy}$": 0, "$4d_{xy}$": 0, "$d_{xz}$": 1, "$3d_{xz}$": 1, "$4d_{xz}$": 1, \
        "$d_{yz}$": 2, "$3d_{yz}$": 2, "$4d_{yz}$": 2, "$d_{x^{2}-y^{2}}$": 3, "$3d_{x^{2}-y^{2}}$": 3, "$4d_{x^{2}-y^{2}}$": 3, \
            "$d_{z^{2}}$": 4, "$3d_{z^{2}}$": 4, "$4d_{z^{2}}$": 4, "$s$": 5, "$1s$": 5, "$2s$": 5, "$3s$": 5, "$4s$": 5, \
                "$p_{x}$": 6, "$1p_{x}$": 6, "$2p_{x}$": 6, "$3p_{x}$": 6, \
                    "$p_{y}$": 7, "$1p_{y}$": 7, "$2p_{y}$": 7, "$3p_{y}$": 7, \
                        "$p_{z}$": 8, "$1p_{z}$": 8, "$2p_{z}$": 8, "$3p_{z}$": 8}
    colors = palettable.colorbrewer.qualitative.Set1_9.mpl_colors
    color_list = {}
    for i, item in enumerate(orb_index.items()):
        color_list[item[0]] = colors[item[1]]
    color_list['total']=(0.7,0.7,0.7)
    return color_list

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

def spd_dos(cdos, structure, atom_label, xlim=None, ylim=None, save_fig=True, save_path=None, set_title=True, title=None):
    #testvar=structure[0][s]
    spd_dos = cdos.get_site_spd_dos(structure[atom_label])
    plotter = DosPlotter(sigma=0.1)
    plotter.add_dos_dict(spd_dos)

    if save_fig:
        fig = pretty_plot(12,8)
        ax=fig.subplot(111)

        if set_title:
            if not title:
                title = "dos %s" % ( structure[atom_label].species)
            fig.suptitle(title,size=30)

        if not xlim:
            xlim = [-2, 2]
        if not ylim:
            ylim = [-10, 10]
        
        plotter.get_plot_ax(xlim=xlim, ylim=ylim, ax=ax)
        plt.tight_layout()
        try:
            creat_dir(save_path,'dos_fig')
        except Exception as err:
            print("please input the path to save")
            print(err)
        try:
            plt.savefig('%s/dos_fig/dos_%s%s.png' %(save_path, atom_label + 1, structure[atom_label].species), format='png', \
                pad_inches=1, bbox_inches='tight')
        except Exception as err:
            print("an error occured when save the figure")
            print(err)
    return plotter
# MUST BE FROM A SINGLE POINT CALCULATION
# f is the calculation's folder position in the path array
# s is the atom you want to print (starts with 0)
# d is the d-orbital label
def total_dos(cdos, structure, atom_label, xlim=None, ylim=None, save_fig=True, save_path=None, set_title=True, title=None):
    t_dos = cdos.get_site_dos(structure[atom_label])
    plotter = DosPlotter(sigma=0.1)
    plotter.add_dos("total", t_dos)

    if save_fig:
        fig = pretty_plot(12,8)
        ax=fig.subplot(111)

        if set_title:
            if not title:
                title = "total dos %s" % ( structure[atom_label].species)
            fig.suptitle(title,size=30)

        if not xlim:
            xlim = [-2, 2]
        if not ylim:
            ylim = [-10, 10]
        
        plotter.get_plot_ax(xlim=xlim, ylim=ylim, ax=ax)
        plt.tight_layout()
        try:
            creat_dir(save_path,'dos_fig')
        except Exception as err:
            print("please input the path to save")
            print(err)
        try:
            plt.savefig('%s/dos_fig/torb_dos_%s%s.png' %(save_path, atom_label + 1, structure[atom_label].species), format='png', \
                pad_inches=1, bbox_inches='tight')
        except Exception as err:
            print("an error occured when save the figure")
            print(err)
    return plotter

def s_orbs_dos(cdos, structure, atom_label, xlim=None, ylim=None, save_fig=True, save_path=None, set_title=True, title=None):
    s_dos = cdos.get_site_orbital_dos(structure[atom_label], Orbital(0))
    plotter = DosPlotter(sigma=0.1)
    plotter.add_dos("s", s_dos)

    if save_fig:
        fig = pretty_plot(12,8)
        ax=fig.subplot(111)

        if set_title:
            if not title:
                title = "s orbitals dos %s" % ( structure[atom_label].species)
            fig.suptitle(title,size=30)

        if not xlim:
            xlim = [-2, 2]
        if not ylim:
            ylim = [-10, 10]
        
        plotter.get_plot_ax(xlim=xlim, ylim=ylim, ax=ax)
        plt.tight_layout()
        try:
            creat_dir(save_path,'dos_fig')
        except Exception as err:
            print("please input the path to save")
            print(err)
        try:
            plt.savefig('%s/dos_fig/sorb_dos_%s%s.png' %(save_path, atom_label + 1, structure[atom_label].species), format='png', \
                pad_inches=1, bbox_inches='tight')
        except Exception as err:
            print("an error occured when save the figure")
            print(err)
    return plotter

def p_orbs_dos(cdos, structure, atom_label, xlim=None, ylim=None, save_fig=True, save_path=None, set_title=True, title=None):
    orbit_tuple = (("$p_{x}$", 3,), ("$p_{y}$", 1,), ("$p_{z}$", 2,),)
    # see https://pymatgen.org/pymatgen.electronic_structure.core.html
    p_dos = {}
    for i in range(len(orbit_tuple)):
        p_dos[orbit_tuple[i][0]] = cdos.get_site_orbital_dos(structure[atom_label], Orbital(orbit_tuple[i][1]))
    plotter = DosPlotter(sigma=0.1)
    plotter.add_dos_dict(p_dos)

    if save_fig:
        color_list=full_color()
        fig = pretty_plot(12,8)
        ax=fig.subplot(111)

        if set_title:
            if not title:
                title = "p orbitals dos %s" % ( structure[atom_label].species)
            fig.suptitle(title,size=30)

        if not xlim:
            xlim = [-2, 2]
        if not ylim:
            ylim = [-10, 10]
        
        plotter.get_plot_ax(xlim=xlim, ylim=ylim, ax=ax, color_list=color_list)
        plt.tight_layout()
        try:
            creat_dir(save_path,'dos_fig')
        except Exception as err:
            print("please input the path to save")
            print(err)
        try:
            plt.savefig('%s/dos_fig/porb_dos_%s%s.png' %(save_path, atom_label + 1, structure[atom_label].species), format='png', \
                pad_inches=1, bbox_inches='tight')
        except Exception as err:
            print("an error occured when save the figure")
            print(err)
    return plotter

def d_orbs_dos(cdos, structure, atom_label, xlim=None, ylim=None, save_fig=True, save_path=None, set_title=True, title=None):
    orbit_tuple = (("$d_{xy}$", 4), ("$d_{xz}$", 7), ("$d_{yz}$", 5), ("$d_{x^{2}-y^{2}}$", 8), ("$d_{z^{2}}$", 6),)
    # see https://pymatgen.org/pymatgen.electronic_structure.core.html
    d_dos = {}
    for i in range(len(orbit_tuple)):
        d_dos[orbit_tuple[i][0]] = cdos.get_site_orbital_dos(structure[atom_label], Orbital(orbit_tuple[i][1]))
    plotter = DosPlotter(sigma=0.1)
    plotter.add_dos_dict(d_dos)
    
    if save_fig:
        color_list=full_color()
        fig = pretty_plot(12,8)
        ax=fig.subplot(111)

        if set_title:
            if not title:
                title = "d orbitals dos %s" % ( structure[atom_label].species)
            fig.suptitle(title,size=30)

        if not xlim:
            xlim = [-2, 2]
        if not ylim:
            ylim = [-10, 10]
        
        plotter.get_plot_ax(xlim=xlim, ylim=ylim, ax=ax, color_list=color_list)
        plt.tight_layout()
        try:
            creat_dir(save_path,'dos_fig')
        except Exception as err:
            print("please input the path to save")
            print(err)
        try:
            plt.savefig('%s/dos_fig/dorb_dos_%s%s.png' %(save_path, atom_label + 1, structure[atom_label].species), format='png', \
                pad_inches=1, bbox_inches='tight')
        except Exception as err:
            print("an error occured when save the figure")
            print(err)
    return plotter

def d_orbs_pdos(cdos, structure, atom_label, orbital_label, xlim=None, ylim=None, save_fig=True, save_path=None, set_title=True, \
    title=None):
    orbit_tuple = (("$d_{xy}$", 4), ("$d_{xz}$", 7), ("$d_{yz}$", 5), ("$d_{x^{2}-y^{2}}$", 8), ("$d_{z^{2}}$", 6),)
    # see https://pymatgen.org/pymatgen.electronic_structure.core.html
    d_dos = cdos.get_site_orbital_dos(structure[atom_label], Orbital(orbit_tuple[orbital_label][1]))
    plotter = DosPlotter(sigma=0.05)
    plotter.add_dos(orbit_tuple[orbital_label][0], d_dos)
 
    if save_fig:
        color_list=full_color()
        fig = pretty_plot(12,8)
        ax=fig.subplot(111)

        if set_title:
            if not title:
                title = "d orbitals dos (" + str(atom_label) + ") %s %s " %(save_path, structure[atom_label].species)
            fig.suptitle(title,size=30)
        
        if not xlim:
            xlim = [-2, 2]
        if not ylim:
            ylim = [-10, 10]
        
        plotter.get_plot_ax(xlim=xlim, ylim=ylim, ax=ax, color_list=color_list)
        plt.tight_layout()
        try:
            creat_dir(save_path,'dos_fig')
        except Exception as err:
            print("please input the path to save")
            print(err)
        try:
            plt.savefig('%s/dos_fig/dorb_pdos_%s%s_%s.png' %(save_path,atom_label+1,structure[atom_label].species,orbital_label), \
            format='png', pad_inches=1, bbox_inches='tight')
        except Exception as err:
            print("an error occured when save the figure")
            print(err)
    return plotter

def d_orbs_subdos(path, s):

    
    v = Vasprun(path + '/vasprun.xml')
    cdos = v.complete_dos
    structure = cdos.structure
    
    
    fig = pretty_plot(36, 16)      #This line is from DosPlotter.get.plot(), but tripled in the x-axis and doubled in the y-axis, in order to arraying the subplots as 2*3.

    title = "d orbitals dos %s %s" % (path , structure[s].species)
    fig.suptitle(title, size=30)

    orb_names = ["$d_{xy}$", "$d_{xz}$", "$d_{yz}$", "$d_{x^{2}-y^{2}}$", "$d_{z^{2}}$"]
    orb_indices = [4, 7, 5, 8, 6] # see https://pymatgen.org/pymatgen.electronic_structure.core.html
    d_dos = {}
    color_list=full_color()
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
    v=Vasprun(path+"/vasprun.xml")
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
    color_list=full_color()
    total_plotter.get_plot_ax(xlim=[-2,2],ylim=[-10,10],ax=ax,color_list=color_list,element='Ru')
    pdos_plotter.get_plot_ax(xlim=[-2,2],ylim=[-10,10],ax=ax,color_list=color_list, element='Ru')
    
    plt.tight_layout()
    creat_dir(path,'dos_fig')
    plt.savefig('%s/dos_fig/dorb_pdos_with_total_%s%s_%s.png' %(path,s+1,structure[s].species,d), format='png', pad_inches=1, bbox_inches='tight')
    return plt

def element_dos_with_total(path, s1, s2):
    v=Vasprun(path+"/vasprun.xml")
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

def lob_spd_dos(cdos, structure, atom_label, xlim=None, ylim=None, save_fig=True, save_path=None, set_title=True, title=None):
    spd_dos = cdos.get_site_spd_dos(structure[atom_label])
    plotter = DosPlotter(sigma=0.1)
    plotter.add_dos_dict(spd_dos)

    if save_fig:
        fig = pretty_plot(12,8)
        ax=fig.subplot(111)

        if set_title:
            if not title:
                title = "lobster dos %s" % ( structure[atom_label].species)
            fig.suptitle(title,size=30)

        if not xlim:
            xlim = [-2, 2]
        if not ylim:
            ylim = [-10, 10]
        
        plotter.get_plot_ax(xlim=xlim, ylim=ylim, ax=ax)
        plt.tight_layout()
        try:
            creat_dir(save_path,'lob_dos_fig')
        except Exception as err:
            print("please input the path to save")
            print(err)
        try:
            plt.savefig('%s/lob_dos_fig/spd_dos_%s%s.png' %(save_path, atom_label + 1, structure[atom_label].species), format='png', \
                pad_inches=1, bbox_inches='tight')
        except Exception as err:
            print("an error occured when save the figure")
            print(err)
    return plotter

def lob_total_dos(cdos, structure, atom_label, xlim=None, ylim=None, save_fig=True, save_path=None, set_title=True, title=None):
    t_dos = cdos.get_site_dos(structure[atom_label])
    plotter = DosPlotter(sigma=0.1)
    plotter.add_dos("total", t_dos)

    if save_fig:
        fig = pretty_plot(12,8)
        ax=fig.subplot(111)

        if set_title:
            if not title:
                title = "lobster dos %s" % ( structure[atom_label].species)
            fig.suptitle(title,size=30)

        if not xlim:
            xlim = [-2, 2]
        if not ylim:
            ylim = [-10, 10]
        
        plotter.get_plot_ax(xlim=xlim, ylim=ylim, ax=ax)
        plt.tight_layout()
        try:
            creat_dir(save_path,'lob_dos_fig')
        except Exception as err:
            print("please input the path to save")
            print(err)
        try:
            plt.savefig('%s/lob_dos_fig/total_dos_%s%s.png' %(save_path, atom_label + 1, structure[atom_label].species), format='png', \
                pad_inches=1, bbox_inches='tight')
        except Exception as err:
            print("an error occured when save the figure")
            print(err)
    return plotter

def replace_latex_orb_list(orb_list):
    orb_names = {"s":("$s$", 0,), "p_x":("$p_{x}$", 3,), "p_y":("$p_{y}$", 1,), "p_z":("$p_{z}$", 2,), "d_xy":("$d_{xy}$", 4,), \
        "d_xz":("$d_{xz}$", 7,), "d_yz":("$d_{yz}$", 5,), "d_x^2-y^2":("$d_{x^{2}-y^{2}}$", 8,), "d_z^2":("$d_{z^{2}}$", 6,)}
    new_orb_list = []
    for i in range(len(orb_list)):
        old_orb_name = orb_list[i]
        n = old_orb_name[0]
        old_orb_label = old_orb_name[1:]
        new_orb_label_list = list(orb_names[old_orb_label][0])
        new_orb_label_list.insert(1, n)
        new_orb_name = "".join(new_orb_label_list)
        new_orb_list.append(new_orb_name)
    return new_orb_list

def lob_d_orbs_dos(cdos, structure, atom_label, xlim=None, ylim=None, save_fig=True, save_path=None, set_title=True, title=None):
    orb_list = list(cdos.pdos[structure[atom_label]].keys())
    orb_list_latex = replace_latex_orb_list(orb_list)
    d_orbs = ["d_xy", "d_xz", "d_yz", "d_x^2-y^2", "d_z^2"]
    d_dos = {}
    for i in range(len(orb_list)):
        if orb_list[i][1:] in d_orbs:
            d_dos[orb_list_latex[i]] = cdos.get_site_orbital_dos(structure[atom_label], orb_list[i])
    plotter = DosPlotter(sigma=0.1)
    plotter.add_dos_dict(d_dos)
    
    if save_fig:
        color_list=full_color()
        fig = pretty_plot(12,8)
        ax=fig.subplot(111)

        if set_title:
            if not title:
                title = "d orbitals dos %s" % ( structure[atom_label].species)
            fig.suptitle(title,size=30)

        if not xlim:
            xlim = [-2, 2]
        if not ylim:
            ylim = [-10, 10]
        
        plotter.get_plot_ax(xlim=xlim, ylim=ylim, ax=ax, color_list=color_list)
        plt.tight_layout()
        try:
            creat_dir(save_path,'lob_dos_fig')
        except Exception as err:
            print("please input the path to save")
            print(err)
        try:
            plt.savefig('%s/lob_dos_fig/dorb_dos_%s%s.png' %(save_path, atom_label + 1, structure[atom_label].species), format='png', \
                pad_inches=1, bbox_inches='tight')
        except Exception as err:
            print("an error occured when save the figure")
            print(err)
    return plotter

def lob_d_orbs_pdos(cdos, structure, atom_label, orbital_label, xlim=None, ylim=None, save_fig=True, save_path=None, set_title=True, title=None):
    d_orbs = ["d_xy", "d_xz", "d_yz", "d_x^2-y^2", "d_z^2"]
    orb_list = list(cdos.pdos[structure[atom_label]].keys())
    d_orbs_list = []
    for i in range(len(d_orbs)):
        for j in range(len(orb_list)):
            if orb_list[j][1:] == d_orbs[i]:
                d_orbs_list.append(orb_list[j])
    d_orbs_list_latex = replace_latex_orb_list(d_orbs_list)

    d_dos = cdos.get_site_orbital_dos(structure[atom_label], d_orbs_list[orbital_label])
    plotter = DosPlotter(sigma=0.1)
    plotter.add_dos(d_orbs_list_latex[orbital_label], d_dos)
    
    if save_fig:
        color_list=full_color()
        fig = pretty_plot(12,8)
        ax=fig.subplot(111)

        if set_title:
            if not title:
                title = "d orbitals dos (" + str(atom_label) + ") %s %s " %(save_path, structure[atom_label].species)
            fig.suptitle(title,size=30)

        if not xlim:
            xlim = [-2, 2]
        if not ylim:
            ylim = [-10, 10]
        
        plotter.get_plot_ax(xlim=xlim, ylim=ylim, ax=ax, color_list=color_list)
        plt.tight_layout()
        try:
            creat_dir(save_path,'lob_dos_fig')
        except Exception as err:
            print("please input the path to save")
            print(err)
        try:
            plt.savefig('%s/lob_dos_fig/dorb_pdos_%s%s_%s.png' %(save_path,atom_label+1,structure[atom_label].species,orbital_label), \
            format='png', pad_inches=1, bbox_inches='tight')
        except Exception as err:
            print("an error occured when save the figure")
            print(err)
    return plotter

def lob_p_orbs_dos(path, s):
    
    doscar = Doscar(path + '/DOSCAR.lobster',path + '/POSCAR')
    cdos = doscar.completedos
    structure = cdos.structure.sites
    pdos = cdos.pdos

    orb_names = ["$p_{x}$", "$p_{y}$", "$p_{z}$"]
    old_orb_names=['p_x','p_y','p_z']
    orb_list = list(pdos[structure[s]].keys()) # see https://pymatgen.org/pymatgen.electronic_structure.core.html
    #orb_indices = [3,1,2] # see https://pymatgen.org/pymatgen.electronic_structure.core.html
    p_dos = {}

    fig = pretty_plot(12,8)
    title = "p orbitals dos %s" % ( structure[s].species)
    fig.suptitle(title,size=30)

    color_list=full_color()
    #testvar=Orbital(orb_indices[0])
    
    for i in range(len(orb_names)):
        for j in range(len(orb_list)):
            if old_orb_names[i] in orb_list[j]:
                orb_name_in_latex = orb_names[i]
                p_dos[orb_name_in_latex] = cdos.get_site_orbital_dos(structure[s], orb_list[j])

    
    ax=fig.subplot(111)
    plotter = DosPlotter(sigma=0.1)
    plotter.add_dos_dict(p_dos)
    
    plotter.get_plot_ax(xlim=[-2, 2], ylim=[-10, 10],ax=ax, color_list=color_list)
    
    #invert_plotter = invert_axes(new_plotter)
    plt.tight_layout()
    creat_dir(path,'lob_dos_fig')
    plt.savefig('%s/lob_dos_fig/porb_dos_%s%s.png' %(path,s+1,structure[s].species), format='png', pad_inches=1, bbox_inches='tight')
    return plt

def occu_dos(cdos, structure, atom_label, emin=None, emax=None, nmin=None, nmax=None):
    dos = cdos.get_site_dos(structure[atom_label])
    densities_up = dos.get_densities(Spin(1))
    densities_down = dos.get_densities(Spin(-1))
    energies = dos.energies
    fermi = dos.efermi

    norm_min = energies[0]
    norm_max = energies[-1]
    if nmin:
        norm_min = nmin + fermi
    if nmax:
        norm_max = nmax + fermi

    occu_min = energies[0]
    occu_max = fermi
    if emin:
        occu_min = emin + fermi
    if emax:
        occu_max = emax + fermi

    norm = 0
    occu = 0
    for i in range(len(energies)):
        if energies[i] >= norm_min and energies[i] <= norm_max:
            norm = norm + densities_up[i] + densities_down[i]
        if energies[i] >= occu_min and energies[i] <= occu_max:
            occu = occu + densities_up[i] + densities_down[i]
    occu_norm = occu/norm

    #print(structure[s])
    #print(orb_names[d],'spin:' ,spin)
    #print(occu_norm, occu_min - fermi, 'eV to',occu_max - fermi, 'eV')
    print(f'{occu_norm:.4f}')
    return occu, occu_norm


def occu_d_orbs_pdos(path, s , d, spin, emin=None, emax=None, nmin=None, nmax=None):
    """
    emin, emax : the range to accumulate the occupency, with referencing to fermi level
    nmin, nmax : the range to accumulate the normalized factor, with referencing to fermi level
    """
    v = Vasprun(path + '/vasprun.xml')
    cdos = v.complete_dos
    structure = cdos.structure.sites
    
    orb_names = ["$d_{xy}$", "$d_{xz}$", "$d_{yz}$", "$d_{x^{2}-y^{2}}$", "$d_{z^{2}}$"]
    orb_indices = [4, 7, 5, 8, 6] # see https://pymatgen.org/pymatgen.electronic_structure.core.html
    d_dos = {}
    d_dos = cdos.get_site_orbital_dos(structure[s], Orbital(orb_indices[d]))
    densities = d_dos.get_densities(Spin(spin))
    energies = d_dos.energies
    fermi = d_dos.efermi

    norm_min = energies[0]
    norm_max = energies[-1]
    if nmin:
        norm_min = nmin + fermi
    if nmax:
        norm_max = nmax + fermi

    occu_min = energies[0]
    occu_max = fermi
    if emin:
        occu_min = emin + fermi
    if emax:
        occu_max = emax + fermi

    norm = 0
    occu = 0
    for i in range(len(densities)):
        if energies[i] >= norm_min and energies[i] <= norm_max:
            norm = norm + densities[i]
        if energies[i] >= occu_min and energies[i] <= occu_max:
            occu = occu + densities[i]
    occu_norm = occu/norm

    #print(structure[s])
    #print(orb_names[d],'spin:' ,spin)
    print(occu_norm, occu_min - fermi, 'eV to',occu_max - fermi, 'eV')
    return occu_norm

def find_max_density(path, s , d, spin, emin, emax):
    v = Vasprun(path + '/vasprun.xml')
    cdos = v.complete_dos
    structure = cdos.structure.sites
    
    orb_names = ["$d_{xy}$", "$d_{xz}$", "$d_{yz}$", "$d_{x^{2}-y^{2}}$", "$d_{z^{2}}$"]
    orb_indices = [4, 7, 5, 8, 6] # see https://pymatgen.org/pymatgen.electronic_structure.core.html
    d_dos = {}
    d_dos = cdos.get_site_orbital_dos(structure[s], Orbital(orb_indices[d]))
    densities = d_dos.get_densities(Spin(spin))
    energies = d_dos.energies
    fermi = d_dos.efermi

    index_list= []
    density_list= []
    energy_list= []
    max_list = []
    for i in range(len(energies)):
        if energies[i] >= emin + fermi and energies[i] <= emax + fermi:
            index_list.append(i)
            density_list.append(densities[i])
            energy_list.append(energies[i])

    for i in range(len(index_list)-2):
        density = densities[index_list[i+1]]
        if density > densities[index_list[i]] and density > densities[index_list[i+2]]:
            max_list.append((index_list[i+1], energies[index_list[i+1]] - fermi, densities[index_list[i+1]]))
    print(max_list)
    return max_list

def find_min_density(path, s , d, spin, emin, emax):
    v = Vasprun(path + '/vasprun.xml')
    cdos = v.complete_dos
    structure = cdos.structure.sites
    
    orb_names = ["$d_{xy}$", "$d_{xz}$", "$d_{yz}$", "$d_{x^{2}-y^{2}}$", "$d_{z^{2}}$"]
    orb_indices = [4, 7, 5, 8, 6] # see https://pymatgen.org/pymatgen.electronic_structure.core.html
    d_dos = {}
    d_dos = cdos.get_site_orbital_dos(structure[s], Orbital(orb_indices[d]))
    densities = d_dos.get_densities(Spin(spin))
    energies = d_dos.energies
    fermi = d_dos.efermi

    index_list= []
    density_list= []
    energy_list= []
    min_list = []
    for i in range(len(energies)):
        if energies[i] >= (emin + fermi) and energies[i] <= (emax + fermi):
            index_list.append(i)
            density_list.append(densities[i])
            energy_list.append(energies[i])

    for i in range(len(index_list)-2):
        density = densities[index_list[i+1]]
        if density < densities[index_list[i]] and density < densities[index_list[i+2]]:
            min_list.append((index_list[i+1], energies[index_list[i+1]] - fermi, densities[index_list[i+1]]))
    print(min_list)
    return min_list

