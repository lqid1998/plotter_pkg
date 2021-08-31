from pymatgen.electronic_structure.core import Spin
import matplotlib.pyplot as plt
import numpy as np

def get_plot_ax_func(self, xlim=None, ylim=None,ax=None, color_list=None):
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
                x.extend(energies)
                y.extend(densities)
        allpts.extend(list(zip(x, y)))      #some plt. were replaced by ax.
        if color_list_flag:
            if self.stack:
                ax.fill(x, y, color=color_list[key], label=str(key))        #use ax.fill instead of plt.fill
            else:
                ax.plot(x, y, color=color_list[key], label=str(key), linewidth=3)      #use ax.plot to plot the dos into the target axes, instead of plt.plot() directly 
        else:
            if self.stack:
                ax.fill(x, y, color=colors[i % ncolors], label=str(key))        #use ax.fill instead of plt.fill
            else:
                ax.plot(x, y, color=colors[i % ncolors], label=str(key), linewidth=3)      #use ax.plot to plot the dos into the target axes, instead of plt.plot() directly 
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

    #ax.set_xlabel("Energies (eV)")      #use ax.set_xlabel ax.set_ylabel instead of plt.xlabel plt.ylabel
    #ax.set_ylabel("Density of states")

    ax.axhline(y=0, color="k", linestyle="--", linewidth=2)     #use ax.axhline instead of plt.axhline
    ax.legend(fontsize=30)         #use ax.legend() insteadof plt.legend

    #These lines can be replaced by ax.legend(fontsize=30)
    #leg = plt.gca().get_legend()
    #ltext = leg.get_texts()  # all the text.Text instance in the legend
    #plt.setp(ltext, fontsize=30)

    #This line is moved to d_orbs_subdos()
    #plt.tight_layout()

    #Modified from invert_axes(plotter) to exchange the x-axis and y-axis of the plot
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

    return ax

def set_color(orb_list):
    import palettable
    orb_index={"$d_{xy}$": 0, "$d_{xz}$": 1, "$d_{yz}$": 2, "$d_{x^{2}-y^{2}}$": 3, "$d_{z^{2}}$": 4\
        ,"$s$": 5,"$p_{x}$": 6,"$p_{y}$": 7,"$p_{z}$": 8,}
    colors = palettable.colorbrewer.qualitative.Set1_9.mpl_colors
    color_list={}
    # pylint: disable=E1101
    for i in range(len(orb_list)):
        orb_name=orb_list[i]
        color_list[orb_name]=colors[orb_index[orb_name]]
    return color_list