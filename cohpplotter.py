import matplotlib.pyplot as plt
from pymatgen.util.plotting import pretty_plot
from plotter_pkg.plotter_main import creat_dir
from pymatgen.electronic_structure.core import Spin
from pymatgen.electronic_structure.cohp import Cohp
from pymatgen.io.lobster.outputs import Cohpcar
from pymatgen.electronic_structure.plotter  import CohpPlotter


def get_plotnew_func(
        self,
        xlim=None,
        ylim=None,
        plot_negative=None,
        integrated=False,
        invert_axes=True,
    ):
        """
        Get a matplotlib plot showing the COHP.

        Args:
            xlim: Specifies the x-axis limits. Defaults to None for
                automatic determination.

            ylim: Specifies the y-axis limits. Defaults to None for
                automatic determination.

            plot_negative: It is common to plot -COHP(E) so that the
                sign means the same for COOPs and COHPs. Defaults to None
                for automatic determination: If are_coops is True, this
                will be set to False, else it will be set to True.

            integrated: Switch to plot ICOHPs. Defaults to False.

            invert_axes: Put the energies onto the y-axis, which is
                common in chemistry.

        Returns:
            A matplotlib object.
        """
        if self.are_coops:
            cohp_label = "COOP"
        else:
            cohp_label = "COHP"

        if plot_negative is None:
            plot_negative = not self.are_coops

        if integrated:
            cohp_label = "I" + cohp_label + " (eV)"

        if plot_negative:
            cohp_label = "-" + cohp_label

        if self.zero_at_efermi:
            energy_label = "$E - E_f$ (eV)"
        else:
            energy_label = "$E$ (eV)"

        ncolors = max(3, len(self._cohps))
        ncolors = min(9, ncolors)

        import palettable

        # pylint: disable=E1101
        colors = palettable.colorbrewer.qualitative.Set1_9.mpl_colors

        plt = pretty_plot(12, 8)

        allpts = []
        keys = self._cohps.keys()
        for i, key in enumerate(keys):
            energies = self._cohps[key]["energies"]
            if not integrated:
                populations = self._cohps[key]["COHP"]       #add a ["COHP"]
            else:
                populations = self._cohps[key]["ICOHP"]
            for spin in [Spin.up, Spin.down]:
                if spin in populations:
                    if invert_axes:
                        x = -populations[spin] if plot_negative else populations[spin]
                        y = energies
                    else:
                        x = energies
                        y = -populations[spin] if plot_negative else populations[spin]
                    allpts.extend(list(zip(x, y)))
                    if spin == Spin.up:
                        plt.plot(
                            x,
                            y,
                            color=colors[i % ncolors],
                            linestyle="-",
                            label=str(key)+ ' spin up',
                            linewidth=2,
                        )
                    else:
                        plt.plot(
                            x, y, color=colors[i+1 % ncolors], linestyle="-", linewidth=2, label=str(key) + ' spin down'
                        )
                    
        if xlim:
            plt.xlim(xlim)
        if ylim:
            plt.ylim(ylim)
        else:
            xlim = plt.xlim()
            relevanty = [p[1] for p in allpts if xlim[0] < p[0] < xlim[1]]
            plt.ylim((min(relevanty), max(relevanty)))

        xlim = plt.xlim()
        ylim = plt.ylim()
        if not invert_axes:
            plt.plot(xlim, [0, 0], "k-", linewidth=2)
            if self.zero_at_efermi:
                plt.plot([0, 0], ylim, "k--", linewidth=2)
            else:
                plt.plot(
                    [self._cohps[key]["efermi"], self._cohps[key]["efermi"]],
                    ylim,
                    color=colors[i % ncolors],
                    linestyle="--",
                    linewidth=2,
                )
        else:
            plt.plot([0, 0], ylim, "k-", linewidth=2)
            if self.zero_at_efermi:
                plt.plot(xlim, [0, 0], "k--", linewidth=2)
            else:
                plt.plot(
                    xlim,
                    [self._cohps[key]["efermi"], self._cohps[key]["efermi"]],
                    color=colors[i % ncolors],
                    linestyle="--",
                    linewidth=2,
                )

        if invert_axes:
            plt.xlabel(cohp_label)
            plt.ylabel(energy_label)
        else:
            plt.xlabel(energy_label)
            plt.ylabel(cohp_label)

        plt.legend()
        leg = plt.gca().get_legend()
        ltext = leg.get_texts()
        plt.setp(ltext, fontsize=25)
        plt.tick_params(labelsize=15)
        plt.tight_layout()
        return plt
CohpPlotter.get_plotnew = get_plotnew_func

def COHP_plot(path, n=0, o=-1, orb_label=None, xlim=[-0.5,0.5], ylim = [-10,10],savefig=False):
    '''
    f: file path
    n: No.# of pair of atoms
    o: orbital label
    '''

    cohpin= Cohpcar( filename='%sCOHPCAR.lobster' % path)
    efermi = cohpin.efermi
    energies = cohpin.energies
    
    cohp_label = list(cohpin.cohp_data.keys())
    data=cohpin.cohp_data
    orb_res_data =cohpin.orb_res_cohp

    if ( o == -1) and (not orb_label):
        cohp=data[cohp_label[n]]['COHP']
        icohp=data[cohp_label[n]]['ICOHP']
        COHP = Cohp(efermi, energies, cohp,icohp=icohp)
        Cohp_dict= {}
        Cohp_dict[cohp_label[n]]= COHP
    else:
        cohp_orb_list = orb_res_data[cohp_label[n]]
        orbital_label=list(cohp_orb_list.keys())
        if not orb_label:
            orb_label=orbital_label[o]
        cohp=cohp_orb_list[orb_label]['COHP']
        icohp=cohp_orb_list[orb_label]['ICOHP']
        COHP = Cohp(efermi, energies, cohp,icohp=icohp)
        Cohp_dict= {}
        Cohp_dict[orb_label]= COHP

    testvar = cohp[Spin.up].max()

    fig_out= False

    if (not cohp[Spin.up].max() < 1e-3) or (not cohp[Spin.down].max() <1e-3) or (not cohp[Spin.up].min() >-1e-3) or (not cohp[Spin.up].min() >-1e-3):
        fig_out= True    


    plotter= CohpPlotter()
    plotter.add_cohp_dict(Cohp_dict)
    new_plotter = plotter.get_plotnew(xlim=xlim, ylim=ylim)
    savefig_label=list(Cohp_dict.keys())[0]
    #new_plotter = plotter.get_plotnew(xlim=[-0.2, 0.2], ylim=[-15, 5],integrated= True)
    #testvar= orbital_indices[0][0]
    plt.tight_layout()
    
    if savefig and fig_out:
        creat_dir(path,'COHP_fig')
        creat_dir(path+'COHP_fig','%s' %n)
        new_plotter.savefig('%sCOHP_fig/%s/%sCOHP.png' %(path,n,savefig_label), format='png', pad_inches=1, bbox_inches='tight')
    elif not savefig and fig_out:
        plt.show()
    

    #return testvar
    return new_plotter

def ICOHP_plot(path, n=0, o=-1, orb_label=None, xlim=[-0.5,0.5], ylim = [-10,10],savefig=False):
    '''
    f: file path
    n: No.# of pair of atoms
    o: orbital label
    '''

    cohpin= Cohpcar( filename='%sCOHPCAR.lobster' % path)
    efermi = cohpin.efermi
    energies = cohpin.energies
    
    cohp_label = list(cohpin.cohp_data.keys())
    data=cohpin.cohp_data
    orb_res_data =cohpin.orb_res_cohp

    if ( o == -1) and (not orb_label):
        cohp=data[cohp_label[n]]['COHP']
        icohp=data[cohp_label[n]]['ICOHP']
        COHP = Cohp(efermi, energies, cohp,icohp=icohp)
        Cohp_dict= {}
        Cohp_dict[cohp_label[n]]= COHP
    else:
        cohp_orb_list = orb_res_data[cohp_label[n]]
        orbital_label=list(cohp_orb_list.keys())
        if not orb_label:
            orb_label=orbital_label[o]
        cohp=cohp_orb_list[orb_label]['COHP']
        icohp=cohp_orb_list[orb_label]['ICOHP']
        COHP = Cohp(efermi, energies, cohp,icohp=icohp)
        Cohp_dict= {}
        Cohp_dict[orb_label]= COHP

    testvar = cohp[Spin.up].max()

    fig_out= False

    if (not icohp[Spin.up].max() < 1e-3) or (not icohp[Spin.down].max() <1e-3) or (not icohp[Spin.up].min() >-1e-3) or (not icohp[Spin.up].min() >-1e-3):
        fig_out= True    


    plotter= CohpPlotter()
    plotter.add_cohp_dict(Cohp_dict)
    new_plotter = plotter.get_plotnew(xlim=xlim, ylim=ylim, integrated=True)
    savefig_label=list(Cohp_dict.keys())[0]
    #new_plotter = plotter.get_plotnew(xlim=[-0.2, 0.2], ylim=[-15, 5],integrated= True)
    #testvar= orbital_indices[0][0]
    plt.tight_layout()
    

    if savefig and fig_out:
        creat_dir(path,'ICOHP_fig')
        creat_dir(path+'ICOHP_fig','%s' %n)
        new_plotter.savefig('%sICOHP_fig/%s/%sICOHP.png' %(path,n,savefig_label), format='png', pad_inches=1, bbox_inches='tight')
    elif not savefig and not fig_out:
        plt.show()
    

    #return testvar
    return new_plotter
