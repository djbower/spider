#!/usr/bin/env python

import logging
import spider_utils as su
import matplotlib.transforms as transforms
import matplotlib.pyplot as plt
import numpy as np
import os
import matplotlib.ticker as ticker
import json

logger = su.get_my_logger(__name__)

#====================================================================
def plot_atmosphere():

    # reads data dumped by spider_time_data

    logger.info( 'building atmosphere' )

    width = 4.7747 / 1.0
    height = 4.7747 / 2.0

    with open('data.json', 'r') as fp:
        data_d = json.load(fp)

    fig_o = su.FigureData( 1, 2, width, height, 'Case1_Case9_atmosphere', units='Myr' )
    fig_o.fig.subplots_adjust(wspace=0.3,hspace=0.3)

    ax0 = fig_o.ax[0]
    ax1 = fig_o.ax[1]

    title = r'\textbf{(a) Case 1, }' + r'$\bf{r_{f,CO_2}=0.1}$'
    label_align_l = ('right','left','left','center','center')
    ylabeloff = 0.5
    plot_figure( title, fig_o, ax0, 1, label_align_l, ylabeloff, data_d )

    title = r'\textbf{(b) Case 9, }' + r'$\bf{r_{f,CO_2}=0.9}$'
    label_align_l = ('right','center','left','center','center')
    ylabeloff = 0.54
    plot_figure( title, fig_o, ax1, 9, label_align_l, ylabeloff, data_d )

    fig_o.savefig(1)

#====================================================================
def plot_figure( title, fig_o, ax, case, label_align_l, ylabeloff, data_d ):

    strcase = 'case'+str(case)
    Case_time_myrs = np.array(data_d[strcase]['time_l']) * 1.0E-6
    Case_CO2_atmos_a = np.array(data_d[strcase]['CO2_atmosphere_bar'])
    Case_H2O_atmos_a = np.array(data_d[strcase]['H2O_atmosphere_bar'])
    Case_phi_global = np.array(data_d[strcase]['phi_l'])

    # to plot melt fraction contours on figure (a)
    # compute time at which desired melt fraction is reached
    phi_cont_l = [0.75, 0.5,0.25,0.1,0.01]
    phi_time_l = [] # contains the times for each contour
    for cont in phi_cont_l:
        time_temp_l = su.find_xx_for_yy( Case_time_myrs, Case_phi_global, cont )
        index = su.get_first_non_zero_index( time_temp_l )
        if index is None:
            out = 0.0
        else:
            out = Case_time_myrs[index]
        phi_time_l.append( out )

    xticks = [1.0E-2, 1.0E-1, 1.0E0, 1.0E1, 1.0E2,1.0E3] #[1E-6,1E-4,1E-2,1E0,1E2,1E4,1E6]#,1]
    xlabel = 'Time (Myr)'
    xlim = (1.0E-2, 4550)

    red = (0.5,0.1,0.1)
    blue = (0.1,0.1,0.5)
    black = 'black'

    ylabel = '$p$\n(bar)'
    trans = transforms.blended_transform_factory( ax.transData, ax.transAxes)
    h1, = ax.semilogx( Case_time_myrs, Case_CO2_atmos_a, color=red, linestyle='-', label=r'CO$_2$')
    h2, = ax.semilogx( Case_time_myrs, Case_H2O_atmos_a, color=blue, linestyle='-', label=r'H$_2$O')
    #axb = ax.twinx()
    #h3, = axb.semilogx( Case_time_myrs, Case_phi_global, color=black, linestyle='--', label=r'Melt, $\phi_g$')
    for cc, cont in enumerate(phi_cont_l):
        ax.axvline( phi_time_l[cc], ymin=0.02, ymax=0.95, color='0.25', linestyle=':' )
        label = cont #int(cont*100) # as percent
        label = str(round(label,2))
        if cc == 0:
            label= r'$\phi_g=$ ' + str(label)
        ha = label_align_l[cc]
        #ax.text( phi_time_l[cc], 0.5, '{0:.2f}'.format(label), va='top', ha=ha, rotation=90, bbox=dict(facecolor='white', edgecolor='none', pad=2), transform=trans )
        ax.text( phi_time_l[cc], 0.6, label, va='top', ha=ha, rotation=90, bbox=dict(facecolor='white', edgecolor='none', pad=2), transform=trans )
    handle_l = [h1,h2]#,h3]
    fig_o.set_myaxes( ax, title=title, ylabel=ylabel, xlabel=xlabel, xticks=xticks )
    if case==1 or case==9:
        fig_o.set_mylegend( ax, handle_l, loc='upper left', ncol=1 )
        #fig_o.set_mylegend( ax, handle_l, loc=(0.02,0.1), ncol=1 )
    ax.xaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=20) )
    ax.xaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=(0.2,0.4,0.6,0.8), numticks=20))
    ax.xaxis.set_minor_formatter(ticker.NullFormatter())
    ax.set_xlim( *xlim )
    ax.yaxis.set_label_coords(-0.1,ylabeloff)
    #axb.set_ylabel( r'$\phi_g$', rotation=0 )
    #axb.yaxis.set_label_coords(1.1,0.525)

#====================================================================
def main():

    plot_atmosphere()
    plt.show()

#====================================================================

if __name__ == "__main__":

    main()
