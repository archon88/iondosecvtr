# -*- coding: utf-8 -*-
"""
################################################################################
#                                                                              #
# iondosecvtr                                                                  #
#                                                                              #
################################################################################
#                                                                              #
# LICENCE INFORMATION                                                          #
#                                                                              #
# This program provides a simple means of converting the relative ionisation   #
# at depth to the relative dose at depth for a clinical electron beam,         #
# following the Code of Practice.                                              #
#                                                                              #
# (C) 2018 Gavin Donald Kirby                                                  #
# Intial version created 2018-04-16T16:20:00Z                                  #
#                                                                              #
# This software is released under the terms of the GNU General Public License  #
# version 3 (GPLv3).                                                           #
#                                                                              #
# This program is free software: you can redistribute it and/or modify it      #
# under the terms of the GNU General Public License as published by the Free   #
# Software Foundation, either version 3 of the License, or (at your option)    #
# any later version.                                                           #
#                                                                              #
# This program is distributed in the hope that it will be useful, but WITHOUT  #
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or        #
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for     #
# more details.                                                                #
#                                                                              #
# For a copy of the GNU General Public License, see                            #
# <http://www.gnu.org/licenses/>.                                              #
#                                                                              #
################################################################################

"""

import math as m
import os

import pandas as pd
import matplotlib.pyplot as plt

__version__ = '2019-07-07T23:00:00Z'


def calculate_spr(zw, r50d):
    """
    This function returns the water-to-air mass stopping power ratio
    for a clinical electron beam of quality index r50d (cm) in water
    of depth zw (cm).
    """
    # Exception handling -- if provided data fall outside
    # the range of validity of the empirical formula.
    if not 1 < r50d < 20:
        raise ValueError('R_(50,D) must be between 1 cm and 20 cm!')    
    if not 0.02*r50d < zw < 1.2*r50d:
        raise ValueError(f'z_w must be between 0.02*R_(50,D) and 1.2*R_(50,D)!\
        Current z_w is {zw}, which falls outside the range ({0.02*r50d}, {1.2*r50d}).')
    
    x = m.log(r50d)
    y = zw/r50d
    # Use empirical formula from Code of Practice.
    return (1.075 - 0.5087*x + 0.0887*pow(x, 2) - 0.084*y) / (1 - 0.4281*x + 0.0646*pow(x, 2) + 0.00309*pow(x, 3) - 0.125*y)    


def convert_ionisation_dose(file, r50d):
    """
    Convert ionisation measurements at depth to dose at depth, using the
    empirical formula in calculate_spr. A readable CSV containing the
    ionisation at depth data is required, as is an r50d in centimetres.
    """
    try:
        df = pd.read_csv(file)
    except Exception:
        raise ValueError(f'File {file} is not a readable csv file!')
    
    # Normalise ionisation to a maximum of 100%.
    df.ion = 100*(df.ion / df.ion.max())
    # Remove inappropriate depths then pass depths and
    # relative ionisation at depth values to calculate_spr.
    df = df[0.02*r50d < 0.1*df['depth']]
    df = df[0.1*df['depth'] < 1.2*r50d]
    df['spr']=[calculate_spr(z, r50d) for z in 0.1*df.depth]
    df['dose']=df.ion*[calculate_spr(z, r50d) for z in 0.1*df.depth]
    
    # Normalise dose to 100% at its maximum.
    df.dose = 100*(df.dose/df.dose.max())
    # Round each of the dataframe's columns appropriately.
    df.ion = round(df.ion, 1)
    df.spr = round(df.spr, 3)
    df.dose = round(df.dose, 1)
    # Now write the calculated dose data to file.
    outfile = os.path.splitext(file)[0] + '_withdose.csv'    
    df.to_csv(outfile, index=False)


def make_plots(file, energy):
    """
    Accept a processed CSV file with added dose at depth information,
    together with a nominal energy in MeV, and create plots of
    ionisation, SPR, and dose at depth.
    """
    try:
        df = pd.read_csv(file)
    except Exception:
        raise ValueError(f'File {file} is not a readable csv file!')
        
    depth_dose_fig = df.plot(x='depth', y='dose', legend=False)
    depth_dose_fig.set_title(f'{energy} MeV electron beam percentage depth dose')
    depth_dose_fig.set_xlabel("Water depth (mm)")
    depth_dose_fig.set_ylabel("Dose at depth (%)")    
    plotname = os.path.splitext(file)[0] + '_PDD.png'    
    plt.savefig(plotname, dpi=900)
    
    depth_ion_fig = df.plot(x='depth', y='ion', legend=False)
    depth_ion_fig.set_title(f'{energy} MeV electron beam percentage depth ionisation')
    depth_ion_fig.set_xlabel("Water depth (mm)")
    depth_ion_fig.set_ylabel("Ionisation at depth (%)")    
    plotname = os.path.splitext(file)[0] + '_PDI.png'    
    plt.savefig(plotname, dpi=900)
    
    depth_spr_fig = df.plot(x='depth', y='spr', legend=False)
    depth_spr_fig.set_title(f'{energy} MeV electron beam mass stopping power ratio (SPR) at depth')
    depth_spr_fig.set_xlabel("Water depth (mm)")
    depth_spr_fig.set_ylabel("SPR at depth")    
    plotname = os.path.splitext(file)[0] + '_SPR.png'    
    plt.savefig(plotname, dpi=900)
 

if __name__ == "__main__":
    print(f'This is iondosecvtr version {__version__}.')
    exit()
