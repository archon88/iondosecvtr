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

import numpy as np
import math as m
import pandas as pd
import os
import matplotlib.pyplot as plt

__version__ = '2018-04-16T16:20:00Z'

def sprcalc(zw,r50d):
    
    if not 1 < r50d < 20:
        raise ValueError('R_(50,D) must be between 1 cm and 20 cm!')
    
    if not 0.02*r50d < zw < 1.2*r50d:
        raise ValueError('z_w must be between 0.02*R_(50,D) and 1.2*R_(50,D)! Current zw is {zw}, which falls outside the range ({lower},{upper}).'.format(zw=zw,lower=0.02*r50d,upper=1.2*r50d))
    
    x = m.log(r50d)
    y = zw/r50d
    
    spr = (1.075 -0.5087*x + 0.0887*np.power(x,2) - 0.084*y)/(1 - 0.4281*x + 0.0646*np.power(x,2) + 0.00309*np.power(x,3)  - 0.125*y)
    
    return spr

def cvtiondose(file,r50d):
    
    try:
        df = pd.read_csv(file)
    except Exception as e:
        raise ValueError('File {file} is not a readable csv file!'.format(file=file))
    
    #Ensure ionisation is normalised to a max of 100
    df.ion = 100*(df.ion/df.ion.max())
    
    #Remove inappropriate depths then pass depths and relative ionisation at depth values to sprcalc
    df = df[0.02*r50d < 0.1*df['depth']]
    df = df[0.1*df['depth'] < 1.2*r50d]
    df['spr']=[sprcalc(z,r50d) for z in 0.1*df.depth]
    df['dose']=df.ion*[sprcalc(z,r50d) for z in 0.1*df.depth]
    
    #Normalise dose to 100 at max
    df.dose = 100*(df.dose/df.dose.max())
    
    df.ion = round(df.ion,1)
    df.spr = round(df.spr,3)
    df.dose = round(df.dose,1)
    
    outfile=os.path.splitext(file)[0]+'_withdose.csv'
    
    df.to_csv(outfile, index=False)
    
    return None

def makeplots(file, energy):
    
    try:
        df = pd.read_csv(file)
    except Exception as e:
        raise ValueError('File {file} is not a readable csv file!'.format(file=file))
        
    fig1=df.plot(x='depth',y='dose',legend=False)
    fig1.set_title('{energy} MeV electron beam percentage depth dose'.format(energy=energy))
    fig1.set_xlabel("Water depth (mm)")
    fig1.set_ylabel("Dose at depth (%)")
    
    plotname=os.path.splitext(file)[0]+'_PDD.png'
    
    plt.savefig(plotname,dpi=900)
    
    fig2=df.plot(x='depth',y='ion',legend=False)
    fig2.set_title('{energy} MeV electron beam percentage depth ionisation'.format(energy=energy))
    fig2.set_xlabel("Water depth (mm)")
    fig2.set_ylabel("Ionisation at depth (%)")
    
    plotname=os.path.splitext(file)[0]+'_PDI.png'
    
    plt.savefig(plotname,dpi=900)
    
    fig3=df.plot(x='depth',y='spr',legend=False)
    fig3.set_title('{energy} MeV electron beam mass stopping power ratio (SPR) at depth'.format(energy=energy))
    fig3.set_xlabel("Water depth (mm)")
    fig3.set_ylabel("SPR at depth")
    
    plotname=os.path.splitext(file)[0]+'_SPR.png'
    
    plt.savefig(plotname,dpi=900)
    
 
if __name__ == "__main__":
    print('This is iondosecvtr version {version}'.format(version=__version__))
    exit()
