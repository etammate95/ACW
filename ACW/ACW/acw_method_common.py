#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 11:10:59 2022

@author: raxis
"""

import os
import math
import numpy as np
import platform
import subprocess
import shutil
from itertools import combinations

import pymol
from pymol import cmd
from pymol import stored
from pymol.wizard import Wizard

from .basic import *
from .pdbdata import *
from .pdbcoord import *

#notes for calc methods

#areaimol_input for acw, acwperf, f1,f2,f3,f4,f5
#areaimol_calc_1 for acw, acwperf, f1,f2
#areaimol_calc_2 for f3,f4,f5
#areaimol_out_1 for acw, acwperf
#areaimol_out_2 for f1, f2
#areaimol_out_3 for f3, f4, f5
#sc_input for all
#sc_calc for all
#sc_out for all
#sc_calc_out subprocess style sc
#areaimol_calc_3 subprocess style area calc_1

#A-D
#B-E
#C-F

#A-G
#B-H
#C-I
#D-J
#E-K
#F-L


#A-K
#B-L
#C-M
#D-N
#E-O
#F-P
#G-Q
#H-R
#I-S
#J-T




def get_data(sele, var):
    """gets data of named selection from pymol as list"""

    stored.data = []
    cmd.iterate_state(-1, sele, 'stored.data.append(%s)'%var)

    return stored.data


def draw_axis():
    """Draws the main coordinate system"""

    for n, a in enumerate(['X_axis', 'Y_axis', 'Z_axis']):
        coord = [0, 0, 0]
        coord[n] = 30
        cmd.pseudoatom('axis', pos=coord, name=a[0])
        cmd.label('axis and name %s'%a[0], 'name')
        coord[n] = -30
        cmd.pseudoatom('axis', pos=coord, name=a[0]*2)
        cmd.distance(a, 'axis and name '+a[0], 'axis and name '+a[0]*2)
        cmd.color(['red', 'green', 'blue'][n], a)
        cmd.hide('labels', a)
    cmd.hide('wire', 'axis')



def straigthen_pca_check(name, name2):
    """Checks if the current selection is aligned to Y."""

    """
    #sheet 0 middle CA
    sheet = pdb_reader('Sheet_0.pdb').splice_same()

    seq_len = len(sheet.splice()[0].residues())

    atoms = []
    for atom in sheet.select('resi', str(int(seq_len/2.0))).select('name', ['CA']).lista:
        atoms.append(atom.xyz())
    """

    #both sheet centrums
    sheet_0 = pdb_reader('Sheet_0.pdb').splice_same().splice()
    sheet_1 = pdb_reader('Sheet_1.pdb').splice_same().splice()
    atoms = []
    for n in range(len(sheet_0)):
        atoms.append(get_midpoint_coord(sheet_0[n].center(),sheet_1[n].center()))

    #calc for rotate axis to Y
    data , eigenvector_subset = PCA(atoms, 1)

    vector = [x[0] for x in eigenvector_subset]

    v = np.array(vector)
    v = v / np.sqrt(np.sum(v ** 2)) #norm
    angle = - math.degrees(math.acos(np.dot((0, 1, 0), v))) #angle from dot
    v_cross = list(np.cross((0, 1, 0), v)) #axis from cross

    #limit to skip
    if angle > -0.5 or angle < -179.5:
        return None, None

    cmd.rotate(v_cross, angle=angle, selection='all', camera=0, object=None, origin=[0, 0, 0])
    cmd.save('Sheet_0.pdb', name)
    cmd.save('Sheet_1.pdb', name2)
    cmd.rotate(v_cross, angle=-angle, selection='all', camera=0, object=None, origin=[0, 0, 0])
    return v_cross, angle


def config_reader(self):
    """Reads the config file."""

    try:
        with open("config.txt", 'r+') as f:
            for line in f.read().splitlines():
                if line:
                    if line[0] != '#' and len(line.split()) == 3:
                        if line.split()[1] == '=':
                            line = line.split()
                            setattr(self, line[0], line[2])

    except FileNotFoundError:
        pass


def sc_input(count, name):
    # combine the sheet PDB files into one

    with open("%s_sheets.pdb"%name, 'w+') as f:

        for n in range(count+1):
            #each sheet is a different chain
            with open("Sheet_%s.pdb"%n, 'r+') as g:
                for line in g.read().splitlines():
                    #this deletes ANISOU
                    if line[0:6] == 'ATOM  ' or line[0:6] == 'HETATM':
                        f.write(line[0:21] + ch_n(n) + line[22:])
                        f.write('\n')
                f.write('TER\n')

        f.write('END')


def sc_calc_out(count, name, path):

    if platform.system() == 'Windows':
        p = subprocess.Popen(['cmd.exe'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True, shell=True)
        p.stdin.write('%s/ccp4.setup.bat\n'%path)
        result = p.communicate('%s/bin/sc.exe XYZIN %s_sheets.pdb > sc_%s_sheet%s.txt\nMOLECULE 1\nCHAIN A\nMOLECULE 2\nCHAIN %s\nEND\n'%(path, name, name, ch_n(count), ch_n(count)))
    elif platform.system() == 'Linux':
        p = subprocess.Popen(['/bin/sh'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
        result = p.communicate('%s/bin/sc XYZIN %s_sheets.pdb > sc_%s_sheet%s.txt\nMOLECULE 1\nCHAIN A\nMOLECULE 2\nCHAIN %s\nEND\n'%(path, name, name, ch_n(count), ch_n(count)))

    sc = 'NaN'
    area1 = 'NaN'
    area2 = 'NaN'

    with open("sc_%s_sheet%s.txt"%(name, ch_n(count)), 'r+') as f:

        for line in f.read().splitlines():
            if line[0:41] == ' Shape complementarity statistic Sc =    ':
                sc = line[41:46]
            elif line[0:50] == '  Total area left after trim for this molecule is ':
                if area1 == 'NaN':
                    area1 = line[50:58].strip()
                else:
                    area2 = line[50:58].strip()

    return (sc, area1, area2)



