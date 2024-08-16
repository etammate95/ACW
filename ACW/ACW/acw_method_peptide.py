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
from .acw_method_common import *

def selector(start, end, d):
    """creates text for sheet picking, unused"""

    txt = '(bychain resi %s-%s and (name O or name N) and polymer '%(start, end)
    txt += 'within %s of sele '%(d)
    txt += 'and resi %s-%s and (name O or name N)) and polymer'%(start, end)

    return txt


def straigthen_pca(name):
    """Rotates main axis to Y and orient iface with PCA."""

    cmd.select(name, name + ' and polymer', 0, 1, 0, 1)
    cmd.save('%s.pdb'%name, name)

    #find vector with parallel vector method
    sheet = pdb_reader('%s.pdb'%name)
    resi = sheet.lista[0].resi #first residue
    A = sheet.select('resi', resi).select('name', ['C']).lista[0].xyz()
    vectors = []
    for atom in sheet.select('resi', resi).select('name', ['C']).lista[1:]:
        vec=[A[0]-atom.x, A[1]-atom.y, A[2]-atom.z]
        vectors.append(vec)
    for comb in combinations(vectors, 2):
        if get_lc(np.cross(comb[0], comb[1]), [0, 0, 0]) < 0.001:
            vector = comb[0][:]
            vector_p = comb[0][:]
            break
    else:
        #find vector with pca
        sheet = pdb_reader('%s.pdb'%name).splice_same()
        seq_len = len(sheet.splice()[0].residues())
        atoms = []
        for atom in sheet.select('resi', str(int(seq_len/2.0))).select('name', ['CA']).lista:
            atoms.append(atom.xyz())
        data , eigenvector_subset = PCA(atoms, 1)
        vector = [x[0] for x in eigenvector_subset]

    #rotate axis to Y: middle amino acid CA
    v = np.array(vector)
    v = v / np.sqrt(np.sum(v ** 2)) #norm
    angle = - math.degrees(math.acos(np.dot((0, 1, 0), v))) #angle from dot
    v_cross = list(np.cross((0, 1, 0), v)) #axis from cross

    cmd.rotate(v_cross, angle=angle, selection='all', camera=0, object=None, origin=[0, 0, 0])

    #orient backbones to Z
    cmd.select('sele_x', 'sele_x and polymer', 0, 1, 0, 1)
    cmd.save('Sheet_X.pdb', 'sele_x')
    sheet = pdb_reader('Sheet_X.pdb').splice_mid()

    sheet_flat = flatten(sheet, 'y')

    atoms = []
    for atom in sheet_flat.backbone().lista:
        atoms.append(atom.xyz())

    data , eigenvector_subset = PCA(atoms, 1)
    vector = [x[0] for x in eigenvector_subset]

    vector[1] = 0
    v = np.array(vector)
    v = v / np.sqrt(np.sum(v ** 2)) #norm
    angle = - math.degrees(math.acos(np.dot((0, 0, 1), v))) #angle from dot
    v_cross = list(np.cross((0, 0, 1), v)) #axis from cross

    cmd.rotate(v_cross, angle=angle, selection='all', camera=0, object=None, origin=[0, 0, 0])

    """
    #display
    coord = [atoms[-1][n]-10*eigenvector_subset[n] for n in [0, 1, 2]]
    cmd.pseudoatom('aaa', pos=coord, name='aa')
    coord = [atoms[-1][n]+10*eigenvector_subset[n] for n in [0, 1, 2]]
    cmd.pseudoatom('bbb', pos=coord, name='bb')
    cmd.distance('line', 'aaa and name aa', 'bbb and name bb')
    """

def perf_calc_fast(count, S=None, B=None):

    #import acwperf
    #cleverly find points of possible water positions
    #adds a circle of points to the surface points, then with union finds the output

    #checks the orientation
    if S is not None:
        v_cross, angle = straigthen_pca_check('Bsheet_%s'%S, 'Bsheet_%s'%B)
    else:
        angle = None

    name='tmp'

    modelA = pdb_reader('Sheet_0.pdb')
    modelB = pdb_reader('Sheet_%s.pdb'% str(count))

    axis = find_axis(modelA)
    project = {'x':0,'y':1,'z':2}[axis]
    avg_height = modelA.average(axis)

    if len(modelA.splice()) == 3:
        model1 = modelA.splice()[1]
        model2 = modelB.splice()[1]
        height = get_length_coord(modelA.splice()[1].center(), modelA.splice()[2].center())
    elif len(modelA.splice()) == 6:
        model1 = modelA.splice()[2]
        model1.addmodel(modelA.splice()[3])
        model2 = modelB.splice()[2]
        model2.addmodel(modelB.splice()[3])
        height = get_length_coord(modelA.splice()[2].center(), modelA.splice()[4].center())/2
    else:
        return 0, 0

    flat1 = flatten(model1, axis)
    flat2 = flatten(model2, axis)

    radius = 3.2 # 1.4+1.8=3.2
    limit = 0.01 #rounding limit was set to 0.05
    print('fast')
    grid_1d = np.arange(-1*radius - limit, radius + limit*2, limit)
    X, Y = np.meshgrid(grid_1d, grid_1d)
    grid_points = np.column_stack((X.ravel(), Y.ravel()))
    grid_distance = np.sum(np.square(grid_points), axis=1)
    mask_d = grid_distance > (radius-limit)**2
    mask_D = grid_distance < (radius+limit)**2
    mask = np.all([mask_d, mask_D], axis=0)

    circle = grid_points[mask]

    mask_radius = (radius - limit)**2
    #for each surface generate the outer surface points
    surface1 = np.empty((0,2))
    for atom1 in flat1.lista:
        x1 = round(atom1.x * 100) / float(100)
        y1 = round(atom1.y * 100) / float(100)
        z1 = round(atom1.z * 100) / float(100)
        c_atom1 = np.delete(np.array([x1, y1, z1]), project)
        closest = []
        for atom2 in flat1.lista: #find close atoms
            d = ((x1-atom2.x)**2 + (y1-atom2.y)**2 + (z1-atom2.z)**2 )** 0.5
            if d < 2*radius+limit:
                closest.append(atom2.xyz(axis))
                closest[-1].append(d)

        closest.sort(key=lambda x:x[2])
        closest = np.array([x[0:2] for x in closest])
        shifted_circle = circle + c_atom1 #put the cicrle around the atom

        for point in closest:
            #print(shifted_circle.shape)
            #distance = np.linalg.norm(shifted_circle - point, axis=1)
            #distance = np.sqrt(np.sum((shifted_circle - point) ** 2, axis=0))
            a_min_b = shifted_circle - point
            distance = np.einsum("ij,ij->i", a_min_b, a_min_b)
            mask = distance > mask_radius
            shifted_circle = shifted_circle[mask]
            if shifted_circle.size == 0:
                break
        else:
            surface1 = np.append(surface1, shifted_circle, axis=0)

    surface2 = np.empty((0,2))
    for atom1 in flat2.lista:
        x1 = round(atom1.x * 100) / float(100)
        y1 = round(atom1.y * 100) / float(100)
        z1 = round(atom1.z * 100) / float(100)
        c_atom1 = np.delete(np.array([x1, y1, z1]), project)
        closest = []
        for atom2 in flat2.lista: #find close atoms
            d = ((x1-atom2.x)**2 + (y1-atom2.y)**2 + (z1-atom2.z)**2 )** 0.5
            if d < 2*radius+limit:
                closest.append(atom2.xyz(axis))
                closest[-1].append(d)

        closest.sort(key=lambda x:x[2])
        closest = np.array([x[0:2] for x in closest])
        shifted_circle = circle + c_atom1 #put the cicrle around the atom

        for point in closest:
            a_min_b = shifted_circle - point
            distance = np.einsum("ij,ij->i", a_min_b, a_min_b)
            mask = distance > mask_radius
            shifted_circle = shifted_circle[mask]
            if shifted_circle.size == 0:
                break
        else:
            surface2 = np.append(surface2, shifted_circle, axis=0)

    tangents = []
    surface1 = set([tuple(x) for x in np.round(surface1,2)])
    surface2 = set([tuple(x) for x in np.round(surface2,2)])
    for point in surface1: #find common points
        if point in surface2:
            tangents.append(point)
    tangents = list(set(tangents))

    """
    #with savefile
    dist = []
    with open("perf_coord_%s_fast.txt" %str(count), 'w+') as f:
        for n, coord1 in enumerate(tangents):
            atom1 = sAtom(name, coord1[0], coord1[1], coord1[2])
            f.write(str(atom1.x)+' '+str(atom1.y)+' '+str(atom1.z)+'\n')
            for coord2 in tangents[n+1:]:
                atom2 = sAtom(name, coord2[0], coord2[1], coord2[2])
                dist.append(get_length(atom1, atom2, 'lengths'))

        if len(tangents) == 0:
            length = 0
        else:
            length = max(dist)-2.8  #testing 2.8
        f.write(str(round(length, 3)))
    """

    #with output
    dist = []
    dist_c = []
    for n, coord1 in enumerate(tangents):
        coord1 = list(coord1)
        coord1.insert(project, 0)
        atom1 = sAtom(name, coord1[0], coord1[1], coord1[2])
        #(str(atom1.x)+' '+str(atom1.y)+' '+str(atom1.z)+'\n')
        for coord2 in tangents[n+1:]:
            coord2 = list(coord2)
            coord2.insert(project, 0)
            atom2 = sAtom(name, coord2[0], coord2[1], coord2[2])
            dist.append(get_length(atom1, atom2))
            dist_c.append([atom1, atom2])

    if len(tangents) == 0:
        length = 0
        length_coord = ""
    else:
        #length = max(dist)-2.8  #two times water 1.4*2

        atom1, atom2 = dist_c[dist.index(max(dist))]
        #cmd.pseudoatom('water2', pos=atom2.xyz())
        for atom in (atom1, atom2):
            flat1_c = coord_reader(vector_segment(atom, get_nearest(atom, flat1, True)[2], 1.4))
            flat2_c = coord_reader(vector_segment(atom, get_nearest(atom, flat2, True)[2], 1.4))
            contact_p = coord_reader(vector_segment(atom, get_midpoint(flat1_c, flat2_c), 1.4))
            for xyz in ('x', 'y', 'z'):
                setattr(atom, xyz, getattr(contact_p, xyz))

        #cmd.pseudoatom('contact_p', pos=contact_p.xyz())
        #cmd.pseudoatom('flat1_c', pos=flat1_c.xyz())
        #cmd.pseudoatom('flat2_c', pos=flat2_c.xyz())
        length = get_length(atom1, atom2)

        #save coordinates
        setattr(atom1, axis, avg_height)
        setattr(atom2, axis, avg_height)
        cmd.pseudoatom('tmp1', pos=atom1.xyz())
        cmd.pseudoatom('tmp2', pos=atom2.xyz())
        cmd.set('suspend_updates', 0, quiet=1)
        if angle is not None:
            cmd.rotate(v_cross, angle=-angle, selection='tmp1 or tmp2', camera=0, object=None, origin=[0, 0, 0])

        atom1 = coord_reader(get_data('tmp1','[x, y, z]')[0], name)
        atom2 = coord_reader(get_data('tmp2','[x, y, z]')[0], name)
        length_coord = '%s_%s_%s_%s_%s_%s'%(round(atom1.x, 3), round(atom1.y, 3), round(atom1.z, 3),
                                            round(atom2.x, 3), round(atom2.y, 3), round(atom2.z, 3))

        #drawing a line
        cmd.distance('lengths','tmp1', 'tmp2')
        cmd.delete('tmp1')
        cmd.delete('tmp2')
        cmd.hide('label', 'lengths')





    return length, height, length_coord




def perf_calc_fast_old(count, S=None, B=None):
    """Calculates SDi with a faster method."""


    #cleverly find points of possible water positions
    #adds a circle of points to the surface points, then with union finds the output

    #checks the orientation (saved in common methods)
    if S is not None:
        v_cross, angle = straigthen_pca_check('Bsheet_%s'%S, 'Bsheet_%s'%B)
    else:
        angle = None

    name='tmp'

    modelA = pdb_reader('Sheet_0.pdb')
    modelB = pdb_reader('Sheet_%s.pdb'% str(count))

    axis = find_axis(modelA)
    avg_height = modelA.average(axis)

    if len(modelA.splice()) == 3:
        model1 = modelA.splice()[1]
        model2 = modelB.splice()[1]
        height = get_length_coord(modelA.splice()[1].center(), modelA.splice()[2].center())
    elif len(modelA.splice()) == 6:
        model1 = modelA.splice()[2]
        model1.addmodel(modelA.splice()[3])
        model2 = modelB.splice()[2]
        model2.addmodel(modelB.splice()[3])
        height = get_length_coord(modelA.splice()[2].center(), modelA.splice()[4].center())/2
    else:
        return 0, 0

    flat1 = flatten(model1, axis)
    flat2 = flatten(model2, axis)

    radius = 3.2 # 1.4+1.8=3.2
    limit = 0.01 #rounding limit was set to 0.05

    """
    sphere = [] #generate a sphere around 0, 0
    for x in ranger(-1*radius -limit*2, radius+limit*2, limit):
        for y in ranger(-1*radius -limit*2, radius+limit*2, limit):
            for z in ranger(-1*radius -limit*2, radius+limit*2, limit):
                d = ((x)**2 + (y)**2 + (z)**2 )** 0.5
                if d<radius+limit and d>radius-limit:
                    sphere.append([x, y, z])
    """

    circle = [] #generate a cicrle around 0, 0
    for x in ranger(-1*radius -limit*2, radius+limit*2, limit):
        for y in ranger(-1*radius -limit*2, radius+limit*2, limit):
            d = ((x)**2 + (y)**2)** 0.5
            if d<radius+limit and d>radius-limit:
                project = {'x':0, 'y':1, 'z':2}[axis]
                point = [x, y]
                point.insert(project, 0)
                circle.append(tuple(point))

    calcs = 0

    #for each surface generate the outer surface points
    surface1 = []
    for atom1 in flat1.lista:
        x1 = round(atom1.x * 100) / float(100)
        y1 = round(atom1.y * 100) / float(100)
        z1 = round(atom1.z * 100) / float(100)
        closest = []
        for atom2 in flat1.lista: #find close atoms
            d = ((x1-atom2.x)**2 + (y1-atom2.y)**2 + (z1-atom2.z)**2 )** 0.5
            calcs += 1
            if d < 2*radius+limit:
                closest.append(atom2)
        for point in circle: #put the cicrle around the atom
            x = x1+point[0]
            y = y1+point[1]
            z = z1+point[2]
            for atom2 in closest: #chek if its too close to other atoms
                d = ((x-atom2.x)**2 + (y-atom2.y)**2 + (z-atom2.z)**2 )** 0.5
                calcs += 1
                if d < radius-limit:
                    break
            else:
                surface1.append((round(x * 100) / 100, round(y * 100) / 100, round(z * 100) / 100))

    surface2 = []
    for atom1 in flat2.lista:
        x1 = round(atom1.x * 100) / float(100)
        y1 = round(atom1.y * 100) / float(100)
        z1 = round(atom1.z * 100) / float(100)
        closest = []
        for atom2 in flat2.lista:
            d = ((x1-atom2.x)**2 + (y1-atom2.y)**2 + (z1-atom2.z)**2 )** 0.5
            calcs += 1
            if d < 2*radius+limit:
                closest.append(atom2)
        for point in circle:
            x = x1+point[0]
            y = y1+point[1]
            z = z1+point[2]
            for atom2 in closest:
                d = ((x-atom2.x)**2 + (y-atom2.y)**2 + (z-atom2.z)**2 )** 0.5
                calcs += 1
                if d < radius-limit:
                    break
            else:
                surface2.append((round(x * 100) / 100, round(y * 100) / 100, round(z * 100) / 100))

    #print('calcs:', calcs)

    surface1 = set(surface1)
    surface2 = set(surface2)
    tangents = []
    for point in surface1: #find common points
        if point in surface2:
            tangents.append(point)

    tangents = list(set(tangents))

    """
    #with savefile
    dist = []
    with open("perf_coord_%s_fast.txt" %str(count), 'w+') as f:
        for n, coord1 in enumerate(tangents):
            atom1 = sAtom(name, coord1[0], coord1[1], coord1[2])
            f.write(str(atom1.x)+' '+str(atom1.y)+' '+str(atom1.z)+'\n')
            for coord2 in tangents[n+1:]:
                atom2 = sAtom(name, coord2[0], coord2[1], coord2[2])
                dist.append(get_length(atom1, atom2, 'lengths'))

        if len(tangents) == 0:
            length = 0
        else:
            length = max(dist)-2.8  #testing 2.8
        f.write(str(round(length, 3)))
    """

    #with output
    dist = []
    dist_c = []
    for n, coord1 in enumerate(tangents):
        atom1 = sAtom(name, coord1[0], coord1[1], coord1[2])
        #(str(atom1.x)+' '+str(atom1.y)+' '+str(atom1.z)+'\n')
        for coord2 in tangents[n+1:]:
            atom2 = sAtom(name, coord2[0], coord2[1], coord2[2])
            dist.append(get_length(atom1, atom2))
            dist_c.append([atom1, atom2])

    if len(tangents) == 0:
        length = 0
        length_coord = ""
    else:
        #length = max(dist)-2.8  #two times water 1.4*2

        atom1, atom2 = dist_c[dist.index(max(dist))]
        #cmd.pseudoatom('water2', pos=atom2.xyz())
        for atom in (atom1, atom2):
            flat1_c = coord_reader(vector_segment(atom, get_nearest(atom, flat1, True)[2], 1.4))
            flat2_c = coord_reader(vector_segment(atom, get_nearest(atom, flat2, True)[2], 1.4))
            contact_p = coord_reader(vector_segment(atom, get_midpoint(flat1_c, flat2_c), 1.4))
            for xyz in ('x', 'y', 'z'):
                setattr(atom, xyz, getattr(contact_p, xyz))

        #cmd.pseudoatom('contact_p', pos=contact_p.xyz())
        #cmd.pseudoatom('flat1_c', pos=flat1_c.xyz())
        #cmd.pseudoatom('flat2_c', pos=flat2_c.xyz())
        length = get_length(atom1, atom2)

        #save coordinates
        setattr(atom1, axis, avg_height)
        setattr(atom2, axis, avg_height)
        cmd.pseudoatom('tmp1', pos=atom1.xyz())
        cmd.pseudoatom('tmp2', pos=atom2.xyz())
        cmd.set('suspend_updates', 0, quiet=1)
        if angle is not None:
            cmd.rotate(v_cross, angle=-angle, selection='tmp1 or tmp2', camera=0, object=None, origin=[0, 0, 0])

        atom1 = coord_reader(get_data('tmp1','[x, y, z]')[0], name)
        atom2 = coord_reader(get_data('tmp2','[x, y, z]')[0], name)
        length_coord = '%s_%s_%s_%s_%s_%s'%(round(atom1.x, 3), round(atom1.y, 3), round(atom1.z, 3),
                                            round(atom2.x, 3), round(atom2.y, 3), round(atom2.z, 3))

        #drawing a line
        cmd.distance('lengths','tmp1', 'tmp2')
        cmd.delete('tmp1')
        cmd.delete('tmp2')
        cmd.hide('label', 'lengths')

    return length, height, length_coord


def class_identifier(count):
    """Identifies the class of the interface"""

    modelA = pdb_reader('Sheet_0.pdb')
    modelB = pdb_reader('Sheet_%s.pdb'% str(count))
    axis = find_axis(modelA)

    if len(modelA.splice()) == 3:
        model1 = modelA.splice()[1]
        model2 = modelB.splice()[1]

        #finds the endpoints of the peptide N, C
        try:
            Nterm1 = model1.residues()[0].select('name', 'CA').lista[0]
        except IndexError:
            Nterm1 = model1.residues()[0].lista[0]
        try:
            Cterm1 = model1.residues()[-1].select('name', 'CA').lista[-1]
        except IndexError:
            Cterm1 = model1.residues()[-1].lista[-1]

        #the atom with the largest distance from the N-C line
        Side1 = get_mdistline(model1, get_line(Nterm1, Cterm1, axis), axis)[-1][2]
        #get_length(Nterm1, Cterm1, show='a')
        #get_length(Nterm1, Side1, show='a')

        #same for the other chain
        Nterm2 = model2.select('resi', [Nterm1.resi]).select('name', [Nterm1.name]).lista[0]
        Cterm2 = model2.select('resi', [Cterm1.resi]).select('name', [Cterm1.name]).lista[0]
        Side2 = model2.select('resi', [Side1.resi]).select('name', [Side1.name]).lista[0]
        #get_length(Nterm2, Cterm2, show='a')
        #get_length(Nterm2, Side2, show='a')

        v1, d1 = triangulate(Nterm1, Cterm1, Side1, axis)
        v2, d2 = triangulate(Nterm2, Cterm2, Side2, axis)

        #test for NtoN then FacetoFace

        if v1[0]*v2[0]+v1[1]*v2[1] > 0:
            if d1*d2 > 0:
                c = '2'
            else:
                c = '3'
        else:
            if d1*d2 > 0:
                c = '1'
            else:
                c = '4'

    elif len(modelA.splice()) == 6:
        model1A = modelA.splice()[2]
        model1B = modelA.splice()[3]
        model2A = modelB.splice()[2]
        model2B = modelB.splice()[3]

        #finds the endpoints of the peptide N, C
        try:
            Nterm1A = model1A.residues()[0].select('name', 'CA').lista[0]
        except IndexError:
            Nterm1A = model1A.residues()[0].lista[0]
        try:
            Cterm1A = model1A.residues()[-1].select('name', 'CA').lista[0]
        except IndexError:
            Cterm1A = model1A.residues()[-1].lista[0]

        #the atom with the largest distance from the N-C line
        Side1A = get_mdistline(model1A, get_line(Nterm1A, Cterm1A, axis), axis)[-1][2]
        #get_length(Nterm1A, Cterm1A, show='a')
        #get_length(Nterm1A, Side1A, show='a')

        #same for the other 3 chain
        Nterm1B = model1B.select('resi', [Nterm1A.resi]).select('name', [Nterm1A.name]).lista[0]
        Cterm1B = model1B.select('resi', [Cterm1A.resi]).select('name', [Cterm1A.name]).lista[0]
        Side1B = model1B.select('resi', [Side1A.resi]).select('name', [Side1A.name]).lista[0]
        #get_length(Nterm1B, Cterm1B, show='a')
        #get_length(Nterm1B, Side1B, show='a')

        Nterm2A = model2A.select('resi', [Nterm1A.resi]).select('name', [Nterm1A.name]).lista[0]
        Cterm2A = model2A.select('resi', [Cterm1A.resi]).select('name', [Cterm1A.name]).lista[0]
        Side2A = model2A.select('resi', [Side1A.resi]).select('name', [Side1A.name]).lista[0]
        #get_length(Nterm2A, Cterm2A, show='bb')
        #get_length(Nterm2A, Side2A, show='bb')

        Nterm2B = model2B.select('resi', [Nterm1A.resi]).select('name', [Nterm1A.name]).lista[0]
        Cterm2B = model2B.select('resi', [Cterm1A.resi]).select('name', [Cterm1A.name]).lista[0]
        Side2B = model2B.select('resi', [Side1A.resi]).select('name', [Side1A.name]).lista[0]
        #get_length(Nterm2B, Cterm2B, show='bb')
        #get_length(Nterm2B, Side2B, show='bb')

        v1A, d1A = triangulate(Nterm1A, Cterm1A, Side1A, axis)
        v1B, d1B = triangulate(Nterm1B, Cterm1B, Side1B, axis)
        v2A, d2A = triangulate(Nterm2A, Cterm2A, Side2A, axis)
        v2B, d2B = triangulate(Nterm2B, Cterm2B, Side2B, axis)

        #test for parallel then NtoN then FacetoFace
        if v1A[0]*v1B[0]+v1A[1]*v1B[1] > 0:
            if d1A*d1B > 0:
                if v1A[0]*v2A[0]+v1A[1]*v2A[1] > 0:
                    if d1A*d2A > 0:
                        c = '2'
                    else:
                        c = '3'
                else:
                    if d1A*d2A > 0:
                        c = '1'
                    else:
                        c = '4'
            else:
                if d1A*d2A > 0:
                    c = '9'
                else:
                    c = '10'
        else:
            if d1A*d1B > 0:
                if d1A*d2A > 0:
                    c = '7'
                else:
                    c = '8'
            else:
                if (d1A*d2A)*(v1A[0]*v2A[0]+v1A[1]*v2A[1]) > 0:
                    c = '6'
                else:
                    c = '5'
    return c


def atomic_area_fast(count, name):
    """Calculates atomic area values from the output of areaimol, and sum it into contributions."""

    A = pdb_reader('areaimol_%s_sheet_%s_atomic.pdb'%(name, 0))
    B = pdb_reader('areaimol_%s_sheet_%s_atomic.pdb'%(name, count))
    AB = pdb_reader('areaimol_%s_sheet_0_%s_atomic.pdb'%(name, count))

    chain_count = len(A.splice())

    if chain_count == 3:
        midA = A.select('chain', ['B'])
        midB = B.select('chain', ['B'])
        midAB = AB.select('chain', ['B', 'E'])
        chains = ['B', 'E']
    else:
        midA = A.select('chain', ['C', 'D'])
        midB = B.select('chain', ['C', 'D'])
        midAB = AB.select('chain', ['C', 'D', 'I', 'J'])
        chains = ['C', 'D', 'I', 'J']

    #calculta area in place of B factor
    for atomA in midA.lista:
        for atomB in midAB.lista:
            if is_same(atomA, atomB):
                atomB.b = str(round(float(atomA.b)-float(atomB.b), 1))
    for atomA in midB.lista:
        for atomB in midAB.lista:
            if is_same(atomA, atomB):
                atomB.b = str(round(float(atomA.b)-float(atomB.b), 1))

    #output
    #list format
    #chain = {'B':'1', 'C':'1', 'D':'2', 'E':'2', 'I':'3', 'J':'4'}
    contribution = []
    b_sum = 0
    fingerprint = [] #interacting sidechainlist
    for chainpair in (['E', 'I', 'J'], ['B', 'C', 'D']):
        sheetprint = []
        sheet_contr = []
        for chain in chainpair:
            if chain not in chains:
                continue
            chainprint = ""
            chain_contr = []
            for res in midAB.select('chain', chain).residues():

                ba = 0
                for atom in res.backbone(False).lista:
                    ba += float(atom.b)
                    b_sum += float(atom.b)

                bb = 0
                for atom in res.backbone(True).lista:
                    bb += float(atom.b)
                    b_sum += float(atom.b)

                if res.lista[0].resn == 'GLY':
                    b = float(res.select('name', ['CA']).lista[0].b)
                    ba += b
                    bb -= b

                if ba > 10 or (ba > 5 and res.lista[0].resn == 'GLY'): #size limit
                    chainprint += res.lista[0].resi+lettercodes[res.lista[0].resn]

                #complete contribution: residue number, chain name, 1letter name, area aa, area bb
                chain_contr.append([res.lista[0].resi, res.lista[0].chain, lettercodes[res.lista[0].resn], str(round(ba, 3)), str(round(bb, 3)), str(round(ba+bb, 3))])
            sheetprint.append(chainprint)
            sheet_contr.append([chainprint]+chain_contr)
        sheetprint.sort()
        fingerprint.append(','.join(sheetprint))
        sheet_contr.sort()
        for n in range(len(sheet_contr)):
            contribution.append(sheet_contr[n][1:])

    contribution_output = []
    if fingerprint[0] < fingerprint[1]:
        order = 0
        if len(contribution) == 2:
            contribution_output = contribution[0] + contribution[1]
        else:
            contribution_output = contribution[0] + contribution[1] + contribution[2] + contribution[3]
    else:
        fingerprint.sort()
        order = 1
        if len(contribution) == 2:
            contribution_output = contribution[1] + contribution[0]
        else:
            contribution_output = contribution[2] + contribution[3] + contribution[0] + contribution[1]

    fingerprint = tuple(fingerprint)

    return contribution_output, fingerprint, order


def bsheet_editor(name):
    """Edits the chain letter in Bsheet files."""

    with open("Bsheet_%s.pdb"%name, 'r+') as f:

        sheet = f.read().splitlines()
        res = int(sheet[0][22:26].strip())
        n = 0
        f.seek(0)
        for line in sheet:
            if line[0:6] == 'ATOM  ' or line[0:6] == 'HETATM' or line[0:6] == 'ANISOU':

                if int(line[22:26].strip())<res:
                    n +=1
                    f.write('TER\n')
                res = int(line[22:26].strip())

                f.write(line[0:21]+ch_n(n)+line[22:])
                f.write('\n')

        f.write('TER\nEND')



def areaimol_input(count, name):
    """Writes files as seperate files, chains to use for areaimol."""

    with open("%s_slicedsheet_%s.pdb"%(name, count), 'w+') as f:

        chains = [chain.center() for chain in pdb_reader("Sheet_%s.pdb"%count).splice(sort=False)]
        chains_ord = [chain.center() for chain in pdb_reader("Sheet_%s.pdb"%count).splice(sort=True)]
        chain_count = len(chains)

        with open('Sheet_%s.pdb'%count, 'r+') as g:
            pdb = g.read().splitlines()
            resi = 0
            last = 0
            pdbout = []
            for cur, line in enumerate(pdb):
                if not (line[0:6] == 'ATOM  ' or line[0:6] == 'HETATM'):
                    continue
                if int(line[22:26].strip()) < resi:
                    resi = 0
                    pdbout.append(pdb[last:cur])
                    last = cur
                resi = int(line[22:26].strip())
            pdbout.append(pdb[last:cur])

        pdbout = [tup[0] for chain in chains_ord for tup in zip(pdbout, chains) if chain == tup[1]]

        for n, chain in enumerate(pdbout):
            for x, line in enumerate(chain):
                if line[0:6] == 'ATOM  ' or line[0:6] == 'HETATM':
                    f.write(''.join([line[0:6], str(x+1).rjust(5), line[11:21], ch_n(n), line[22:], '\n']))
            f.write('TER\n')
        f.write('END')

    if count > 0:
        with open("%s_slicedsheet_0_%s.pdb"%(name, count), 'w+') as f:

            with open("%s_slicedsheet_%s.pdb"%(name, count), 'r+') as g:
                for line in g.read().splitlines()[:-1]:
                    f.write(line+'\n')

            with open("%s_slicedsheet_0.pdb"%name, 'r+') as g:
                n = 0
                for line in g.read().splitlines()[:-1]:
                    if line[0:6] == 'ATOM  ' or line[0:6] == 'HETATM':
                        f.write(''.join([line[0:21], ch_n(n+chain_count), line[22:], '\n']))
                    else:
                        n += 1
                        f.write(line+'\n')
            f.write('END')


def areaimol_calc(count, name, path, atomic=False):
    """write and execute Areaimol calculation via CCP4 trough subprocess."""
    #old name: areaimol_calc_3

    if atomic:
        output = 'OUTPUT\n'
    else:
        output = ''

    if platform.system() == 'Windows':
        p = subprocess.Popen(['cmd.exe'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True, shell=True)
        p.stdin.write('%s/ccp4.setup.bat\n'%path)
        result = p.communicate('%s/bin/areaimol.exe XYZIN %s_slicedsheet_%s.pdb > areaimol_%s_sheet_%s.txt\n%sEND\n'%(path, name, count, name, count, output))
    elif platform.system() == 'Linux':
        p = subprocess.Popen(['/bin/sh'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
        result = p.communicate('%s/bin/areaimol XYZIN %s_slicedsheet_%s.pdb > areaimol_%s_sheet_%s.txt\n%sEND\n'%(path, name, count, name, count, output))
    if atomic:
        shutil.move('XYZOUT', 'areaimol_%s_sheet_%s_atomic.pdb'%(name, count))

    if count == 0:
        pass
    else:
        if platform.system() == 'Windows':
            p = subprocess.Popen(['cmd.exe'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True, shell=True)
            p.stdin.write('%s/ccp4.setup.bat\n'%path)
            result = p.communicate('%s/bin/areaimol.exe XYZIN %s_slicedsheet_0_%s.pdb > areaimol_%s_sheet_0_%s.txt\n%sEND\n'%(path, name, count, name, count, output))
        elif platform.system() == 'Linux':
            p = subprocess.Popen(['/bin/sh'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
            result = p.communicate('%s/bin/areaimol XYZIN %s_slicedsheet_0_%s.pdb > areaimol_%s_sheet_0_%s.txt\n%sEND\n'%(path, name, count, name, count, output))
        if atomic:
            shutil.move('XYZOUT', 'areaimol_%s_sheet_0_%s_atomic.pdb'%(name, count))


def areaimol_out(count, name, sheet_number):
    """Read the output file of areaimol to get area."""
    #old name: areaimol_out_1
    #specific middle chains from 3 and 6 long sheets

    area1 = -1
    area1X = -1
    areaX = -1
    carea1 = -1
    carea1X = -1
    careaX = -1

    if sheet_number == 3:
        #parallel

        with open("areaimol_%s_sheet_%s.txt"%(name, count), 'r+') as f:
            for line in f.read().splitlines():
                if line[0:21] == 'Total area of chain B':
                    areaX = float(line.split()[-1])
                if line[0:29] == 'Total contact area of chain B':
                    careaX = float(line.split()[-1])

        with open("areaimol_%s_sheet_0_%s.txt"%(name, count), 'r+') as f:
            for line in f.read().splitlines():
                if line[0:21] == 'Total area of chain B':
                    area1X = float(line.split()[-1])
                if line[0:29] == 'Total contact area of chain B':
                    carea1X = float(line.split()[-1])
                if line[0:21] == 'Total area of chain E':
                    area1X = area1X + float(line.split()[-1])
                if line[0:29] == 'Total contact area of chain E':
                    carea1X = carea1X + float(line.split()[-1])

        with open("areaimol_%s_sheet_0.txt"%(name), 'r+') as f:
            for line in f.read().splitlines():
                if line[0:21] == 'Total area of chain B':
                    area1 = float(line.split()[-1])
                if line[0:29] == 'Total contact area of chain B':
                    carea1 = float(line.split()[-1])

    else:
        #anti

        with open("areaimol_%s_sheet_%s.txt"%(name, count), 'r+') as f:
            for line in f.read().splitlines():
                if line[0:21] == 'Total area of chain C':
                    areaX = float(line.split()[-1])
                if line[0:29] == 'Total contact area of chain C':
                    careaX = float(line.split()[-1])

                if line[0:21] == 'Total area of chain D':
                    areaX = areaX + float(line.split()[-1])
                if line[0:29] == 'Total contact area of chain D':
                    careaX = careaX + float(line.split()[-1])

        with open("areaimol_%s_sheet_0_%s.txt"%(name, count), 'r+') as f:
            for line in f.read().splitlines():
                if line[0:21] == 'Total area of chain C':
                    area1X = float(line.split()[-1])
                if line[0:29] == 'Total contact area of chain C':
                    carea1X = float(line.split()[-1])

                if line[0:21] == 'Total area of chain D':
                    area1X = area1X + float(line.split()[-1])
                if line[0:29] == 'Total contact area of chain D':
                    carea1X = carea1X + float(line.split()[-1])

                if line[0:21] == 'Total area of chain I':
                    area1X = area1X + float(line.split()[-1])
                if line[0:29] == 'Total contact area of chain I':
                    carea1X = carea1X + float(line.split()[-1])

                if line[0:21] == 'Total area of chain J':
                    area1X = area1X + float(line.split()[-1])
                if line[0:29] == 'Total contact area of chain J':
                    carea1X = carea1X + float(line.split()[-1])

        with open("areaimol_%s_sheet_0.txt"%(name), 'r+') as f:
            for line in f.read().splitlines():
                if line[0:21] == 'Total area of chain C':
                    area1 = float(line.split()[-1])
                if line[0:29] == 'Total contact area of chain C':
                    carea1 = float(line.split()[-1])

                if line[0:21] == 'Total area of chain D':
                    area1 = area1 + float(line.split()[-1])
                if line[0:29] == 'Total contact area of chain D':
                    carea1 = carea1 + float(line.split()[-1])

    return (area1, area1X, areaX, carea1, carea1X, careaX)



def clash_detector(count):
    """Detects if atoms are too close to each other."""

    clash = False
    a_radii = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80}

    model1 = pdb_reader('Sheet_0.pdb').splice_mid()
    model2 = pdb_reader('Sheet_%s.pdb'% str(count))
    for a1 in model1.lista:
        for a2 in model2.lista:
            if get_length(a1, a2) < a_radii[a1.elem] + a_radii[a2.elem] - 1:
                clash = True
                return clash

    model1 = pdb_reader('Sheet_0.pdb')
    model2.splice_mid()
    for a1 in model1.lista:
        for a2 in model2.lista:
            if get_length(a1, a2) < a_radii[a1.elem] + a_radii[a2.elem] - 1:
                clash = True
                return clash

    return clash
