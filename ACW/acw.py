# -*- coding: utf-8 -*-
"""
Created on Thu Dec 30 18:15:56 2021
Last edited on Mon Dec 04 12:14:30 2023
Final pre-release version
@author: mate
"""

import os
import sys
import glob
import shutil
import threading
import time
import inspect
import platform
from pymol import cmd
from pymol.wizard import Wizard

# adding install folder to the system path
acw_database_path = os.path.dirname(inspect.getfile(inspect.currentframe()))
if acw_database_path == '':
    acw_database_path = os.getcwd()
sys.path.insert(0, acw_database_path)
from ACW import *


class AcWizardMultiSheet(Wizard):
    """Pymol wizard for the evaluation and display of amyloid peptide structures.

    The wizard offers two mode:
    Evaluation mode to analyze a new amyloid peptide structure.
    Database mode to display PDB and personally saved structures.
    """
    def __init__(self, location):
        Wizard.__init__(self)

        self.main_path = location
        try:

            cmd.cd(self.main_path)
        except Exception as e:
            print('Installation folder is missing')
            raise OSError

        self.advanced = False
        self.re_calc = False
        self.re_calc_all = False

        try:
            self.pdb_list = read_save('acw_personal_database.json')
            print('Personal ACW database is now loaded')
        except Exception:
            try:
                self.pdb_list = read_save('acw_database.json')
                print('Original ACW database is now loaded')
            except Exception:
                print('Database files are missing or corrupted')
                self.pdb_list = amy_struct_list()

        self.path = None
        config_reader(self)
        if self.path:
            if not glob.glob('%s/bin/areaimol*'%self.path):
                self.path = None
                print('Areaimol not found at CCP4 path, the evaluation mode is not available')
            elif not glob.glob('%s/bin/sc*'%self.path):
                self.path = None
                print('Sc not found at CCP4 path, the evaluation mode is not available')
            if self.path:
                print('CCP4 path found, evaluation mode is available')
        else:
            print('CCP4 not detected, the evaluation mode is not available')

        # Pymol settings
        cmd.set('mouse_selection_mode', 4) # set selection mode to object
        cmd.set('orthoscopic')
        cmd.set('label_size', '20')
        cmd.set('dash_color', 'brightorange')
        try:
            cmd.undo_disable()
        except Exception:
            cmd.set('suspend_undo', 1)

        # menu buttons for advanced
        self.sec_class_types = ['a', 'b', '1', '2', '3', '4', '5', '6', '7', '8', 'O', 'L']
        sec_class = [[2, '2nd classes', '']]
        for c in self.sec_class_types:
            sec_class.append([1, c, 'cmd.get_wizard().set_value("sec_class","%s","amyloid")'%c])
        self.menu['sec_class'] = sec_class

        self.packing_types = ['1D', '2D', '3D']
        packing = [[2, 'Packings', '']]
        for p in self.packing_types:
            packing.append([1, p, 'cmd.get_wizard().set_value("packing","%s","amyloid")'%p])
        self.menu['packing'] = packing

        self.picking_mode_types = ['Select', 'ON', 'Color M', 'Color Y', 'Delete', 'OFF']

        # selection criteria possible values
        self.sele_clas_types = [str(a) for a in range(1, 9)] + ['L', 'O']
        self.sele_secclas_types = [str(a) for a in range(1, 9)] + ['-']
        self.sele_seq_len_types = [str(a) for a in range(3, 12)]
        self.sele_packing_types = ['1D', '2D', '3D']
        self.sele_source_types = ['personal', 'PDB']
        self.sele_aa_types = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
        # selection criteria selection containers
        self.sele_clas = []
        self.sele_secclas = []
        self.sele_seq_len = []
        self.sele_packing = []
        self.sele_source = []
        self.sele_aa = []

        # inital wizard parameters
        self.picking_mode = 'OFF'
        self.need_prompt = 'ON'
        self.labels = 'OFF'
        self.axis = 'OFF'
        self.error_message = ''
        self.busy = 0
        self.sheet_count = -1 # number of chain in sheet
        self.pick_count = 0 # number of sheet picked
        self.interface_count = 0 # which interface to display
        self.view_count = 0
        self.n_Msheet = 0 # sheet counters for automatic
        self.n_Bsheet = 0
        self.favor_param = None
        self.seq = ''
        self.seq_len = 0
        self.main_axis = 0
        self.overwrite_state = 0

        # set up evaluation state
        self.mode = 'evaluate'
        self.calc_state = [0]
        self.main_object = None
        self.pdb_count = -1
        self.pdb_list.add_blank()
        self.menu_update()

        cmd.refresh_wizard()
        if not self.path:
            self.mode_change('database')
            if self.mode == 'evaluate':
                print('ACW could not start')
                raise OSError


    def mode_change(self, new_mode):
        """Changes between the two main modes."""

        if new_mode == 'database':
            if self.busy or self.calc_state[0] in [2, 5]:
                print('Wait for completion')
                self.error_message = 'Wait for completion'
                return
            if len(self.pdb_list.lista) == 1:
                print('Database is missing, database mode is not available')
                self.error_message = 'Database is missing, database mode is not available'
                return
            self.picking_mode = 'OFF'
            cmd.delete('all')
            self.delete_all(['everything'])
            self.reset_counters()
            self.mode = new_mode
            try:
                cmd.cd(self.main_path)
                cmd.cd('%s'%self.pdb_list[0].code)
            except Exception:
                print('Database is missing, database mode is not available')
                self.error_message = 'Database is missing, database mode is not available'
                self.menu_update()
                self.mode = 'evaluate'
                return
            self.pdb_count = 0
            self.pdb_list.lista.pop(-1)
            self.reset_selection()
            self.menu_update()

        elif new_mode == 'evaluate':
            if not self.path:
                print('Without CCP4, evaluation mode is not available')
                self.error_message = 'Without CCP4, evaluation mode is not available'
                return

            self.picking_mode = 'OFF'
            cmd.delete('all')
            self.delete_all(['everything'])
            cmd.cd(self.main_path)
            write_save('acw_personal_database.json', self.pdb_list)
            self.reset_counters()
            self.pdb_list.add_blank()
            self.pdb_count = -1
            self.calc_state = [0]
            self.mode = new_mode
            self.menu_update()


    def next_peptide(self, n=1, x=None, goto=None):
        """Cleans up and loads in next peptide from the selection"""

        if self.mode == 'evaluate':
            cmd.get_wizard().mode_change("database")

        cmd.set('suspend_updates', 1, quiet=1)
        cmd.delete('all')
        self.delete_all(['everything'])
        self.reset_counters()
        self.picking_mode = 'OFF'
        cmd.set('mouse_selection_mode', 4) # set selection mode to object

        cmd.cd(self.main_path)
        write_save('acw_personal_database.json', self.pdb_list)
        cmd.cd('%s'%self.cur_pdb().code)

        while True:
            if x is not None:
                self.pdb_count = x
                break
            if goto is not None:
                if isinstance(goto, int):
                    print('To go to a given structure write its pdb code or its number in the database')
                    self.error_message = 'To go to a given structure write its pdb code or its number in the database'
                    self.pdb_count = 0
                else:
                    for pos, amy in enumerate(self.pdb_list.lista):
                        if amy.code == goto:
                            self.pdb_count = pos
                            break
                    else:
                        try:
                            goto = int(goto)-1
                            if 0 <= goto < len(self.pdb_list.lista):
                                self.pdb_count = goto
                            else:
                                print('Given  number or PDB-code not in range, going to first structure')
                                self.error_message = 'Given  number or PDB-code not in range, going to first structure'
                                self.pdb_count = 0
                        except ValueError:
                            print('File name not found, going to first')
                            self.error_message = 'File name not found, going to first'
                            self.pdb_count = 0
                    break

            # check next structure
            self.pdb_count += n

            if self.pdb_count == -2:
                self.pdb_count = len(self.pdb_list.lista)-1

            if self.pdb_count in (len(self.pdb_list.lista), -1):
                self.pick_count = -2
                self.pdb_count = -1
                cmd.cd('../%s'%self.pdb_list[0].code)
                cmd.set('suspend_updates', 0, quiet=1)
                print('No more structures.\nPress next or previous peptide to continue')
                cmd.refresh_wizard()
                return

            # criteria to open:
            selection = True
            if self.sele_clas:
                if self.cur_pdb().clas not in self.sele_clas:
                    selection = False
            if self.sele_secclas:
                if self.cur_pdb().sec_class not in self.sele_secclas:
                    selection = False
            if self.sele_seq_len:
                if self.cur_pdb().seq_len not in self.sele_seq_len:
                    selection = False
            if self.sele_packing:
                if self.cur_pdb().packing not in self.sele_packing:
                    selection = False
            if self.sele_source:
                if self.cur_pdb().source not in self.sele_source:
                    selection = False
            if self.sele_aa:
                if not all(aa in self.cur_pdb().seq for aa in self.sele_aa):
                    selection = False
            if selection:
                break

        # change directory
        cmd.cd('..')
        if not os.path.isdir(self.cur_pdb().code):
            os.makedirs(self.cur_pdb().code)
        cmd.cd(self.cur_pdb().code)

        if self.re_calc:
            cmd.load('%s.pdb'%self.cur_pdb().code)
            self.main_object = cmd.get_names()[0]
            self.setup()
            self.automatic_calc()
            cmd.set('suspend_updates', 0, quiet=1)
            self.reset()
            if self.re_calc_all:
                cmd.get_wizard().next_peptide(1)

        else:
            try:
                self.seq_len = int(self.cur_pdb().seq_len)
                self.read_in()
                if len(self.cur_pdb().interface) > 0:
                    self.next_interface()
                self.menu_update()
                self.toggle_labels()

            except Exception:
                print('Entry is missing files or corrupted')
                self.error_message = 'Entry is missing files or corrupted'
                return

            finally:
                cmd.set('suspend_updates', 0, quiet=1)
                self.reset()


    def next_interface(self, verbose=True):
        """Displays interfaces and their descriptors."""
        if len(self.cur_pdb().interface) == 0 or self.pick_count == -2:
            print('No interface to display')
            self.error_message = 'No interface to display'
            return

        if self.picking_mode == 'ON':
            print('Not available in manual mode')
            self.error_message = 'Not available in manual mode'
            return

        if self.interface_count == len(self.cur_pdb().interface):
            self.interface_count = 0

        # print information
        if verbose:
            print('')
            print('Interface %s'%(self.interface_count+1))
            sheets = ''
            for pairs in self.cur_iface().sheets.split('+'):
                sheets += '+'.join(pairs.split('_'))+', '
            print('Interacting sheet pair(s): '+sheets[:-2])
            print('Sidechains in the interface: '+self.cur_iface().composition)
            print('I.Class: %s, Sc: %s'%(self.cur_iface().i_class, self.cur_iface().sc))
            print('B.Area: %s, SDi: %s'%(self.cur_iface().areaimol_chain, self.cur_iface().perfi))
            #print(self.cur_iface().perfi_old)

            clash = self.cur_iface().clash.split('+')
            for n, c in enumerate(clash):
                if c == 'yes':
                    pair = self.cur_iface().sheets.split('+')[n].split('_')
                    print('Steric clash between %s and %s'%(pair[0], pair[1]))

        try:
            # recolor the interacting sheets
            cmd.color('green', 'all and elem C')
            cmd.delete('lengths')
            for pairs in self.cur_iface().sheets.split('+'):
                for pair in pairs.split('_'):
                    sheet = pdb_reader('Bsheet_%s.pdb'%pair).splice()
                    cmd.select('sele', 'none')
                    for chain in sheet:
                        x, y, z = chain.lista[0].xyz()
                        position = (x-0.005, x+0.005, y-0.005, y+0.005, z-0.005, z+0.005)
                        cmd.select('sele', 'sele or ( X>%s and X<%s and Y>%s and Y<%s and Z>%s and Z<%s )'%position, 0, 1, 0, 1)
                    cmd.select('sele', 'bychain sele', 0, 1, 0, 1)
                    cmd.color('lightmagenta', 'sele and polymer and elem C')

            # show lengths
            for coord_pair in self.cur_iface().length_coord.split('+'):

                coords = list(map(float, coord_pair.split('_')))
                #get_length_coord(coords[0:3], coords[3:6], show='lengths')

                #rectangle:    02
                #              13
                height = float(self.cur_iface().height)/2
                rectangle = [coords[0:3], coords[0:3], coords[3:6], coords[3:6]]
                rectangle[0][1] += height
                rectangle[1][1] -= height
                rectangle[2][1] += height
                rectangle[3][1] -= height

                get_length_coord(rectangle[0], rectangle[1], show='lengths')
                get_length_coord(rectangle[0], rectangle[2], show='lengths')
                get_length_coord(rectangle[3], rectangle[1], show='lengths')
                get_length_coord(rectangle[2], rectangle[3], show='lengths')

                cmd.hide('label', 'lengths')

        except Exception:
            pass

        self.interface_count += 1
        self.reset()


    def read_in(self):
        """Open Bsheet files from setup."""

        self.main_object = self.cur_pdb().code

        for n in range(int(self.cur_pdb().n_Bsheet)):
            cmd.load('Bsheet_%s.pdb'%(n+1))
        cmd.load('Nonprotein.pdb')

        cmd.color('green', 'elem C')
        cmd.hide('everything', 'polymer')
        cmd.show('wire', 'sidechain')
        cmd.show('stick', 'backbone')
        cmd.show('nb_spheres')
        self.pick_count = 0
        cmd.select('none')
        cmd.select('sele', 'Bsheet_1', 0, 1, 0, 1)

        # find main axis
        self.pick_sheet('sele_x', False)
        cmd.select('sele_x', 'sele_x and polymer', 0, 1, 0, 1)
        cmd.save('Sheet_X.pdb', 'sele_x')
        self.main_axis = find_axis(pdb_reader('Sheet_X.pdb'))
        cmd.select('sele', 'none')
        self.view_count = 0
        self.view()

        cmd.refresh_wizard()


    def setup(self):
        """Setups and finds the sheets from an open structure."""

        if self.re_calc:
            pass

        else:
            if self.mode == 'evaluate' and self.calc_state[0] in [4]:
                print("Error detected, press '5: Reset' to continue")
                self.error_message = "Error detected, '5: Reset' to continue"
                return

            if self.mode == 'evaluate' and self.calc_state[0] not in [0]:
                print('Not available after sheets were assigned')
                self.error_message = 'Not available after sheets were assigned'
                return

            try:
                self.main_object = cmd.get_names()[0]
            except IndexError:
                print('No structure detected')
                self.error_message = 'No structure detected'
                cmd.refresh_wizard()
                return

            try:
                self.pdb_list[self.main_object]
                print('This name is already in the database, please rename the input file')
                self.error_message = 'This name is already in the database, please rename the input file'
                cmd.refresh_wizard()
                return
            except KeyError:
                pass

            cmd.cd(self.main_path)
            if os.path.isdir(self.main_object):
                shutil.rmtree(self.main_object)

            try:
                os.makedirs(self.main_object)
            except OSError:
                print('Cannot create a folder')
                self.calc_state = [4, 'Cannot create a folder']
                cmd.refresh_wizard()
                return
            cmd.cd(self.main_object)
            print('Searching for sheets...')

            self.calc_state[0] = 5
            cmd.set('suspend_updates', 1, quiet=1)

        # remove hydrogens
        cmd.remove('elem H')
        # select & remove all non A altlocs
        cmd.remove('not (alt ""+A)')
        # remove the alt identifier
        cmd.alter('all', 'alt=""')

        # renumber chains if needed
        try:
            chains = list(dict.fromkeys(get_data('polymer', 'chain'))) #fails if there are unmodeled atoms
            for chain in chains:
                offset = int(get_data('first (chain %s and polymer)'%chain, 'resi')[0]) -1

                if offset == 0:
                    continue
                else:
                    cmd.alter('chain %s and polymer'%chain, 'resi=str(int(resi)-%s)' % str(offset))
        except IndexError:
            print('The model contains unmodeled residues in the sequence, please edit the PDB file')
            self.calc_state = [4, 'The model contains unmodeled residues in the sequence, please edit the PDB file']
            cmd.set('suspend_updates', 0, quiet=1)
            cmd.refresh_wizard()
            return

        # setup sheets
        sym = cmd.get_symmetry(self.main_object)
        if sym:
            if sym[0] * sym[1] * sym[2] > 1:
                cmd.symexp('sym', self.main_object, self.main_object, 2000)
        cmd.color('green', 'elem C')
        cmd.hide('everything', 'polymer')
        cmd.show('wire', 'sidechain')
        cmd.show('stick', 'backbone')
        cmd.show('nb_spheres')

        self.pick_count = 0

        cmd.select('none')
        cmd.select('sele', 'bychain first (%s and polymer)'%self.main_object, 0, 1, 0, 1)

        try:
            cmd.save('Sheet_X.pdb', 'sele and polymer')
            model = pdb_reader('Sheet_X.pdb')
            self.seq_len = len(model.residues())
            self.seq = "".join([lettercodes[a.lista[0].resn] for a in model.residues()])

            if str(self.seq_len) != model.residues()[-1].lista[0].resi:
                print('The numbering of the residues is not continuous, please edit the PDB file')
                self.calc_state = [4, 'The numbering of the residues is not continuous, please edit the PDB file']
                cmd.set('suspend_updates', 0, quiet=1)
                cmd.refresh_wizard()
                return

            # straigthen and rotate to main axis
            self.pick_sheet('sele_x', False)
            straigthen_pca('sele_x')
            self.main_axis = 'y'

            cmd.select('sele', 'none')
            self.view_count = 0
            self.view()

        except IndexError:
            print('No suitable sheets found, see README for more information')
            self.error_message = 'No suitable sheets found, see README for more information'
            self.calc_state = [0]
            cmd.set('suspend_updates', 0, quiet=1)
            cmd.refresh_wizard()
            return

        self.suspend_refresh()

        if self.labels == 'ON':
            self.toggle_labels()

        # find all sheets in asymetric unit
        cmd.select('sele_main', '%s and polymer'%self.main_object, 0, 1, 0, 1)
        cmd.color('green', 'all and elem C')
        cmd.color('lightmagenta', 'sele_main and elem C')
        n = 0
        x = 0
        while len(get_data('sele_main', 'chain')) > 0:
            n += 1
            x += 1

            cmd.select('sele', 'bychain first sele_main', 0, 1, 0, 1)
            self.suspend_refresh()
            try:
                self.pick_sheet('Bsheet_%s'%n, False, f_p=self.favor_param)
            except IndexError:
                cmd.select('sele_main', 'sele_main and not Bsheet_%s'%n, 0, 1, 0, 1)
                cmd.color('green', 'Bsheet_%s and elem C'%n)
                cmd.delete('Bsheet_%s'%n)
                n -= 1
                continue

            cmd.select('sele_main', 'sele_main and not Bsheet_%s'%n, 0, 1, 0, 1)
            cmd.color('green', 'Bsheet_%s and elem C'%n)
            cmd.save('Bsheet_%s.pdb'%n, 'Bsheet_%s'%n)
            bsheet_editor(n)

            if x == 500:
                print('Error in locating sheets')
                cmd.set('suspend_updates', 0, quiet=1)
                self.calc_state = [4, 'Error in locating sheets']
                cmd.refresh_wizard()
                return

        cmd.delete('sele_main')
        self.n_Msheet = n
        self.suspend_refresh()

        # find all sheets everywhere else
        cmd.select('sele_all', 'all and polymer', 0, 1, 0, 1)
        for x in range(1, n+1):
            cmd.select('sele_all', 'sele_all and not Bsheet_%s'%x, 0, 1, 0, 1)
        cmd.color('lightmagenta', 'sele_all and elem C')

        while len(get_data('sele_all', 'chain')) > 0:
            n += 1
            x += 1
            cmd.select('sele', 'bychain first sele_all', 0, 1, 0, 1)
            self.suspend_refresh()
            try:
                self.pick_sheet('Bsheet_%s'%n, False, f_p=self.favor_param)
            except IndexError:
                cmd.select('sele_all', 'sele_all and not Bsheet_%s'%n, 0, 1, 0, 1)
                cmd.color('green', 'Bsheet_%s and elem C'%n)
                cmd.delete('Bsheet_%s'%n)
                n -= 1
                continue
            cmd.select('sele_all', 'sele_all and not Bsheet_%s'%n, 0, 1, 0, 1)
            cmd.color('green', 'Bsheet_%s and elem C'%n)
            cmd.save('Bsheet_%s.pdb'%n, 'Bsheet_%s'%n)
            bsheet_editor(n)

            if x == 500:
                print('Error in locating sheets')
                self.calc_state = [4, 'Error in locating sheets']
                cmd.set('suspend_updates', 0, quiet=1)
                cmd.refresh_wizard()
                return

        cmd.color('green', 'all and elem C')
        cmd.delete('sele_all')
        self.n_Bsheet = n

        # save waters, ions, leftover chains...
        cmd.select('sele_nonprot', 'all', 0, 1, 0, 1)
        for x in range(self.n_Bsheet):
            cmd.select('sele_nonprot', 'sele_nonprot and not Bsheet_%s'%(x+1), 0, 1, 0, 1)
        cmd.save('Nonprotein.pdb', 'sele_nonprot')

        # save and print results
        if self.n_Msheet > 0 and self.n_Bsheet > 1:
            self.cur_pdb().code = self.main_object
            self.cur_pdb().seq = self.seq
            self.cur_pdb().seq_len = str(self.seq_len)
            self.cur_pdb().n_Msheet = self.n_Msheet
            self.cur_pdb().n_Bsheet = self.n_Bsheet
            cmd.delete('all')
            self.read_in()
            print('%s sheets were found'%self.n_Bsheet)
            self.toggle_labels()
            self.calc_state = [1]
            cmd.set('suspend_updates', 0, quiet=1)
            cmd.refresh_wizard()

        else:
            print('Not enough sheets were found for calculation')
            self.error_message = 'Not enough sheets were found for calculation'
            self.calc_state = [0]
            cmd.set('suspend_updates', 0, quiet=1)
            cmd.refresh_wizard()


    def reset(self):
        """Reset everything picking related."""
        cmd.delete('sele')
        cmd.select('none')
        cmd.unpick()
        cmd.refresh_wizard()


    def suspend_refresh(self, wait=0.05):
        """Wait to allow pymol to draw a frame."""
        cmd.set('suspend_updates', 0, quiet=1)
        time.sleep(wait)
        cmd.set('suspend_updates', 1, quiet=1)


    def reset_counters(self):
        """Resets all state and data variables"""
        self.sheet_count = -1 # number of chain in sheet
        self.pick_count = 0 # number of sheet picked
        self.interface_count = 0 # which interface to display
        self.view_count = 0
        self.labels = 'OFF'
        self.axis = 'OFF'
        self.n_Msheet = 0 # sheet counters for automatic
        self.n_Bsheet = 0
        self.favor_param = None
        self.seq = ''
        self.seq_len = 0
        self.main_axis = 0
        self.overwrite_state = 0
        cmd.refresh_wizard()


    def reset_selection(self):
        """Resets selection, also goes to first structure"""
        self.sele_clas = []
        self.sele_secclas = []
        self.sele_seq_len = []
        self.sele_packing = []
        self.sele_source = []
        self.sele_aa = []
        self.next_peptide(x=0)
        self.menu_update()
        cmd.refresh_wizard()


    def delete_all(self, to_rm=[]):
        """Delete sheets, files and interfaces"""

        if not self.re_calc:
            if (self.busy or self.calc_state[0] in [2, 5]) and 'everything' in to_rm :
                print('Not available during calculation')
                self.error_message = 'Not available during calculation'
                return

        cmd.deselect()
        cmd.delete('sele*')
        cmd.delete('dist*')
        self.pick_count = 0

        if 'everything' in to_rm:
            cmd.delete('all')
            self.reset_counters()

        if 'everything' in to_rm or 'files' in to_rm:
            for name in ('distance*', 'sc*', 'areaimol*', 'Sheet_*', 'sele*'):
                for file in glob.glob(name):
                    os.remove(file)
            if self.main_object:
                for file in glob.glob(self.main_object + '_*.pdb'):
                    os.remove(file)

        if 'interface' in to_rm:
            cmd.util.cbag('all')
            cmd.delete('lengths')
            self.cur_pdb().interface = []
            self.interface_count = 0

        if 'last_amy' in to_rm:
            self.pdb_list.lista.pop(-1)
            self.pdb_list.add_blank()

        if 'sheets' in to_rm: # unused, deletes sheets and folder
            for file in glob.glob('Bsheet_*'):
                os.remove(file)
            if self.mode == 'evaluate' and self.calc_state[0] > 0:
                cmd.cd(self.main_path)
                shutil.rmtree(self.main_object)

        if 'everything' in to_rm:
            if self.mode == 'evaluate' and self.calc_state[0] > 0:
                cmd.cd(self.main_path)
            self.calc_state = [0]
            self.main_object = None

        cmd.refresh_wizard()


    def cleanup(self):
        """Clean up, executed on quitting"""
        self.delete_all(['everything'])
        self.reset()
        cmd.delete('sele*')

        if self.mode == 'database':
            cmd.cd(self.main_path)
        elif self.mode == 'evaluate':
            self.pdb_list.lista.pop(-1)

        if len(self.pdb_list.lista) > 0:
            write_save('acw_personal_database.json', self.pdb_list)

        print('Cleanup finished')


    def get_prompt(self):
        """Text to be displayed at top left of display"""
        self.prompt = []

        if self.need_prompt == 'OFF':
            self.prompt = None

        elif self.picking_mode == 'ON':
            if self.pick_count == 0:
                self.prompt = ['Manual mode', 'Choose the first sheet']
            if self.pick_count == 1:
                self.prompt = ['Manual mode', 'Choose the second sheet']

        elif self.mode == 'database':
            if self.pick_count == -2:
                self.prompt = ['No more structures.', 'Press next or previous peptide to continue']

            else:
                self.prompt = ['%s %s %s/%s'%(self.cur_pdb().code, self.cur_pdb().seq, str(self.pdb_count+1), str(len(self.pdb_list.lista)))]
                self.prompt.append('Packing %s, Class: %s'%(str(self.cur_pdb().packing), str(self.cur_pdb().clas)))
                if self.cur_pdb().sec_class != '-':
                    self.prompt[-1] += '(%s)'%self.cur_pdb().sec_class
                if len(self.cur_pdb().interface) > 0:
                    self.prompt.append('Sidechains: '+self.cur_iface(-1).composition)
                    self.prompt.append('I.Class: %s, Sc: %s'%(self.cur_iface(-1).i_class, self.cur_iface(-1).sc))
                    self.prompt.append('B.A.: %s, SDi: %s'%(round(float(self.cur_iface(-1).areaimol_chain), 1), self.cur_iface(-1).perfi))

                    if 'yes' in self.cur_iface(-1).clash.split('+'):
                        self.prompt.append('Steric clash found in the interface(s)')

        elif self.mode == 'evaluate':
            if self.main_object:
                self.prompt.append(self.main_object)

            if self.calc_state[0] == 0: # starting ready for setup
                self.prompt.append('Open an amyloid crystal structure and set up the sheets for evaluation')

            elif self.calc_state[0] == 1: # after setup ready for auto calc
                self.prompt.append('%s sheets were found'%self.n_Bsheet)

            elif self.calc_state[0] == 2: # busy auto calc
                self.prompt.append('Calculating...')

            elif self.calc_state[0] == 3: # after successful auto calc
                self.prompt.append('%s interface(s) were found'%len(self.cur_pdb().interface))
                if len(self.cur_pdb().interface) > 0:
                    self.prompt.append('Sidechains: '+self.cur_iface(-1).composition)
                    self.prompt.append('Class: %s, Sc: %s'%(self.cur_iface(-1).i_class, self.cur_iface(-1).sc))
                    self.prompt.append('B.A.: %s, SDi: %s'%(self.cur_iface(-1).areaimol_chain, self.cur_iface(-1).perfi))
                    if 'yes' in self.cur_iface(-1).clash.split('+'):
                        self.prompt.append('Steric clash found in the interface(s)')

            elif self.calc_state[0] == 4: # error state
                self.prompt.append(self.calc_state[1])

            elif self.calc_state[0] == 5: # busy setup
                self.prompt.append('Searching for sheets...')

        if self.error_message:
            self.prompt.append(self.error_message)

        return self.prompt


    def toggle_labels(self):
        """Shows and hides the sheet labels"""
        if self.labels == 'ON':
            cmd.delete('Labels')
            self.labels = 'OFF'
        else:
            try:
                if self.pick_count == -2:
                    raise IndexError

                if self.mode == 'database':
                    n_labels = int(self.cur_pdb().n_Bsheet)
                else:
                    n_labels = int(self.n_Bsheet)

                for n in range(n_labels):
                    xyz = pdb_reader('Bsheet_%s.pdb'%(n+1)).splice()[-1].center()
                    cmd.pseudoatom('Labels', pos=xyz, name=str(n+1))
                    cmd.label('Labels and name %s'%(n+1), str(n+1))
                cmd.hide('wire', 'Labels')
            except Exception:
                print('No sheets were detected')
                self.error_message = 'No sheets were detected'

            self.labels = 'ON'


    def do_select(self, name):
        """This activates first when clicking"""
        if self.busy:
            return
        cmd.get_wizard().launcher("multi", "do_pick", 0)


    def do_pick(self, bondFlag):
        """Commands after clicking"""
        if bondFlag or self.busy:
            return

        self.busy = 1
        try:
            self.error_message = ''
            if self.picking_mode == 'ON':

                if self.pick_count == 0:
                    cmd.util.cbag('all')
                    cmd.delete('lengths')
                    try:
                        self.reference_sheet = cmd.get_names('objects', 0, 'sele')[0].split('_')[1]
                        self.prepare_sheet()
                    except IndexError:
                        print('No sheets were found in selection')
                        return
                    self.pick_count += 1
                    self.sheet_count = len(pdb_reader('Sheet_0.pdb').splice())

                    try:
                        areaimol_calc(self.pick_count-1, self.main_object, self.path, True)
                        print('New reference sheet is set')
                    except Exception:
                        print('Error during surface calculation')
                        self.error_message = 'Error during surface calculation'
                        return
                    self.reset()

                elif self.pick_count == 1:
                    try:
                        curent_sheet = cmd.get_names('objects', 0, 'sele')[0].split('_')[1]
                        if curent_sheet == self.reference_sheet:
                            print('Choose a different sheet')
                            self.reset()
                            return
                        self.prepare_sheet()
                    except IndexError:
                        print('No sheets were found in selection')
                        return

                    self.pick_count += 1

                    try:
                        sc, area1, area2 = sc_calc_out(self.pick_count-1, self.main_object, self.path)
                        if sc != 'NaN':
                            sc = round(float(sc), 2)
                        areaimol_calc(self.pick_count-1, self.main_object, self.path, True)
                    except Exception:
                        print('Error during calculation')
                        self.error_message = 'Error during calculation'
                        return
                    area1, area1x, areax, carea1, carea1x, careax = areaimol_out(self.pick_count-1, self.main_object, self.sheet_count)
                    contribution, fingerprint, order = atomic_area_fast(self.pick_count-1, self.main_object)
                    length, height, length_coord = perf_calc_fast(self.pick_count-1, self.reference_sheet, curent_sheet)
                    i_class = class_identifier(self.pick_count-1)

                    if clash_detector(1):
                        clash = 'yes'
                    else:
                        clash = '-'

                    areaimol_chain = round((area1+areax-area1x) / (2*self.sheet_count/3), 3)
                    if height == 0 or length == 0:
                        perfi = 'NaN'
                    else:
                        perfi = round(areaimol_chain/height/length, 2)
                    areaimol_chain = round(areaimol_chain, 1)
                    contr = '_'.join(contribution[0])
                    for res in contribution[1:]:
                        contr += '-' + '_'.join(res)

                    print('')
                    print('Sheets: %s+%s'%(self.reference_sheet, curent_sheet))
                    print('Sidechains: (%s)+(%s)'%(fingerprint[0], fingerprint[1]))
                    print('Sc: %s'%sc)
                    print('B.area: %s'%str(areaimol_chain))
                    print('SDi: ' + str(perfi))
                    print('Class: ' + i_class)
                    if clash == 'yes':
                        print('Steric clash found between sheets')

                    # save file
                    n_output = 1
                    while os.path.isfile('ACW_manual_output_%s_%s.txt'%(self.main_object, n_output)):
                        n_output += 1

                    with open('ACW_manual_output_%s_%s.txt'%(self.main_object, n_output), 'w') as f:
                        f.write('Interface %s\n'%(n_output))
                        f.write('Sheets: %s+%s\n'%(self.reference_sheet, curent_sheet))
                        f.write('Sidechains: (%s)+(%s)\n'%(fingerprint[0], fingerprint[1]))
                        f.write('Sc: %s\n'%str(sc))
                        f.write('B.Area: %s\n'%str(areaimol_chain))
                        f.write('SDi: %s\n'%str(perfi))
                        f.write('Class: %s\n'%i_class)
                        if clash == 'yes':
                            f.write('Steric clash found between sheets\n')
                        f.write('\n')
                        f.write('B.area details: residue, sidechain, backbone, total \n')
                        for n, c1 in enumerate(contr.split('+')):
                            for c2 in c1.split('-'):
                                c3 = c2.split('_')
                                f.write('%s %s %s %s %s\n'%(c3[0], c3[2], c3[3], c3[4], c3[5]))
                            f.write('\n')
                        f.write('\n')

                    self.error_message = 'Interface %s is saved'%n_output
                    self.delete_all(['files'])
                    self.reset()

                else:
                    self.delete_all(['files'])
                    self.reset()

            elif self.picking_mode in ['Color M', 'Color Y']:
                try:
                    #self.pick_sheet('to_color')
                    cmd.select('to_color', 'bychain sele', 0, 1, 0, 1)
                except IndexError:
                    print('No sheet found')
                    return

                for n in range(2, self.seq_len + 1, 2):
                    cmd.color('green', 'to_color and resi %s'%n)
                for n in range(1, self.seq_len + 1, 2):
                    cmd.color('brown', 'to_color and resi %s'%n)
                if self.picking_mode[-1] == 'M':
                    cmd.color('lightmagenta', 'to_color and backbone ')
                else:
                    cmd.color('yellow', 'to_color and backbone ')
                cmd.color('cyan', 'to_color and resi 1 and name N')
                cmd.color('blue', 'to_color and resi %s and (name C or name OXT or name O)'%self.seq_len)

            elif self.picking_mode == 'Delete':
                try:
                    self.pick_sheet('to_rm')
                    cmd.select('to_rm', 'bychain to_rm', 0, 1, 0, 1)
                    cmd.remove('to_rm')
                    #cmd.hide('everything', 'sele')

                except IndexError:
                    print('No sheet found')

            else:
                pass
        finally:
            self.busy = 0


    def pick_sheet(self, name, color=False, f_p=None):
        """expand selection to a whole sheet"""

        cmd.select('sele_click', 'sele', 0, 1, 0, 1) # saving to preserve original selection
        success = False

        # favored parameter if available []
        if f_p:
            p, d = f_p
            cmd.select('sele', 'sele_click and polymer', 0, 1, 0, 1)
            selector_text = ''.join(['(bychain resi ', p[0], '-', p[1],
                                     ' and (name O or name N) and polymer within ',
                                     d, ' of sele and resi ', p[0], '-', p[1],
                                     ' and (name O or name N)) and polymer'])
            atoms = cmd.count_atoms('sele')
            for n in range(12):
                cmd.select('sele', selector_text, 0, 1, 0, 1)
                atoms2 = cmd.count_atoms('sele')
                if atoms == atoms2:
                    break
                atoms = atoms2
            cmd.save('sele.pdb', 'sele')
            if len(pdb_reader('sele.pdb').splice(False)) in [3, 6]:
                success = True

        # combinations of settings for finding betasheet
        if not success:
            parameters = [['2', str(self.seq_len-1)], ['1', str(self.seq_len-1)],
                          ['3', str(self.seq_len-1)], ['2', str(self.seq_len-2)],
                          ['2', str(self.seq_len)], ['3', str(self.seq_len-2)]]
            for p in parameters:
                for d in ['2.9', '2.7', '2.8', '3.0', '3.2', '3.4']:
                    cmd.select('sele', 'sele_click and polymer', 0, 1, 0, 1)
                    selector_text = ''.join(['(bychain resi ', p[0], '-', p[1],
                                             ' and (name O or name N) and polymer within ',
                                             d, ' of sele and resi ', p[0], '-', p[1],
                                             ' and (name O or name N)) and polymer'])

                    #selector_text = ''.join(['bychain (resi ',p[0],'-',p[1],' and name O and polymer within ',d,' of sele and resi ',p[0],'-',p[1],' and name N and polymer)'])
                    #selector_text += ''.join([' or (resi ',p[0],'-',p[1],' and name N within ',d,' of sele and resi ',p[0],'-',p[1],' and name O and polymer)'])

                    atoms = cmd.count_atoms('sele')
                    for n in range(12):
                        cmd.select('sele', selector_text, 0, 1, 0, 1)
                        atoms2 = cmd.count_atoms('sele')
                        if atoms == atoms2:
                            break
                        atoms = atoms2
                    cmd.save('sele.pdb', 'sele')
                    if len(pdb_reader('sele.pdb').splice(False)) in [3, 6]:
                        self.favor_param = (p, d)
                        break
                else:
                    continue
                break
            else:
                cmd.select(name, 'sele_click', 0, 1, 0, 1)
                cmd.delete('sele_click')
                self.reset()
                raise IndexError

        if color:
            try:
                colors = ['lightmagenta', 'cyan', 'salmon', 'yellow', 'cyan', 'hydrogen', 'slate']
                cmd.color(colors[self.pick_count], 'sele and elem C')
            except IndexError:
                cmd.color('lightmagenta', 'sele and elem C')

        cmd.select(name, 'sele', 0, 1, 0, 1)
        cmd.delete('sele_click')
        self.reset()


    def prepare_sheet(self):
        """Select sheet, color it and save it."""

        self.pick_sheet('sele_%s'%self.pick_count, True)

        cmd.remove('sele_%s and not polymer'%self.pick_count)
        cmd.save('Sheet_%s.pdb'%self.pick_count, 'sele_%s'%self.pick_count)
        sc_input(self.pick_count, self.main_object)
        areaimol_input(self.pick_count, self.main_object)
        cmd.refresh_wizard()


    def exp_image(self):
        """Exports a png picture."""

        cmd.ray('3000', '3000')
        n_output = 1
        while os.path.isfile('%s_%s.png'%(self.main_object, n_output)):
            n_output += 1
        cmd.save('%s_%s.png'%(self.main_object, n_output))
        print('%s_%s.png saved'%(self.main_object, n_output))

        self.reset()


    def set_value(self, value_type, value, object_type):
        """Change a saved value to an other."""

        if object_type == 'wizard':
            setattr(self, value_type, value)
        elif object_type == 'amyloid':
            setattr(self.pdb_list[self.pdb_count], value_type, value)
        elif object_type == 'interface':
            setattr(self.pdb_list[self.pdb_count].interface[self.interface_count-1], value_type, value)

        self.menu_update()
        self.reset()


    def toggle_value(self, value_type, values):
        """Toggle a saved value in lists."""
        for value in values:
            if value in getattr(self, value_type + '_types'):
                if value in getattr(self, value_type):
                    getattr(self, value_type).remove(value)
                else:
                    getattr(self, value_type).append(value)

        self.menu_update()
        self.reset()


    def def_changer(self, setting):
        """Toggles default settings."""

        setattr(self, setting, {'ON':'OFF', 'OFF':'ON'}[getattr(self, setting)])

        self.menu_update()
        self.reset()


    def cur_pdb(self):
        """Current pdb structure."""
        return self.pdb_list[self.pdb_count]


    def cur_iface(self, i=0):
        """Current interface of pdb structure."""
        return self.pdb_list[self.pdb_count].interface[self.interface_count+i]


    def automatic_calc(self):
        """Calculates all possible interfaces of the asymetric unit."""
        if not self.picking_mode == 'OFF':
            self.error_message = 'Not available in manual mode'
        elif self.mode == 'evaluate':
            if self.calc_state[0] == 4:
                self.error_message = "Error detected, '5: Reset' to continue"
            elif self.calc_state[0] in [0]:
                self.error_message = "Not available, use '1: Set up the sheets' first"
            elif self.calc_state[0] not in [1]:
                self.error_message = 'Not available'
        elif not self.path:
            self.error_message = 'CCP4 is not available, calculations cannot be performed'

        if self.error_message:
            print(self.error_message)
            return

        cmd.set('suspend_updates', 1, quiet=1)
        cmd.get_wizard().delete_all(['interface'])

        # calculate all possible interfaces
        print('Searching for interfaces...')
        self.calc_state = [2]
        cmd.refresh_wizard()
        interfaces = {}
        variables = ['areaimol_chain', 'areaimol_c', 'sc', 'perfi', 'length', 'height']
        self.sheet_count = len(pdb_reader('Bsheet_1.pdb').splice())

        for S in range(1, self.n_Msheet+1):
            cmd.save('Sheet_0.pdb', 'Bsheet_%s'%S)
            sheet_0 = pdb_reader('Sheet_0.pdb').splice_mid()
            center_0 = coord_reader(sheet_0.center())
            dist_0 = get_nearest(center_0, sheet_0)[-1][0]
            sc_input(0, self.main_object)
            areaimol_input(0, self.main_object)

            try:
                areaimol_calc(0, self.main_object, self.path, True)
            except Exception:
                print('Error during surface calculation')
                cmd.set('suspend_updates', 0, quiet=1)
                self.calc_state = [4, 'Error during surface calculation']
                cmd.refresh_wizard()
                return

            for B in range(S+1, self.n_Bsheet+1):
                cmd.color('green', 'all and elem C')
                cmd.color('lightmagenta', 'Bsheet_%s and elem C'%S)
                cmd.color('lightmagenta', 'Bsheet_%s and elem C'%B)
                self.suspend_refresh()
                cmd.save('Sheet_0.pdb', 'Bsheet_%s'%S)
                cmd.save('Sheet_1.pdb', 'Bsheet_%s'%B)

                # test for too far away sheets
                sheet_1 = pdb_reader('Sheet_1.pdb')
                center_1 = coord_reader(sheet_1.center())
                if get_length(center_0, center_1) - get_nearest(center_1, sheet_1)[-1][0]-dist_0 > 6.5:
                    continue
                nearest = get_mnearest(sheet_0, sheet_1, fast=True)
                if nearest[0] > 6.5:
                    continue

                #surface area calculation
                areaimol_input(1, self.main_object)
                areaimol_calc(1, self.main_object, self.path, True)
                area1, area1x, areax, carea1, carea1x, careax = areaimol_out(1, self.main_object, self.sheet_count)
                areaimol_chain = round((area1+areax-area1x)/(2*self.sheet_count/3), 3)
                areaimol_c = round((carea1+careax-carea1x)/(2*self.sheet_count/3), 3)

                # surface area requirement
                if areaimol_chain < 50:
                    continue

                contribution, fingerprint, order = atomic_area_fast(1, self.main_object)
                # number of amino acid requirement
                fingerprint_test = fingerprint[0].split(',') + fingerprint[1].split(',')
                if sum(list(map(len, fingerprint_test)))/len(fingerprint_test) < 3: # 1+2 interaction
                    continue
                if not all(len(fp) >= 2 for fp in fingerprint_test): # each chain must contact
                    continue

                print('Interface found between %s and %s'%(str(S), str(B)))
                if clash_detector(1):
                    print('Clash detcted between sheets %s and %s'%(S, B))
                    clash = 'yes'
                else:
                    clash = '-'

                i_class = class_identifier(1)
                fingerprint = tuple(list(fingerprint)+[i_class])

                sc_input(1, self.main_object)
                sc, area1, area2 = sc_calc_out(1, self.main_object, self.path)
                sc = float(sc)

                length, height, length_coord = perf_calc_fast(1, S, B)
                if length == 0:
                    perfi = 0
                else:
                    perfi = round(areaimol_chain/height/length, 3)

                # save interface to dictionary based on fingerprint
                try:
                    if order:
                        interfaces[fingerprint]['sheets'].append(str(S)+'_'+str(B))
                    else:
                        interfaces[fingerprint]['sheets'].append(str(B)+'_'+str(S))
                    for var in variables:
                        interfaces[fingerprint][var].append(vars()[var])
                    interfaces[fingerprint]['contribution'].append(contribution[:])
                    interfaces[fingerprint]['length_coord'].append(length_coord)
                    interfaces[fingerprint]['i_class'].append(i_class)
                    interfaces[fingerprint]['clash'].append(clash)
                except KeyError:
                    if order:
                        interfaces[fingerprint] = {'sheets':[str(S)+'_'+str(B)]}
                    else:
                        interfaces[fingerprint] = {'sheets':[str(B)+'_'+str(S)]}
                    for var in variables:
                        interfaces[fingerprint][var] = [vars()[var]]
                    interfaces[fingerprint]['contribution'] = [contribution[:]]
                    interfaces[fingerprint]['length_coord'] = [length_coord]
                    interfaces[fingerprint]['i_class'] = [i_class]
                    interfaces[fingerprint]['clash'] = [clash]

        for n, fp in enumerate(interfaces):
            # saving all data from dictionary
            # txt joined together
            self.cur_pdb().add_interface(amy_interface())
            self.cur_pdb().interface[-1].sheets = '+'.join(interfaces[fp]['sheets'])
            self.cur_pdb().interface[-1].composition = '(%s)+(%s)'%(fp[0], fp[1])
            self.cur_pdb().interface[-1].length_coord = '+'.join(interfaces[fp]['length_coord'])
            self.cur_pdb().interface[-1].clash = '+'.join(interfaces[fp]['clash'])
            self.cur_pdb().interface[-1].i_class = ''.join(set(interfaces[fp]['i_class']))
            contribution = ''
            for interface in interfaces[fp]['contribution']:
                contribution += '_'.join(interface[0])
                for res in interface[1:]:
                    contribution += '-'+'_'.join(res)
                contribution += '+'
            self.cur_pdb().interface[-1].contribution = contribution[:-1]
            # variables averaged
            for var in variables:
                final_value = str(round(sum(interfaces[fp][var])/len(interfaces[fp][var]), 2))
                setattr(self.cur_pdb().interface[-1], var, final_value)
            self.cur_pdb().interface[-1].areaimol_chain = str(round(float(self.cur_pdb().interface[-1].areaimol_chain), 1))
            self.cur_pdb().interface[-1].areaimol_c = str(round(float(self.cur_pdb().interface[-1].areaimol_c), 1))

        self.cur_pdb().interface.sort(key=lambda entry: -float(entry.areaimol_chain))

        print('')
        print('Number of unique interface(s) found: '+str(len(self.cur_pdb().interface)))

        # determine class
        if len(self.cur_pdb().interface) > 0:
            self.cur_pdb().clas = self.cur_pdb().interface[0].i_class
            for interface in self.cur_pdb().interface:
                if self.cur_pdb().clas != interface.i_class:
                    self.cur_pdb().sec_class = interface.i_class
                    print('Main class is %s with secondary class %s'%(self.cur_pdb().clas, self.cur_pdb().sec_class))
                    break
            else:
                print('Main class is %s'%self.cur_pdb().clas)
                self.cur_pdb().sec_class = '-'
        else:
            print('Unknown class')
            self.cur_pdb().clas = '-'
            self.cur_pdb().sec_class = '-'

        # determine packing
        chain_pairs = []
        for interface in self.cur_pdb().interface:
            chain_pairs += interface.sheets.replace('+', '_').split('_')
        interactions = 0
        for n in range(int(self.cur_pdb().n_Msheet)):
            interactions += chain_pairs.count(str(n+1))

        packing = interactions/float(self.cur_pdb().n_Msheet)
        if packing < 1.5:
            self.cur_pdb().packing = '1D'
            print('Packing is 1D')
        elif packing < 2.5:
            self.cur_pdb().packing = '2D'
            print('Packing is 2D')
        elif packing < 3.5:
            self.cur_pdb().packing = '3D'
            print('Packing determination is not certain, but it is likely 3D')
        else:
            self.cur_pdb().packing = '3D'
            print('Packing is 3D')

        # writes output
        with open('ACW_output_%s.txt'%self.main_object, 'w') as f:
            f.write('Evaluation of %s\n'%self.main_object)
            f.write('Main class is %s\n'%self.cur_pdb().clas)
            if self.cur_pdb().sec_class != '-':
                f.write('Secondary class is %s\n'%self.cur_pdb().sec_class)
            f.write('Packing is %s\n'%self.cur_pdb().packing)
            f.write('The number of interfaces is %s\n'%str(len(self.cur_pdb().interface)))
            f.write('\n')

            for n, interface in enumerate(self.cur_pdb().interface):
                f.write('Interface %s\n'%(n+1))
                sheets = ''
                for pairs in interface.sheets.split('+'):
                    sheets += '+'.join(pairs.split('_'))+', '
                f.write('Sheets: %s\n'%sheets[:-2])
                f.write('Sidechains: %s\n'%interface.composition)
                f.write('Sc: %s\n'%interface.sc)
                f.write('B.Area: %s\n'%interface.areaimol_chain)
                f.write('SDi: %s\n'%interface.perfi)
                f.write('Class: %s\n'%interface.i_class)

                clash = interface.clash.split('+')
                for x, c in enumerate(clash):
                    if c == 'yes':
                        pair = interface.sheets.split('+')[x].split('_')
                        f.write('Steric clash between %s and %s\n'%(pair[0], pair[1]))

                f.write('\n')
                f.write('Buried area details:\n')
                for i, c1 in enumerate(interface.contribution.split('+')):
                    f.write('Sheets: '+'+'.join((interface.sheets.split('+')[i].split('_')))+'\n')
                    f.write('residue\ttype\ts.ch.\tb.b.\ttotal\n')
                    for c2 in c1.split('-'):
                        c3 = c2.split('_')
                        f.write('%s\t%s\t%s\t%s\t%s\n'%(c3[0], c3[2], c3[3], c3[4], c3[5]))
                    f.write('\n')
                f.write('\n')

            if len(self.cur_pdb().interface) == 0:
                f.write('No interface found')

        # exit and display interfaces
        self.delete_all(['files'])
        if len(interfaces) > 0:
            self.interface_count = 0
            self.next_interface(verbose=False)
        else:
            cmd.color('green', 'all and elem C')
        cmd.set('suspend_updates', 0, quiet=1)
        if self.mode == 'evaluate':
            self.calc_state = [3]
        self.menu_update()


    def save_to_database(self):
        """Save finished evaluation to the personal database."""
        if self.calc_state[0] != 3:
            print('No calculation was performed - no data  to save to the database')
            self.error_message = 'No calculation was performed - no data  to save to the database'
            cmd.refresh_wizard()
            return

        print('\nSaving to database')

        # changing to database
        self.mode = 'database'
        self.pdb_count = len(self.pdb_list.lista)-1
        self.seq_len = int(self.cur_pdb().seq_len)
        cmd.cd(self.main_path)
        write_save('acw_personal_database.json', self.pdb_list)
        cmd.cd('%s'%self.cur_pdb().code)
        self.interface_count -= 1
        self.next_interface(verbose=False)
        self.menu_update()


    def delete_from_database(self, answer=0):
        """Delete a structure and its folder from the personal database."""
        if self.overwrite_state == 0:
            self.overwrite_state = 1
            self.menu_update()

        else:
            self.overwrite_state = 0
            self.menu_update()
            if self.cur_pdb().source == 'PDB':
                print('Only user data can be deleted')
                self.error_message = 'Only user data can be deleted'
                return

            if answer:
                print('%s was deleted'%self.cur_pdb().code)
                self.error_message = '%s was deleted'%self.cur_pdb().code
                x = self.pdb_count
                self.next_peptide(-1)
                shutil.rmtree('../%s'%self.pdb_list[x].code)
                self.pdb_list.lista.pop(x)
                self.next_peptide(0)


    def overwrite(self, answer=0):
        """Prompt the user if they want to overwrite a previous output folder."""

        self.overwrite_state = 0
        self.menu_update()
        if answer:
            for n, amy in enumerate(self.pdb_list.lista):
                if amy.code == self.main_object:
                    if amy.source != 'PDB':
                        self.pdb_list.lista.pop(n)
                        shutil.rmtree(self.main_object)
                        break
                    else:
                        print('This name is already used in the PDB database, cannot be overwritten')
                        self.error_message = 'This name is already used in the PDB database, cannot be overwritten'
                        return

            else:
                shutil.rmtree(self.main_object)
            print('Files were deleted')
            self.error_message = 'Files were deleted'
            cmd.get_wizard().launcher("multi", "setup")


    def exporter(self):
        """Export current selection to json and txt."""

        # create selection output
        output = amy_struct_list()
        for amy in self.pdb_list.lista:

            selection = True
            if self.sele_clas:
                if amy.clas not in self.sele_clas:
                    selection = False
            if self.sele_secclas:
                if amy.sec_class not in self.sele_secclas:
                    selection = False
            if self.sele_seq_len:
                if amy.seq_len not in self.sele_seq_len:
                    selection = False
            if self.sele_packing:
                if amy.packing not in self.sele_packing:
                    selection = False
            if self.sele_source:
                if amy.source not in self.sele_source:
                    selection = False
            if self.sele_aa:
                if not all(aa in amy.seq for aa in self.sele_aa):
                    selection = False

            if selection:
                output.add(amy)

        # write data
        if len(output.lista) > 0:

            name = ''
            if self.sele_clas:
                name += '(class:%s)'%"_".join(self.sele_clas)
            if self.sele_secclas:
                name += '(s.class:%s)'%"_".join(self.sele_secclas)
            if self.sele_seq_len:
                name += '(length:%s)'%"_".join(self.sele_seq_len)
            if self.sele_packing:
                name += '(packing:%s)'%"_".join(self.sele_packing)
            if self.sele_source:
                name += '(source:%s)'%"_".join(self.sele_source)
            if self.sele_aa:
                name += '(aa:%s)'%"_".join(self.sele_aa)

            cmd.cd(self.main_path)
            write_save('Selection_%s.json'%name, output)

            with open('Selection_%s_structures.txt'%name, 'w') as f:
                f.write('Code\tSeq.\tLength\tClass\tS.class\tPacking\n')
                for amy in output.lista:

                    f.write('\t'.join([amy.code, amy.seq, amy.seq_len, amy.clas, amy.sec_class, amy.packing])+'\n')

            with open('Selection_%s_interfaces.txt'%name, 'w') as f:
                f.write('Code\tSeq.\tClass\tSc\tB.area\tSDi\tComp.\n')
                for amy in output.lista:
                    for i_f in amy.interface:
                        f.write('\t'.join([amy.code, amy.seq, i_f.i_class, i_f.sc, i_f.areaimol_chain, i_f.perfi, i_f.composition])+'\n')

            cmd.cd('%s'%self.cur_pdb().code)

            print('The selection is saved')
            self.error_message = 'The selection is saved'

        else:
            print('The selection is empty')
            self.error_message = 'The selection is empty'


    def test(self):
        """Advanced test method 1."""


        pass



    def test2(self):
        """Advanced test method 2."""
        cmd.set('line_width', '2.5')
        cmd.set('nb_spheres_size', '0.350')
        cmd.set('stick_radius', '0.350')


    def test3(self):
        """Advanced test method 3."""


        pass



    def re_calc_method(self, only_this=False):
        """Advanced method for recalculating the database."""
        self.re_calc = True
        if not only_this:
            self.re_calc_all = True
        try:
            self.next_peptide(x=self.pdb_count)
        finally:
            cmd.set('suspend_updates', 0, quiet=1)
            self.re_calc = False
            self.re_calc_all = False


    def view(self):
        """Reorients the pymol viewer."""
        cmd.orient('all')

        my_view = list(cmd.get_view())
        views = {'x':[0, 0, 1, 1, 0, 0, 0, 1, 0],
                 'y':[0, 1, 0, 0, 0, 1, 1, 0, 0],
                 'z':[1, 0, 0, 0, 1, 0, 0, 0, 1],
                 0:[0, 0, 1, 1, 0, 0, 0, 1, 0]}

        for n in range(9):
            my_view[n] = views[self.main_axis][n]
        cmd.set_view(tuple(my_view))

        if self.view_count == 1:
            cmd.turn('z', 90)
        if self.view_count == 2:
            cmd.turn('x', 90)
        if self.view_count == 3:
            cmd.turn('y', 90)
            self.view_count = 0
            return

        self.view_count += 1


    def draw_axis_method(self):
        """Draws the coordinate axis."""
        if self.axis == 'ON':
            for n in ['X_axis', 'Y_axis', 'Z_axis', 'axis']:
                cmd.delete(n)
            self.axis = 'OFF'
        else:
            draw_axis()
            self.axis = 'ON'


    def manual_mode(self):
        """Turns on or of the manual mode."""
        if self.picking_mode == 'ON':
            self.pick_count = 0
            self.picking_mode = 'OFF'
            if len(self.cur_pdb().interface) > 0:
                self.next_interface()
        else:
            if self.mode == 'evaluate' and self.calc_state[0] not in [1, 3]:
                print('Not available without proper data')
                self.error_message = 'Not available without proper data'
                return

            if self.path:
                cmd.set('mouse_selection_mode', 4)
                cmd.util.cbag('all')
                cmd.delete('lengths')
                self.picking_mode = 'ON'
            else:
                print('CCP4 is not available, calculations cannot be performed')
                self.error_message = 'CCP4 is not available, calculations cannot be performed'
        cmd.refresh_wizard()


    def launcher(self, thread, name, *args, **kwargs):
        """Launch the other methods with single or multithreading."""
        try:
            self.error_message = ''
            self_fn = getattr(self, name)

            if thread == 'multi':
                t = threading.Thread(target=self_fn, args=args, kwargs=kwargs)
                t.setDaemon(1)
                t.start()

            elif thread == 'single':
                self_fn(*args, **kwargs)

            cmd.refresh_wizard()

        except Exception:
            cmd.set('suspend_updates', 0, quiet=1)
            print('Unknown error occurred')
            self.error_message = 'Unknown error occurred'


    def menu_update(self):
        """Updates the menu buttons."""
        for s in [self.sele_clas, self.sele_secclas, self.sele_seq_len, self.sele_packing]:
            s.sort(key=lambda x : str(len(x))+x)

        selection = [[2, 'List selection:', ''], [1, 'Class: ' + ','.join(self.sele_clas), []]]
        for t in self.sele_clas_types:
            selection[-1][-1].append([1, 'Class %s '%t, 'cmd.get_wizard().launcher("single","toggle_value","sele_clas",["%s"])'%t])
        selection[-1][-1].append([1, 'Parallel', 'cmd.get_wizard().launcher("single","toggle_value","sele_clas",["1","2","3","4"])'])
        selection[-1][-1].append([1, 'Antip.', 'cmd.get_wizard().launcher("single","toggle_value","sele_clas",["5","6","7","8"])'])

        selection.append([1, 'Sec. class: ' + ', '.join(self.sele_secclas), []])
        for t in self.sele_secclas_types:
            selection[-1][-1].append([1, 'S.c. %s'%t, 'cmd.get_wizard().launcher("single","toggle_value","sele_secclas",["%s"])'%t])
        selection[-1][-1].append([1, 'Parallel', 'cmd.get_wizard().launcher("single","toggle_value","sele_secclas",["1","2","3","4"])'])
        selection[-1][-1].append([1, 'Antip.', 'cmd.get_wizard().launcher("single","toggle_value","sele_secclas",["5","6","7","8"])'])

        selection.append([1, 'Length: ' + ','.join(self.sele_seq_len), []])
        for t in self.sele_seq_len_types:
            selection[-1][-1].append([1, ' %s '%t, 'cmd.get_wizard().launcher("single","toggle_value","sele_seq_len",["%s"])'%t])

        selection.append([1, 'Packing: ' + ','.join(self.sele_packing), []])
        for t in self.sele_packing_types:
            selection[-1][-1].append([1, ' %s '%t, 'cmd.get_wizard().launcher("single","toggle_value","sele_packing",["%s"])'%t])

        selection.append([1, 'Source: ' + ','.join(self.sele_source), []])
        for t in self.sele_source_types:
            selection[-1][-1].append([1, ' %s '%t, 'cmd.get_wizard().launcher("single","toggle_value","sele_source",["%s"])'%t])

        selection.append([1, 'Amino acid: ' + ','.join(self.sele_aa), []])
        for t in self.sele_aa_types:
            selection[-1][-1].append([1, ' %s '%t, 'cmd.get_wizard().launcher("single","toggle_value","sele_aa",["%s"])'%t])

        selection.append([1, 'Reset criteria ', 'cmd.get_wizard().launcher("single","reset_selection")'])
        self.menu['selection'] = selection

        self.menu['Advanced'] = [[2, 'Advanced', ''],
                                 [1, 'Prompt: %s' %self.need_prompt, 'cmd.get_wizard().def_changer("need_prompt")'],
                                 [1, 'test1', 'cmd.get_wizard().launcher("multi","test")'],
                                 [1, 'test2', 'cmd.get_wizard().launcher("multi","test2")'],
                                 [1, 'test3', 'cmd.get_wizard().launcher("multi","test3")'],
                                 [1, 'Toggle axes', 'cmd.get_wizard().draw_axis_method()'],
                                 [1, 'Extend the last sheet', 'cmd.get_wizard().extend_sheet()'],
                                 [1, 'Picking: Delete', 'cmd.get_wizard().set_value("picking_mode","Delete","wizard")'],
                                 [1, 'Picking: Color M', 'cmd.get_wizard().set_value("picking_mode","Color M","wizard")'],
                                 [1, 'Picking: Color Y', 'cmd.get_wizard().set_value("picking_mode","Color Y","wizard")'],
                                 [1, 'Recalculation of database', 'cmd.get_wizard().launcher("multi","re_calc_method")'],
                                 [1, 'Recalculation of this', 'cmd.get_wizard().launcher("multi","re_calc_method",True)']
                                 ]

        cmd.refresh_wizard()


    def get_panel(self):
        """Configures the buttons on right panel."""

        if self.pdb_count == -1 and self.mode == 'database':
            view_interface = 'View interfaces (0/0)'
        elif len(self.cur_pdb().interface) == 0:
            view_interface = 'View interfaces (0/0)'
        else:
            view_interface = 'View interfaces (%s/%s)'%(self.interface_count, len(self.cur_pdb().interface))

        if self.picking_mode == 'OFF':
            pick_mode = 'ON'
        else:
            pick_mode = 'OFF'

        if self.mode == 'database':
            panel = [
                [1, 'ACW Wizard: Database', ''],
                [3, 'Selection criteria', 'selection'],
                [2, 'Next peptide', 'cmd.get_wizard().launcher("single","next_peptide")'],
                [2, 'Previous peptide', 'cmd.get_wizard().launcher("single","next_peptide",-1)'],

                [1, 'Viewing:', ''],
                [2, view_interface, 'cmd.get_wizard().launcher("single","next_interface")'],
                [2, 'Export image', 'cmd.get_wizard().launcher("single","exp_image")'],
                [2, 'Toggle labels', 'cmd.get_wizard().launcher("single","toggle_labels")'],
                [2, 'Reset/rotate view', 'cmd.get_wizard().launcher("single","view")'],

                [1, 'Edit:', ''],
                [2, 'Manual Mode: %s'%self.picking_mode, 'cmd.get_wizard().launcher("single","manual_mode")'],
                [2, 'Export current selection', 'cmd.get_wizard().launcher("single","exporter")'],

                [1, 'Menu:', ''],
                [2, 'Change to ACW Evaluation', 'cmd.get_wizard().launcher("single","mode_change","evaluate")'],
                [2, 'Quit', 'cmd.set_wizard()']
                ]

            if self.cur_pdb().source == 'personal':
                panel.insert(11, [2, 'Delete from database', 'cmd.get_wizard().launcher("single","delete_from_database")'])

            if self.overwrite_state:
                panel.insert(0, [1, 'This deletes its folder,', 1])
                panel.insert(1, [1, 'are you sure?', 1])
                panel.insert(2, [2, 'Yes', 'cmd.get_wizard().launcher("single","delete_from_database",1)'])
                panel.insert(3, [2, 'No', 'cmd.get_wizard().launcher("single","delete_from_database",0)'])

        elif self.mode == 'evaluate':

            panel = [
                [1, 'ACW Wizard: Evaluation', 1],
                [2, '1: Set up the sheets', 'cmd.get_wizard().launcher("multi","setup")'],
                [2, '2: Automatic calculation', 'cmd.get_wizard().launcher("multi","automatic_calc")'],
                [2, '3: '+view_interface, 'cmd.get_wizard().launcher("single","next_interface")'],
                [2, '4: Save to database', 'cmd.get_wizard().launcher("single","save_to_database")'],
                [2, '5: Reset', 'cmd.get_wizard().launcher("single","delete_all",["everything","interface","last_amy"])'],

                [1, 'Viewing:', ''],
                [2, 'Export image', 'cmd.get_wizard().launcher("single","exp_image")'],
                [2, 'Toggle labels', 'cmd.get_wizard().launcher("single","toggle_labels")'],
                [2, 'Reset/rotate view', 'cmd.get_wizard().launcher("single","view")'],

                [1, 'Edit:', ''],
                [2, 'Manual Mode: %s'%self.picking_mode, 'cmd.get_wizard().launcher("single","manual_mode")'],

                [1, 'Menu:', ''],
                [2, 'Change to ACW Database', 'cmd.get_wizard().launcher("single","mode_change","database")'],
                [2, 'Quit', 'cmd.set_wizard()']
                ]

            if self.overwrite_state:
                panel.insert(0, [1, 'Overwrite existing data?', 1])
                panel.insert(1, [2, 'Yes', 'cmd.get_wizard().launcher("single","overwrite",1)'])
                panel.insert(2, [2, 'No', 'cmd.get_wizard().launcher("single","overwrite",0)'])

        if self.advanced:
            panel.append([3, 'Advanced', 'Advanced'])

        return panel


def acw_starter():
    try:
        cmd.set_wizard(AcWizardMultiSheet(acw_database_path))
    except Exception:
        cmd.set_wizard()

cmd.extend('acw', acw_starter)


def acw_goto(goto=0):
    try:
        cmd.get_wizard().cmd.get_wizard().launcher('single', 'next_peptide', goto=goto)
    except AttributeError:
        cmd.set_wizard(AcWizardMultiSheet(acw_database_path))
        cmd.get_wizard().cmd.get_wizard().launcher('single', 'next_peptide', goto=goto)

cmd.extend('acwgoto', acw_goto)

print("You can now type 'acw' to start the Amyloid Coordinate Wizard")
