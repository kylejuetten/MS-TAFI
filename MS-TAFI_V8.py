# -*- coding: utf-8 -*-

# Import Libraries

from pandastable import Table, config
import numpy as np
import math
import plotly.graph_objects as go
from plotly.offline import plot
import plotly.express as px
import tkinter as tk
from tkinter.filedialog import askopenfilename, askopenfile
import pandas as pd
#pd.options.mode.chained_assignment = None
import re
#

LARGEFONT = ("Verdana", 35)


class tkinterApp(tk.Tk):

    # __init__ function for class tkinterApp
    def __init__(self, *args, **kwargs):

        # __init__ function for class Tk
        tk.Tk.__init__(self, *args, **kwargs)
        self.title('MS-TAFI')

        # creating a container
        container = tk.Frame(self)
        container.pack(side="top", fill="both", expand=True)
        global ws
        global hs
        ws = self.winfo_screenwidth()
        hs = self.winfo_screenheight()

        # define window size
        w = ws/2
        h = hs*0.935
        x = (ws) - (ws/0.995)
        y = (hs) - (hs)
        self.geometry('%dx%d+%d+%d' % (w, h, x, y))

        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)

        # initializing frames to an empty array
        self.frames = {}

        # iterating through a tuple consisting
        # of the different page layouts
        for F in (StartPage, ResultsPage):

            frame = F(container, self)

            # initializing frame of that object from
            # startpage, page1, page2 respectively with
            # for loop
            self.frames[F] = frame

            frame.grid(row=0, column=0, sticky="nsew")

        self.show_frame(StartPage)

    # to display the current frame passed as
    # parameter
    def show_frame(self, cont):
        frame = self.frames[cont]
        frame.tkraise()

# first window frame startpage


class StartPage(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.controller = controller
        # Initialize Window
        self.columnconfigure(tuple(range(6)), weight=1)
        self.rowconfigure(tuple(range(27)), weight=1)

        # App Title
        title_label = tk.Label(
            self, text='Mass Spectrometry Tool for the Analysis of Fragment Ions')
        title_label.config(font=('helvetica', 14, 'bold'))
        title_label.grid(column=0, row=0, columnspan=6)

        # Open Deconvoluted Mass + Intensity List
        csv_label = tk.Label(
            self, text='Upload Deconvoluted Mass + Intensity List:')
        csv_label.config(font=('helvetica', 10, 'bold'))
        csv_label.grid(column=0, row=1, columnspan=6)

        def clicked():
            global filename
            filename = askopenfilename(filetypes=(("CSV Files", "*.csv"),))
            file_label['text'] = filename
        open_csv_button = tk.Button(
            self, text="Upload", command=clicked, bg='brown', fg='white', font=('helvetica', 9, 'bold'))
        open_csv_button.grid(column=0, row=2, columnspan=6)
        file_label = tk.Label(self, text='')
        file_label.grid(column=0, row=3, columnspan=6)

        # Protein Sequence Entry
        seq_label = tk.Label(self, text='Protein Sequence: ')
        seq_label.config(font=('helvetica', 10, 'bold'))
        seq_label.grid(column=0, row=4, columnspan=6)
        seq_entry = tk.Text(self, height=10, wrap=tk.WORD)
        seq_entry.configure(font=('helvetica', 10))
        seq_entry.grid(column=0, row=5, columnspan=6)

        # PTM position Entry
        PTM_position_var = tk.StringVar()
        PTM_position_label = tk.Label(
            self, text='PTM position(s): \n (comma separated)')
        PTM_position_label.config(font=('helvetica', 10, 'bold'))
        PTM_position_label.grid(column=0, row=6, columnspan=2)
        PTM_position_entry = tk.Entry(
            self, textvariable=PTM_position_var, justify='center', width=30)
        PTM_position_entry.grid(column=0, row=7, columnspan=2)

        # PTM mass Entry
        PTM_mass_var = tk.StringVar()
        PTM_mass_label = tk.Label(
            self, text='PTM mass(es): \n (comma separated)')
        PTM_mass_label.config(font=('helvetica', 10, 'bold'))
        PTM_mass_label.grid(column=2, row=6, columnspan=2)
        PTM_mass_entry = tk.Entry(
            self, textvariable=PTM_mass_var, justify='center', width=30)
        PTM_mass_entry.grid(column=2, row=7, columnspan=2)

        # Ligand mass Entry
        ligand_mass_var = tk.StringVar()
        ligand_mass_label = tk.Label(self, text='Ligand mass: \n')
        ligand_mass_label.config(font=('helvetica', 10, 'bold'))
        ligand_mass_label.grid(column=4, row=6, columnspan=2)
        ligand_mass_entry = tk.Entry(
            self, textvariable=ligand_mass_var, justify='center', width=30)
        ligand_mass_entry.grid(column=4, row=7, columnspan=2)

        # Activation Entry
        menu = tk.StringVar()
        menu.set('Activation method')
        drop = tk.OptionMenu(self, menu, 'CID', 'HCD', 'ETD', 'EThcD', 'UVPD')
        drop.config(bg='brown', fg='white', activebackground='brown',
                    activeforeground='white', font=('helvetica', 9, 'bold'))
        drop['menu'].config(bg='brown', fg='white', activebackground='brown',
                            activeforeground='white', font=('helvetica', 9, 'bold'))
        activation_label = tk.Label(self, text='Activation:')
        activation_label.config(font=('helvetica', 10, 'bold'))
        activation_label.grid(column=0, row=8, columnspan=2)
        drop.grid(column=0, row=9, columnspan=2)

        # PPM Error Entry
        PPM_var = tk.StringVar(value='10')
        PPM_label = tk.Label(self, text='ppm error:')
        PPM_label.config(font=('helvetica', 10, 'bold'))
        PPM_label.grid(column=2, row=8, columnspan=2)
        PPM_entry = tk.Entry(self, textvariable=PPM_var, justify='center')
        PPM_entry.grid(column=2, row=9, columnspan=2)

        # TIC Entry
        tic_var = tk.StringVar(value='1')
        tic_label = tk.Label(self, text='TIC:')
        tic_label.config(font=('helvetica', 10, 'bold'))
        tic_label.grid(column=4, row=8, columnspan=2)
        tic_entry = tk.Entry(self, textvariable=tic_var, justify='center')
        tic_entry.grid(column=4, row=9, columnspan=2)

        # Holo Ion mode
        holo_menu = tk.StringVar()
        holo_menu.set('Apo')
        drop = tk.OptionMenu(self, holo_menu, 'Apo', 'Holo')
        drop.config(bg='brown', fg='white', activebackground='brown',
                    activeforeground='white', font=('helvetica', 9, 'bold'))
        drop['menu'].config(bg='brown', fg='white', activebackground='brown',
                            activeforeground='white', font=('helvetica', 9, 'bold'))
        holo_label = tk.Label(self, text='Apo or Holo:')
        holo_label.config(font=('helvetica', 10, 'bold'))
        holo_label.grid(column=3, row=10, columnspan=2)
        drop.grid(column=3, row=11, columnspan=2)
        
        # Neutral loss search
        neutral_loss_menu = tk.StringVar()
        neutral_loss_menu.set('No')
        neutral_loss_drop = tk.OptionMenu(self, neutral_loss_menu, 'Yes', 'No')
        neutral_loss_drop.config(bg='brown', fg='white', activebackground='brown',
                    activeforeground='white', font=('helvetica', 9, 'bold'))
        neutral_loss_drop['menu'].config(bg='brown', fg='white', activebackground='brown',
                            activeforeground='white', font=('helvetica', 9, 'bold'))
        neutral_loss_label = tk.Label(self, text='Search for Netural Losses:')
        neutral_loss_label.config(font=('helvetica', 10, 'bold'))
        neutral_loss_label.grid(column=1, row=10, columnspan=2)
        neutral_loss_drop.grid(column=1, row=11, columnspan=2)

        
        # Results Label
        results_label = tk.Label(self, text='Results \n \n Protein Mass: ' + '         ' +
                                 '\t Sequence Length: ' + '   ' + '    ' +
                                 '\t Sequence Coverage: ' + '     ' + ' ' +
                                 '\t Fragments Explained: ' + '     ' + '  \n')
        results_label.config(font=('helvetica', 12))
        results_label.grid(column=0, row=13, columnspan=6)
        self.rowconfigure(13, weight=3)
        
        
        

        # Functions

        def user_ions(activation):

            # Define ion types according to activation method
            HCD = ['b', 'b+H2O', 'b-H2O', 'b-NH3', 'y', 'y-H2O', 'y-NH3']
            ETD = ['c', 'z']
            EThcD = ['b', 'b+H2O', 'b-H2O', 'b-NH3','c', 'y', 'y-H2O', 'y-NH3', 'z']
            UVPD = ['a', 'a-NH3', 'a+1', 'd', 'b', 'b+H2O', 'b-H2O', 'b-NH3', 'c', 'x', 'x+1', 'v', 'w', 'y', 'y-1', 'y-2', 'y-H2O', 'y-NH3', 'z']

            ion_types = []

            # Assign ion types for search
            if activation == 'HCD' or activation == 'CID':
                ion_types = HCD

            if activation == 'ETD':
                ion_types = ETD

            if activation == 'ETHCD':
                ion_types = EThcD

            if activation == 'UVPD':
                ion_types = UVPD

            return ion_types

        def user_TIC(tic2):

            # Control if no value is entered for TIC
            if tic2 == '':
                tic2 = '1'

            if tic2 == '0':
                tic2 = '1'

            return int(tic2)

        def user_PTMs(x, y):

            # Control if no values are entered for PTM position/mass

            if x == '':
                x = [0]

            else:
                x = [int(x) for x in x.split(',')]

            if y == '':
                y = [0]

            else:
                y = [float(x) for x in y.split(',')]

            PTM_dict = {x[i]: y[i] for i in range(len(x))}

            return (PTM_dict)

        def user_Ligand(x):

            if x == '':
                x = [0]

            else:
                x = [float(x) for x in x.split(',')]

            return (x)

        def theoretical_mass(seq, ptm_dict, Ligand_mass):

            # Define amino acid monoisotopic masses
            AA_masses = {"A": 71.03711, "C": 103.00919, "D": 115.02694, "E": 129.04259,
                         "F": 147.06841, "G": 57.02146, "H": 137.05891, "I": 113.08406,
                         "K": 128.09496, "L": 113.08406, "M": 131.04049, "N": 114.04293,
                         "P": 97.05276, "Q": 128.05858, "R": 156.10111, "S": 87.03203,
                         "T": 101.04768, "V": 99.06841, "W": 186.07931, "Y": 163.06333}

            # Split protein sequence into list
            protein_aa_list = [aa for aa in protein_sequence]

            protein_monoisotopic_mass = 0

            # Sum protein monoisotopic mass
            for i in protein_aa_list:
                if i in AA_masses:
                    protein_monoisotopic_mass += AA_masses[i]

            # Correct mass according to PTM's
            for value in ptm_dict.values():
                protein_monoisotopic_mass += value

            # Correct mass according to Ligands
            for i in Ligand_mass:
                protein_monoisotopic_mass += i

            # Correct mass for N/C term
            protein_monoisotopic_mass += 18.01

            return round(protein_monoisotopic_mass, 2)

        def theoretical_fragment_calc(seq, ptm_dict, Ligand_mass):

            # Define amino acid monoisotopic masses
            AA_masses = {"A": 71.03711, "C": 103.00919, "D": 115.02694, "E": 129.04259,
                         "F": 147.06841, "G": 57.02146, "H": 137.05891, "I": 113.08406,
                         "K": 128.09496, "L": 113.08406, "M": 131.04049, "N": 114.04293,
                         "P": 97.05276, "Q": 128.05858, "R": 156.10111, "S": 87.03203,
                         "T": 101.04768, "V": 99.06841, "W": 186.07931, "Y": 163.06333}

            # Split sequence into list
            protein_aa_list = [aa for aa in protein_sequence]

            mass_list = []

            for i in protein_aa_list:
                if i in AA_masses:
                    mass_list.append(AA_masses[i])

            # Add PTM's to corresponding position on sequence
            for index, value in enumerate(mass_list):
                for k, v in ptm_dict.items():
                    if index == (k-1):
                        mass_list[index] += v

            # Dictionary of lists for theoretical fragments
            fragment_dict = {'a': [], 'a+1': [], 'b': [], 'c': [],
                'x': [], 'x+1': [], 'y': [], 'y-1': [], 'y-2': [], 'z': []}

            # N-Term fragment list
            fragment_mass = 0
            N_term_fragment_list = []

            for index, value in enumerate(mass_list):
                N_term_fragment_list.append(fragment_mass)
                fragment_mass += mass_list[index]

            # C-Term fragment list
            fragment_mass = 0
            C_term_fragment_list = []

            for index, value in reversed(list(enumerate(mass_list))):
                C_term_fragment_list.append(fragment_mass)
                fragment_mass += mass_list[index]

            # Correct N-Term fragments according to fragment type
            for i in N_term_fragment_list:
                fragment_dict['a'].append(i-27.9949)
                fragment_dict['a+1'].append(i-26.987075)
                fragment_dict['b'].append(i)
                fragment_dict['c'].append(i+17.026549)

            # Correct C-Term fragments according to fragment type
            for i in C_term_fragment_list:
                fragment_dict['x'].append(i+43.989830)
                fragment_dict['x+1'].append(i+44.997655)
                fragment_dict['y'].append(i+18.010565)
                fragment_dict['y-1'].append(i+17.00274)
                fragment_dict['y-2'].append(i+15.994915)
                fragment_dict['z'].append(i+1.991841)

            # Dictionary of lists for theoretical fragments with ligand
            ligand_fragment_dict = {'a': [], 'a+1': [], 'b': [], 'c': [],
                'x': [], 'x+1': [], 'y': [], 'y-1': [], 'y-2': [], 'z': []}

            # N-Term fragment list
            fragment_mass = 0
            ligand_N_term_fragment_list = []

            for index, value in enumerate(mass_list):
                ligand_N_term_fragment_list.append(fragment_mass)
                fragment_mass += mass_list[index]

            # C-Term fragment list
            fragment_mass = 0
            ligand_C_term_fragment_list = []

            for index, value in reversed(list(enumerate(mass_list))):
                ligand_C_term_fragment_list.append(fragment_mass)
                fragment_mass += mass_list[index]

            # Correct for Ligand mass
            for i in ligand_N_term_fragment_list:
                for j in Ligand_mass:
                    ligand_fragment_dict['a'].append(i-27.9949 + j)
                    ligand_fragment_dict['a+1'].append(i-26.987075 + j)
                    ligand_fragment_dict['b'].append(i + j)
                    ligand_fragment_dict['c'].append(i+17.026549 + j)

            for i in ligand_C_term_fragment_list:
                for j in Ligand_mass:
                    ligand_fragment_dict['x'].append(i+43.989830 + j)
                    ligand_fragment_dict['x+1'].append(i+44.997655 + j)
                    ligand_fragment_dict['y'].append(i+18.010565 + j)
                    ligand_fragment_dict['y-1'].append(i+17.00274 + j)
                    ligand_fragment_dict['y-2'].append(i+15.994915 + j)
                    ligand_fragment_dict['z'].append(i+1.991841 + j)
                    
            # Dictionary of lists for theoretical fragments with neutral losses
            global neutral_loss_fragment_dict 
            neutral_loss_fragment_dict = {'a': [], 'a-NH3': [],
                                          'a+1': [],
                                          'd':[],
                                          'b': [], 'b-NH3': [], 'b-H2O': [], 'b+H2O': [], 
                                          'c': [],
                                          'x': [], 
                                          'x+1': [],
                                          'v':[], 'w':[],
                                          'y': [], 'y-NH3': [], 'y-H2O': [], 
                                          'y-1': [], 
                                          'y-2': [], 
                                          'z': []}
            
            # Dictionary for side chain losses of w ions
            w_side_chain_loss_dict = {
                'K': 57.05785,
                'D': 43.98983,
                'E': 58.00548,
                'R': 85.06399,
                'H': -43.00581, #Correct for 0 loss
                'S': 15.99492,
                'T': 15.99492,
                'N': 43.00581,
                'Q': 57.02146,
                'C': 31.97207,
                'G': -43.00581,  #Correct for 0 loss
                'P': -43.00581, #Correct for 0 loss
                'A': -43.00581, #Correct for 0 loss
                'I': 28.03130,
                'L': 42.04695,
                'M': 60.00337,
                'F': -43.00581, #Correct for 0 loss
                'W': -43.00581, #Correct for 0 loss
                'Y': -43.00581, #Correct for 0 loss
                'V': 14.01565}
                
                # Dictionary for side chain losses of d ions
            d_side_chain_loss_dict = {
                    'K': 57.05785,
                    'D': 43.98983,
                    'E': 58.00548,
                    'R': 85.06399,
                    'H': 0, #Correct for 0 loss
                    'S': 15.99492,
                    'T': 15.99492,
                    'N': 43.00581,
                    'Q': 57.02146,
                    'C': 31.97207,
                    'G': 0,  #Correct for 0 loss
                    'P': 0, #Correct for 0 loss
                    'A': 0, #Correct for 0 loss
                    'I': 28.03130,
                    'L': 42.04695,
                    'M': 60.00337,
                    'F': 0, #Correct for 0 loss
                    'W': 0, #Correct for 0 loss
                    'Y': 0, #Correct for 0 loss
                    'V': 14.01565}
            
            
            # N-Term fragment list
            fragment_mass = 0
            neutral_loss_N_term_fragment_list = []

            for index, value in enumerate(mass_list):
                neutral_loss_N_term_fragment_list.append(fragment_mass)
                fragment_mass += mass_list[index]

            # C-Term fragment list
            fragment_mass = 0
            neutral_loss_C_term_fragment_list = []

            for index, value in reversed(list(enumerate(mass_list))):
                neutral_loss_C_term_fragment_list.append(fragment_mass)
                fragment_mass += mass_list[index]
            
            # Correct for Ligand mass
            for i in neutral_loss_N_term_fragment_list:
                
                    amino_acid = protein_aa_list[neutral_loss_N_term_fragment_list.index(i)-1]
                    d_side_chain = d_side_chain_loss_dict[amino_acid]
                    
                                    
                    neutral_loss_fragment_dict['a'].append(i-27.9949)
                    neutral_loss_fragment_dict['a-NH3'].append(i-27.9949-17.0265)
                    neutral_loss_fragment_dict['a+1'].append(i-26.987075)
                    neutral_loss_fragment_dict['d'].append(i-27.9949- d_side_chain) #a - partial side chain
                    neutral_loss_fragment_dict['b'].append(i)
                    neutral_loss_fragment_dict['b-NH3'].append(i-17.0265)
                    neutral_loss_fragment_dict['b-H2O'].append(i-18.0106)
                    neutral_loss_fragment_dict['b+H2O'].append(i+18.0106)
                    neutral_loss_fragment_dict['c'].append(i+17.026549)

            for i in neutral_loss_C_term_fragment_list:
                
                    amino_acid = protein_aa_list[-(neutral_loss_C_term_fragment_list.index(i))]
                    amino_acid_mass = AA_masses[amino_acid]
                    side_chain = amino_acid_mass-56.04747
                    w_side_chain = w_side_chain_loss_dict[amino_acid]
                    
                    neutral_loss_fragment_dict['x'].append(i+43.989830)
                    neutral_loss_fragment_dict['x+1'].append(i+44.997655)
                    neutral_loss_fragment_dict['v'].append(i+43.989830-side_chain-27.9949) # full side chain - CO
                    neutral_loss_fragment_dict['w'].append(i+43.989830-w_side_chain-43.00581) #x - partial side chain - CHNO
                    neutral_loss_fragment_dict['y'].append(i+18.010565)
                    neutral_loss_fragment_dict['y-NH3'].append(i+18.010565-17.0265)
                    neutral_loss_fragment_dict['y-H2O'].append(i+18.0106-18.0106)
                    neutral_loss_fragment_dict['y-1'].append(i+17.00274)
                    neutral_loss_fragment_dict['y-2'].append(i+15.994915)
                    neutral_loss_fragment_dict['z'].append(i+1.991841)


            return (fragment_dict, ligand_fragment_dict, neutral_loss_fragment_dict)

        def chemical_formula(df):
            # Define a dictionary containing the atomic composition of each amino acid
            aa_compositions = {
                'A': {'C': 3, 'H': 5, 'N': 1, 'O': 1, 'S': 0},
                'R': {'C': 6, 'H': 12, 'N': 4, 'O': 1, 'S': 0},
                'N': {'C': 4, 'H': 6, 'N': 2, 'O': 2, 'S': 0},
                'D': {'C': 4, 'H': 5, 'N': 1, 'O': 3, 'S': 0},
                'C': {'C': 3, 'H': 5, 'N': 1, 'O': 1, 'S': 1},
                'E': {'C': 5, 'H': 7, 'N': 1, 'O': 3, 'S': 0},
                'Q': {'C': 5, 'H': 8, 'N': 2, 'O': 2, 'S': 0},
                'G': {'C': 2, 'H': 3, 'N': 1, 'O': 1, 'S': 0},
                'H': {'C': 6, 'H': 7, 'N': 3, 'O': 1, 'S': 0},
                'I': {'C': 6, 'H': 11, 'N': 1, 'O': 1, 'S': 0},
                'L': {'C': 6, 'H': 11, 'N': 1, 'O': 1, 'S': 0},
                'K': {'C': 6, 'H': 12, 'N': 2, 'O': 1, 'S': 0},
                'M': {'C': 5, 'H': 9, 'N': 1, 'O': 1, 'S': 1},
                'F': {'C': 9, 'H': 9, 'N': 1, 'O': 1, 'S': 0},
                'P': {'C': 5, 'H': 7, 'N': 1, 'O': 1, 'S': 0},
                'S': {'C': 3, 'H': 5, 'N': 1, 'O': 2, 'S': 0},
                'T': {'C': 4, 'H': 7, 'N': 1, 'O': 2, 'S': 0},
                'W': {'C': 11, 'H': 10, 'N': 2, 'O': 1, 'S': 0},
                'Y': {'C': 9, 'H': 9, 'N': 1, 'O': 2, 'S': 0},
                'V': {'C': 5, 'H': 9, 'N': 1, 'O': 1, 'S': 0},
            }
            #Full side chain compositions for v ions
            aa_side_chain_compositions = {
                'A': {'C': 1, 'H': 3, 'N': 0, 'O': 0, 'S': 0},
                'R': {'C': 4, 'H': 10, 'N': 3, 'O': 0, 'S': 0},
                'N': {'C': 2, 'H': 4, 'N': 1, 'O': 1, 'S': 0},
                'D': {'C': 2, 'H': 3, 'N': 0, 'O': 2, 'S': 0},
                'C': {'C': 1, 'H': 3, 'N': 0, 'O': 0, 'S': 1},
                'E': {'C': 3, 'H': 5, 'N': 0, 'O': 2, 'S': 0},
                'Q': {'C': 3, 'H': 6, 'N': 1, 'O': 1, 'S': 0},
                'G': {'C': 0, 'H': 1, 'N': 0, 'O': 0, 'S': 0},
                'H': {'C': 4, 'H': 5, 'N': 2, 'O': 0, 'S': 0},
                'I': {'C': 4, 'H': 9, 'N': 0, 'O': 0, 'S': 0},
                'L': {'C': 4, 'H': 9, 'N': 0, 'O': 0, 'S': 0},
                'K': {'C': 4, 'H': 10, 'N': 1, 'O': 0, 'S': 0},
                'M': {'C': 3, 'H': 7, 'N': 0, 'O': 0, 'S': 1},
                'F': {'C': 7, 'H': 7, 'N': 0, 'O': 0, 'S': 0},
                'P': {'C': 3, 'H': 5, 'N': 0, 'O': 0, 'S': 0},
                'S': {'C': 1, 'H': 3, 'N': 0, 'O': 1, 'S': 0},
                'T': {'C': 2, 'H': 5, 'N': 0, 'O': 1, 'S': 0},
                'W': {'C': 9, 'H': 8, 'N': 1, 'O': 0, 'S': 0},
                'Y': {'C': 7, 'H': 7, 'N': 0, 'O': 1, 'S': 0},
                'V': {'C': 3, 'H': 7, 'N': 0, 'O': 0, 'S': 0},
            }
            #Partial side chain losses for d and w ions
            d_w_aa_side_chain_compositions = {
                'A': {'C': 0, 'H': 0, 'N': 0, 'O': 0, 'S': 0},
                'R': {'C': 3, 'H': 7, 'N': 3, 'O': 0, 'S': 0},
                'N': {'C': 1, 'H': 1, 'N': 1, 'O': 1, 'S': 0},
                'D': {'C': 1, 'H': 0, 'N': 0, 'O': 2, 'S': 0},
                'C': {'C': 0, 'H': 0, 'N': 0, 'O': 0, 'S': 1},
                'E': {'C': 2, 'H': 2, 'N': 0, 'O': 2, 'S': 0},
                'Q': {'C': 2, 'H': 3, 'N': 1, 'O': 1, 'S': 0},
                'G': {'C': 0, 'H': 0, 'N': 0, 'O': 0, 'S': 0},
                'H': {'C': 0, 'H': 0, 'N': 0, 'O': 0, 'S': 0},
                'I': {'C': 2, 'H': 4, 'N': 0, 'O': 0, 'S': 0},
                'L': {'C': 3, 'H': 6, 'N': 0, 'O': 0, 'S': 0},
                'K': {'C': 3, 'H': 7, 'N': 1, 'O': 0, 'S': 0},
                'M': {'C': 2, 'H': 4, 'N': 0, 'O': 0, 'S': 1},
                'F': {'C': 0, 'H': 0, 'N': 0, 'O': 0, 'S': 0},
                'P': {'C': 0, 'H': 0, 'N': 0, 'O': 0, 'S': 0},
                'S': {'C': 0, 'H': 0, 'N': 0, 'O': 1, 'S': 0},
                'T': {'C': 0, 'H': 0, 'N': 0, 'O': 1, 'S': 0},
                'W': {'C': 0, 'H': 0, 'N': 0, 'O': 0, 'S': 0},
                'Y': {'C': 0, 'H': 0, 'N': 0, 'O': 0, 'S': 0},
                'V': {'C': 1, 'H': 2, 'N': 0, 'O': 0, 'S': 0},
                }
            
            C, H, O, N, S = [], [], [], [], []
            
            for i in df.itertuples():
                ion_type = i._1
                
                if ion_type == 'a' or ion_type == 'a+1' or ion_type == 'd' or ion_type == 'b' or ion_type == 'c' or \
                    ion_type == 'a-NH3' or ion_type == 'b-H2O' or ion_type == 'b-H2O' or ion_type == 'b+H2O' or ion_type == 'b-NH3':
                        
                    seq = protein_sequence[:(i.position)]
                    
                else:
                    
                    seq = protein_sequence[(-i.position):]
            
                
                peptide_composition = {}
                for aa in seq:
                    aa_composition = aa_compositions[aa]
                    for element, count in aa_composition.items():
                        peptide_composition[element] = peptide_composition.get(element, 0) + count
                
                
                if ion_type == 'a':
                    peptide_composition['C'] += -1
                    peptide_composition['O'] += -1
                    peptide_composition['H'] += 0
                    peptide_composition['S'] += 0
                        
                if ion_type == 'a+1':
                    peptide_composition['C'] += -1
                    peptide_composition['O'] += -1
                    peptide_composition['H'] += 1
                    peptide_composition['S'] += 0
                    
                if ion_type == 'd':  #a ion - partial side chain loss
                    peptide_composition['C'] += -1 - d_w_aa_side_chain_compositions[aa]['C']
                    peptide_composition['O'] += -1 - d_w_aa_side_chain_compositions[aa]['O']
                    peptide_composition['H'] += 0 - d_w_aa_side_chain_compositions[aa]['H']
                    peptide_composition['S'] += 0 - d_w_aa_side_chain_compositions[aa]['S']
                    peptide_composition['N'] += - d_w_aa_side_chain_compositions[aa]['N']
                        
                if ion_type == 'a-NH3':
                    peptide_composition['C'] += -1
                    peptide_composition['O'] += -1
                    peptide_composition['H'] += -3
                    peptide_composition['N'] += -1
                    peptide_composition['S'] += 0
                        
                if ion_type == 'b':
                    peptide_composition['C'] += 0
                    peptide_composition['O'] += 0
                    peptide_composition['H'] += 0
                    peptide_composition['S'] += 0
                
                if ion_type == 'b-H2O':
                    peptide_composition['C'] += 0
                    peptide_composition['O'] += -1
                    peptide_composition['H'] += -2
                    peptide_composition['S'] += 0
                
                if ion_type == 'b-NH3':
                    peptide_composition['C'] += 0
                    peptide_composition['O'] += 0
                    peptide_composition['H'] += -3
                    peptide_composition['N'] += -1
                    peptide_composition['S'] += 0
            
                if ion_type == 'b+H2O':
                    peptide_composition['C'] += 0
                    peptide_composition['O'] += 1
                    peptide_composition['H'] += 2
                    peptide_composition['S'] += 0
                
                if ion_type == 'c':
                    peptide_composition['C'] += 0
                    peptide_composition['O'] += 0
                    peptide_composition['H'] += 3
                    peptide_composition['N'] += 1
                    peptide_composition['S'] += 0
                
                if ion_type == 'x':
                    peptide_composition['C'] += 1
                    peptide_composition['O'] += 2
                    peptide_composition['H'] += 0
                    peptide_composition['S'] += 0
                
                if ion_type == 'x+1':
                    peptide_composition['C'] += 1
                    peptide_composition['O'] += 2
                    peptide_composition['H'] += 1
                    peptide_composition['S'] += 0
                    
                if ion_type == 'v': #x+1 ion - full side chain - CO
                    peptide_composition['C'] += 1 -1 - aa_side_chain_compositions[seq[0]]['C']
                    peptide_composition['O'] += 2 -1 - aa_side_chain_compositions[seq[0]]['O']
                    peptide_composition['H'] += 1 - 1 - aa_side_chain_compositions[seq[0]]['H']
                    peptide_composition['S'] += 0 -aa_side_chain_compositions[seq[0]]['S']
                    peptide_composition['N'] += -aa_side_chain_compositions[seq[0]]['N']
            
                if ion_type == 'w': #x ion - partial side chain - CHNO
                    peptide_composition['C'] += 1 - 1 - d_w_aa_side_chain_compositions[seq[0]]['C']
                    peptide_composition['O'] += 2 - 1- d_w_aa_side_chain_compositions[seq[0]]['O']
                    peptide_composition['H'] += 0 - 1 - d_w_aa_side_chain_compositions[seq[0]]['H']
                    peptide_composition['S'] += 0 - d_w_aa_side_chain_compositions[seq[0]]['S']
                    peptide_composition['N'] += -1 - d_w_aa_side_chain_compositions[seq[0]]['N']
            
                if ion_type == 'y':
                    peptide_composition['C'] += 0
                    peptide_composition['O'] += 1
                    peptide_composition['H'] += 2
                    peptide_composition['S'] += 0
                
                if ion_type == 'y-H2O':
                    peptide_composition['C'] += 0
                    peptide_composition['O'] += 0
                    peptide_composition['H'] += 0
                    peptide_composition['S'] += 0
            
                if ion_type == 'y-NH3':
                    peptide_composition['C'] += 0
                    peptide_composition['O'] += 1
                    peptide_composition['H'] += -1
                    peptide_composition['N'] += -1
                    peptide_composition['S'] += 0
                
                if ion_type == 'y-1':
                    peptide_composition['C'] += 0
                    peptide_composition['O'] += 1
                    peptide_composition['H'] += 2
                    peptide_composition['S'] += 0
                        
                if ion_type == 'y-2':
                    peptide_composition['C'] += 0
                    peptide_composition['O'] += 1
                    peptide_composition['H'] += 1
                    peptide_composition['S'] += 0
                
                if ion_type == 'z':
                    peptide_composition['C'] += 0
                    peptide_composition['O'] += 1
                    peptide_composition['H'] += 0
                    peptide_composition['N'] += -1
                    peptide_composition['S'] += 0
                        
                C.append(peptide_composition['C'])
                H.append(peptide_composition['H'])
                O.append(peptide_composition['O'])
                N.append(peptide_composition['N'])
                S.append(peptide_composition['S'])
                
            return (C,H,O,N,S)
            
            
            
        def ppm_difference(observed, theoretical):

            # Calculate ppm difference
            ppm_difference = (((observed-theoretical)/theoretical)*1000000)
            ppm_difference = round(ppm_difference, 2)

            return (ppm_difference)

        def match_fragments(deconvoluted_masses, theoretical_fragments, ion_types, tic):

        
            fragment_matches = {'ion type': [], 'position': [], 'Raw Intensity': [], 'TIC adjusted': [],
                                'Observed mass': [], 'Theoretical mass': [], 'ppm error': [], 'Relative Ion Number': [], 'N/C-term': []}
        
            ppm_diff = lambda x, y: abs(x - y) / x * 1e6
        
            ion_types_set = set(ion_types)
            nc_term_conditions = {'a', 'a+1', 'd', 'b', 'c', 'a-NH3', 'b-H2O', 'b-H2O', 'b+H2O', 'b-NH3'}
        
            observed_masses, ion_types_list, positions = [], [], []
            ppm_errors, theoretical_masses, raw_intensities = [], [], []
            tic_adjusted, relative_ion_numbers, nc_terms = [], [], []
        
            for i in deconvoluted_masses.itertuples():
                for key, value in theoretical_fragments.items():
                    mask = (np.array(value) > 0) & (key in ion_types_set) & (abs(ppm_diff(i.Mass, np.array(value))) <= int(PPM_error))
        
                    indices = np.where(mask)[0]
                    observed_masses.extend([i.Mass] * len(indices))
                    ion_types_list.extend([key] * len(indices))
                    positions.extend(indices)
                    ppm_errors.extend(ppm_diff(i.Mass, np.array(value))[mask])
                    theoretical_masses.extend(np.array(value)[mask])
                    raw_intensities.extend([i.Intensity] * len(indices))
                    tic_adjusted.extend([i.Intensity / tic] * len(indices))
        
                    relative_ion_numbers.extend(
                        [len(protein_sequence) - index if key in {'y', 'y-H2O', 'y-NH3', 'y-1', 'y-2', 'x', 'x+1', 'v', 'w', 'z'}
                          else index for index in indices]
                        )
        
                    nc_term = 'N-terminal' if key in nc_term_conditions else 'C-terminal'
                    nc_terms.extend([nc_term] * len(indices))
        
            fragment_matches['Observed mass'].extend(observed_masses)
            fragment_matches['ion type'].extend(ion_types_list)
            fragment_matches['position'].extend(positions)
            fragment_matches['ppm error'].extend(ppm_errors)
            fragment_matches['Theoretical mass'].extend(theoretical_masses)
            fragment_matches['Raw Intensity'].extend(raw_intensities)
            fragment_matches['TIC adjusted'].extend(tic_adjusted)
            fragment_matches['Relative Ion Number'].extend(relative_ion_numbers)
            fragment_matches['N/C-term'].extend(nc_terms)
        
            return fragment_matches
            

        def find_seq_coverage(df, seq):

            # Calculate sequence coverage
            unique = df['Relative Ion Number'].nunique()
            seq_coverage = (unique/(len(seq)-1)*100)
            
            if seq_coverage > 100:
                seq_coverage= 100

            return round(seq_coverage)

        def fragments_explained(matching_fragments, deconvoluted_masses):

            # Calculate Fragments Explained
            frags_explained = (len(matching_fragments) /
                               len(deconvoluted_masses))*100
            if frags_explained > 100:
                frags_explained = 100

            return round(frags_explained)

        def fragment_abundance_plot(df, protein_sequence):

            # Color mapping
            color_map = {'a': '#B4D496', 'a+1': '#B4D496', 'a-NH3': '#B4D496',
                         'b': '#A0BEEE', 'b-NH3': '#A0BEEE', 'b-H2O': '#A0BEEE', 'b+H2O': '#A0BEEE',
                         'c': '#FF8A8A', 'x': '#00B050', 'x+1': '#00B050',
                         'y': '#3764AA', 'y-1': '#3764AA', 'y-2': '#3764AA', 'y-NH3': '#3764AA', 'y-H2O': '#3764AA',
                         'z': '#C00000', 'd': '#35E3E3', 'v': '#35E3E3', 'w': '#35E3E3'}
        
            fig = px.bar(df, x='Relative Ion Number', y='TIC adjusted Intensity', color='ion type',
                         category_orders={'ion type': list(color_map.keys())},
                         hover_data={'ion type': True, 'Ion Number': True, 'TIC adjusted Intensity': True},
                         color_discrete_map=color_map,
                         labels={'TIC adjusted Intensity': 'Relative Intensity', 'Relative Ion Number': 'Backbone Position'})
        
            # Format Plot
            fig.update_xaxes(showline=True, linewidth=2, linecolor='black',
                             mirror=True, title_font=dict(size=40), tickfont=dict(size=40), range=[0, len(protein_sequence)])
            fig.update_yaxes(showline=True, linewidth=2, linecolor='black',
                             mirror=True, title_font=dict(size=40), tickfont=dict(size=40), showticklabels=True)
            
            # Combine layout updates
            fig.update_layout(title_font_family='Arial', title_x=0.5, title_font=dict(size=40),
                              legend_font=dict(size=40), template='seaborn', plot_bgcolor='white')
            
            # Combine trace modifications
            fig.for_each_trace(lambda t: t.update(name=f"{t.name}         "))
            fig.data = fig.data[::-1]
            fig.layout.legend.traceorder = 'reversed'
        
            return fig

        def fragment_abundance_plot_NC_term(df, protein_sequence):

            fig = px.bar(df, x='Relative Ion Number', y='TIC adjusted Intensity', color='N/C-term',

                         category_orders={'': ['N-terminal', 'C-terminal']},

                         hover_data={
                             'Ion Number': True, 'Relative Ion Number': False, 'TIC adjusted Intensity': True},

                         color_discrete_map={'N-terminal': '#008B8B',
                                             'C-terminal': '#CD853F'},

                         labels={'TIC adjusted Intensity': 'TIC Adjusted Intensity', 'Relative Ion Number': 'Backbone Position'})

            # Format Plot
            fig.update_xaxes(showline=True, linewidth=2, linecolor='black',
                             mirror=True, title_font=dict(size=40), tickfont=dict(size=40), range=[0, len(protein_sequence)])
            fig.update_yaxes(showline=True, linewidth=2, linecolor='black',
                             mirror=True, title_font=dict(size=40), tickfont=dict(size=40), showticklabels=True)
            
            # Combine layout updates
            fig.update_layout(title_font_family='Arial', title_x=0.5, title_font=dict(size=40),
                              legend_font=dict(size=40), template='seaborn', plot_bgcolor='white')
            
            # Combine trace modifications
            fig.for_each_trace(lambda t: t.update(name=f"{t.name}         "))
            fig.data = fig.data[::-1]
            fig.layout.legend.traceorder = 'reversed'

            return fig

        def get_results(seq_coverage, frags_explained, protein_monoisotopic_mass, seq_len):

            # Print monoisotopic mass, sequence coverage, fragments explained, sequence length
            Results = 'Results \n \n Protein Mass: ' + str(protein_monoisotopic_mass) + \
                '\t Sequence Length: ' + str(seq_len) + ' A.A.' + \
                                                        '\t Sequence Coverage: ' + str(seq_coverage) + '%' + \
                '\t Fragments Explained: ' + str(frags_explained) + '% \n'

            return Results

        def sequence_coverage_plot(protein_sequence, matching_fragments):
            # Separate protein into list
            protein_aa_list = ['<b>' + aa + '<b>' for aa in protein_sequence]
            
            # Set X and Y coordinates
            y = [i for i in range(math.ceil(len(protein_aa_list)/25))]
            y.reverse()
            y_list = []
            for i in y:
                for j in range(25):
                    y_list.append(i)
            x_list = [(i+1) for i in range(25)] * \
                (math.ceil(len(protein_aa_list)/25))
            aa_number = [(i+1) for i in range(len(protein_sequence))]
            
            # Create dataframe to hold coordinates
            
            coordinates = pd.DataFrame()
            coordinates['protein'] = protein_aa_list
            coordinates['x'] = x_list[0:len(protein_aa_list)]
            coordinates['y'] = y_list[0:len(protein_aa_list)]

            y_max = coordinates['y'].max()
            y_min = coordinates['y'].min()
            
            font_size = 44-y_max*1.5
            trace_data = []
            
            protein_sequence_trace = go.Scatter(
                x=x_list,
                y=y_list,
                mode="text",
                textfont = dict(family = 'Arial', size = font_size),
                text= protein_aa_list ,
                textposition="middle center",
                name='ion type',
                hovertext=aa_number,
                hoverinfo='text',
                hovertemplate='%{hovertext}')

            trace_data.append(protein_sequence_trace)
            

        ###############################################################################

            keys = ['a+1', 'a', 'a-NH3', 'b', 'b-NH3', 'b-H2O', 'b+H2O', 'c', 'x', 'x+1', 'y', 'y-NH3', 'y-H2O', 'y-1', 'y-2', 'z', 'd', 'v', 'w']
            
            ion_dict = {i: {'ion type': [], 'ion mass': [], 'backbone position': [], 'ion location x': [], 'ion location y': []} for i in keys}
            
            key_set = set(keys)
            
            for index, row in enumerate(matching_fragments.itertuples()):
                for i in coordinates.itertuples():
                    key = row._1
                    if key in key_set:
                        if (i.Index + 1) == row._8:
                            ion_type = ion_dict[key]
                            ion_type['ion type'].append(key)
                            ion_type['ion mass'].append(row._5)
                            ion_type['ion location x'].append(coordinates['x'][i.Index])
                            ion_type['ion location y'].append(coordinates['y'][i.Index])
                            if key in ['a+1', 'a', 'a-NH3', 'b', 'b-NH3', 'b-H2O', 'b+H2O', 'c']:
                                ion_type['backbone position'].append(row._8)
                            else:
                                ion_type['backbone position'].append(row._2)

                                    
            ###############################################################################
            global ion_df_dict
            ion_df_dict = {f'{i}' : pd.DataFrame() for i in keys}
            
            for k,v in ion_df_dict.items():
                
                    v['ion type'] = ion_dict[k]['ion type']
                    v['ion position'] = ion_dict[k]['backbone position']
                    v['ion mass'] = ion_dict[k]['ion mass']

            ###############################################################################
            color_dict = {'a': '#008000', 'a+1': '#008000', 'a-NH3': '#008000','x': '#008000', 'x+1': '#008000',
                          'b': '#0000FF', 'b-NH3': '#0000FF','b-H2O': '#0000FF', 'b+H2O': '#0000FF',
                          'y': '#0000FF', 'y-NH3': '#0000FF', 'y-H2O': '#0000FF','y-1': '#0000FF', 'y-2': '#0000FF',
                          'c': '#FF0000', 'z': '#FF0000','d': '#35E3E3', 'v': '#35E3E3', 'w': '#35E3E3'}
            
            for i in keys[::-1]:
                
                
                if i in ['a','a+1','a-NH3','b','b-NH3','b-H2O','b+H2O','c','d']:
                    x_data = [item + increment for item in  ion_dict[i]['ion location x'] for increment in [0.5, 0.5, 0]]
                    y_data = [item + increment for item in  ion_dict[i]['ion location y'] for increment in [0,0.2,0.4]]
                    
                else:
                    x_data = [item + increment for item in  ion_dict[i]['ion location x'] for increment in [0.5, 0.5, 1.0]]
                    y_data = [item + increment for item in  ion_dict[i]['ion location y'] for increment in [0,-0.2,-0.4]]

                for k in range(0, len(x_data),3):
                    x = x_data[k:k+3]
                    y = y_data[k:k+3]
                    
                    trace = (go.Scatter(
                        x = x,
                        y = y,
                        mode = 'markers',
                        name = i,
                        marker = dict(
                            color = 'rgba(0, 0, 255, 0)',
                            size = 10),
                        hoverinfo = 'skip'))
                
                    line_trace = (go.Scatter(
                        x = x,
                        y = y,
                        mode = 'lines',
                        line = dict(color = 'rgba(183, 183, 183, 0.8)',
                                width = 10-y_max*0.5),
                        hoverinfo = 'skip'))
                    
                    trace_data.append(trace)
                    trace_data.append(line_trace)
            
            
            for i in keys[::-1]:
                
                    
                    if i in ['a','a+1','a-NH3','b','b-NH3','b-H2O','b+H2O','c','d']:
                        x_data = [item + increment for item in  ion_dict[i]['ion location x'] for increment in [0.5, 0.5, 0]]
                        y_data = [item + increment for item in  ion_dict[i]['ion location y'] for increment in [0,0.2,0.4]]
                        
                    else:
                        x_data = [item + increment for item in  ion_dict[i]['ion location x'] for increment in [0.5, 0.5, 1.0]]
                        y_data = [item + increment for item in  ion_dict[i]['ion location y'] for increment in [0,-0.2,-0.4]]
                        
                    x_scale = y_scale = 0.0  # Default values for cases not covered in conditions
                    
                    if i in ['a', 'a+1', 'a-NH3']:
                        y_scale = 0.2
                    elif i in ['b', 'b-NH3', 'b-H2O', 'b+H2O']:
                        x_scale, y_scale = -0.167, 0.266
                    elif i == 'c':
                        x_scale, y_scale = -0.333, 0.332
                    elif i == 'd':
                        x_scale, y_scale = -0.5, 0.4
                    elif i in ['x', 'x+1']:
                        y_scale = -0.2
                    elif i in ['y', 'y-NH3', 'y-H2O', 'y-1', 'y-2']:
                        x_scale, y_scale = 0.167, -0.266
                    elif i == 'z':
                        x_scale, y_scale = 0.333, -0.332
                    elif i in ['v', 'w']:
                        x_scale, y_scale = 0.5, -0.4

                    for k in range(0, len(x_data),3):
                        count = int(k/3)
                        x = x_data[k:k+3]
                        y = y_data[k:k+3]
                
                        point_trace = (go.Scatter(
                        x =[x[0]+x_scale],
                        y = [y[0]+y_scale],
                        mode = 'markers',
                        name = i,
                        marker = dict(
                            symbol = 'circle',
                            color = color_dict[i],
                            size = font_size*0.1,
                            line = dict(
                                color = color_dict[i],
                                width = 8)),
                        customdata = [(ion_df_dict[i]['ion type'][count], ion_df_dict[i]['ion position'][count], ion_df_dict[i]['ion mass'][count])],
                        hovertemplate = ('Ion Type: %{customdata[0]}<br>' +
                                        'Ion Number: %{customdata[1]}<br>' +
                                        'Ion Mass: %{customdata[2]}<br><extra></extra>'),
                        hovertext = [(ion_dict[i]['ion type'][count], ion_dict[i]['backbone position'][count],
                                      ion_dict[i]['ion mass'][count])],
                        hoverinfo = 'text')
                        )
                    
                        trace_data.append(point_trace)
                   
                
            fig = go.Figure(data = trace_data)
            
            def check_value(value, my_dict):
                for key, val in my_dict.items():
                    if abs(value - val) <= 0.01:
                        return key
                return False
            
            color_count = 0
            

            for i,row in coordinates.iterrows():
                if (i+1) in PTM_position_list:
                             
                    fig.add_shape(type = 'rect',
                                      xref = 'x', yref = 'y',
                                      x0 = coordinates['x'][i] - 0.3, y0 = coordinates['y'][i] + 0.15,
                                      x1 = coordinates['x'][i] + 0.3, y1 = coordinates['y'][i] - 0.15,
                                      line = dict(color = '#000000',
                                                  width = 3,
                                                  ), fillcolor = px.colors.qualitative.Plotly[color_count])
                    color_count +=1
            

        ###############################################################################
     
            fig.update_layout(template='simple_white',
                              showlegend = False,
                               yaxis = dict(range = [y_min-(y_max*0.2),y_max*1.2]),
                               xaxis = dict(range = [0,26]))
            fig.update_xaxes(ticks='', showticklabels=False, showline=False)
            fig.update_yaxes(ticks='', showticklabels=False, showline=False)
            
            
            return fig

        def get_variables():

            seq = seq_entry.get("1.0", "end-1c")
            tic2 = tic_var.get()
            activation = menu.get().upper()
            PTM_position = PTM_position_var.get()
            PTM_mass = PTM_mass_var.get()
            Ligand_mass = ligand_mass_var.get()
            Apo_Holo_menu = holo_menu.get()
            Neutral_loss_menu = neutral_loss_menu.get()
            global PPM_error
            PPM_error = PPM_var.get()

            return seq, tic2, activation, PTM_position, PTM_mass, Ligand_mass, Apo_Holo_menu, Neutral_loss_menu, PPM_error
        

        def file_read(file):

            with file as f:
                for l in f:
                    if l.startswith("Mass") or l.startswith('mass'): 
                        columns = l.split(',')
                        new_columns = []
                        for x in columns:
                             new_columns.append(x.strip())
                        break
                    
                capitalized = [word.capitalize() for word in new_columns]
                deconvoluted_masses = pd.read_csv(
                    file, float_precision=None, names=capitalized)
                
            deconvoluted_masses = deconvoluted_masses.sort_values(
                'Intensity', ascending=False).drop_duplicates('Mass').sort_index()
            
            return deconvoluted_masses

        # Code to run above functions on button press

        def execute():
            
            # Get variables from app entry boxes
            seq, tic2, activation, PTM_position, PTM_mass, Ligand_mass, Apo_Holo_menu, Neutral_loss_menu, PPM_error = get_variables()
            if PTM_position != '':
                global PTM_position_list
                PTM_position_list = list(map(int,PTM_position.split(',')))

            else:
                PTM_position_list = []


            
            # Open File
            file = open(filename, 'r', encoding='utf-8-sig')
            global protein_sequence
            protein_sequence = ''.join(re.split(r"\s", seq))

            # Read in deconvoluted mass list, remove duplicates, keep highest intensity value
            deconvoluted_masses = file_read(file)
            # Define ion types in search according to activation entry (HCD,CID,ETD,EThcD,UVPD)
            ion_types = user_ions(activation)

            # Normalize according to TIC entry
            tic = user_TIC(tic2)

            # Adjust protein mass and fragments according to PTM entry and Ligand Masses
            ptm_dict = user_PTMs(PTM_position, PTM_mass)
            Ligand_mass = user_Ligand(Ligand_mass)
            
            if Apo_Holo_menu == 'Apo':
                protein_monoisotopic_mass = theoretical_mass(
                    protein_sequence, ptm_dict, [0])
            else:
                protein_monoisotopic_mass = theoretical_mass(
                    protein_sequence, ptm_dict, Ligand_mass)
            
            
            theoretical_fragments = theoretical_fragment_calc(
                protein_sequence, ptm_dict, Ligand_mass)[0]
            ligand_fragments = theoretical_fragment_calc(
                protein_sequence, ptm_dict, Ligand_mass)[1]
          
            neutral_loss_fragments = theoretical_fragment_calc(
                protein_sequence, ptm_dict, Ligand_mass)[2]
            
            global theoretical_fragment_df
            theoretical_fragment_df = pd.DataFrame.from_dict(theoretical_fragments)

            global neutral_loss_fragment_df
            neutral_loss_fragment_df = pd.DataFrame.from_dict(neutral_loss_fragments)

            # Match fragments
            global matching_fragments
            
            # Determine if Apo/Holo
            if Apo_Holo_menu == 'Apo' and Neutral_loss_menu == 'No':
                
                matching_fragments = pd.DataFrame.from_dict(match_fragments(
                    deconvoluted_masses, theoretical_fragments, ion_types, tic))
                matching_fragments = matching_fragments.sort_values(
                    'TIC adjusted', ascending=False).drop_duplicates('Theoretical mass').sort_index()
                
            if Apo_Holo_menu == 'Apo' and Neutral_loss_menu == 'Yes':
                
                matching_fragments = pd.DataFrame.from_dict(match_fragments(
                    deconvoluted_masses, neutral_loss_fragments, ion_types, tic))
                matching_fragments = matching_fragments.sort_values(
                    'TIC adjusted', ascending=False).drop_duplicates('Theoretical mass').sort_index()
                
            if Apo_Holo_menu == 'Holo':
                matching_fragments = pd.DataFrame.from_dict(match_fragments(
                    deconvoluted_masses, ligand_fragments, ion_types, tic))
                matching_fragments = matching_fragments.sort_values(
                    'TIC adjusted', ascending=False).drop_duplicates('Theoretical mass').sort_index()

            C,H,O,N,S = chemical_formula(matching_fragments)
            
            matching_fragments['C'] = C
            matching_fragments['H'] = H
            matching_fragments['O'] = O
            matching_fragments['N'] = N
            matching_fragments['S'] = S
            
            # Format table for display
            matching_fragments.rename(columns={'TIC adjusted': 'TIC adjusted Intensity', 
                                               'position': 'Ion Number',
                                               'Relative Ion Number': 'Relative Ion Number'}, inplace= True)
            matching_fragments['ppm error'] = ['%.2f' %
                                               i for i in matching_fragments['ppm error']]
            matching_fragments['Raw Intensity'] = ['%.2f' %
                i for i in matching_fragments['Raw Intensity']]
            matching_fragments['TIC adjusted Intensity'] = [
                '%.2e' % i for i in matching_fragments['TIC adjusted Intensity']]
            
            # Display table (Set table size, location, float precision, select columns)
            root = tk.Tk()
            ws = root.winfo_screenwidth()
            hs = root.winfo_screenheight()
            w = ws/2
            h = hs/2.3
            x = (ws) - (ws/2) - (ws/128)
            y = (hs) - (hs)
            root.geometry('%dx%d+%d+%d' % (w, h, x, y))
            root.title('Matching Fragments')
            frame = tk.Frame(root)
            frame.pack(fill='both', expand=True)
            options = {'floatprecision': 4, 'cellwidth': 100}
            pt = Table(frame, dataframe=matching_fragments.iloc[:, list(range(0, 8)) + list(range(9, 14))])
            config.apply_options(options, pt)
            pt.show()

            # Print Results in main window
            seq_coverage = find_seq_coverage(
                matching_fragments, protein_sequence)
            frags_explained = fragments_explained(
                matching_fragments, deconvoluted_masses)
            seq_len = len(protein_sequence)
            results = get_results(
                seq_coverage, frags_explained, protein_monoisotopic_mass, seq_len)
            results_label = tk.Label(self, text=results)
            results_label.config(font=('helvetica', 12))
            results_label.grid(column=0, row=13, columnspan = 6)
            self.rowconfigure(13, weight=3)

            # Display table with only ion type, position, and intensity
            tidy_table = pd.pivot_table(matching_fragments[['ion type', 'TIC adjusted Intensity', 'Relative Ion Number']].copy(),
                                        values='TIC adjusted Intensity', index='Relative Ion Number', columns = 'ion type', fill_value = 0)
            protein_len = pd.DataFrame(
                {'Relative Ion Number': [x for x in range(1, len(protein_sequence))]})

            tidy_table = tidy_table.merge(protein_len, on='Relative Ion Number', how='right').fillna(0)
            fragment_types = ['a', 'a+1', 'b', 'c', 'x','x+1','y','y-1','y-2','z']
            for i in fragment_types:
                if i not in tidy_table:
                    tidy_table[i] = 0
            tidy_table = tidy_table[['Relative Ion Number', 'a', 'a+1', 'b', 'c', 'x','x+1','y','y-1','y-2','z']]
            
            # Format tidy table for display
            format_mapping = {'a': '{:.2e}', 'a+1': '{:.2e}', 'b': '{:.2e}', 'c': '{:.2e}',
                              'x': '{:.2e}', 'x+1': '{:.2e}', 'y': '{:.2e}', 'y-1': '{:.2e}', 'y-2': '{:.2e}','z': '{:.2e}'}
            for key, value in format_mapping.items():
                tidy_table[key] = tidy_table[key].apply(value.format)
            
            # Display table
            root2 = tk.Tk()
            x2 = (ws) - (ws/2) - (ws/128)
            y2 = (hs) - (hs/1.9)
            root2.geometry('%dx%d+%d+%d' % (w, h, x2, y2))
            root2.title("Fragment's position relative to N-Term")
            frame = tk.Frame(root2)
            frame.pack(fill='both', expand=True)
            pt2 = Table(frame, dataframe=tidy_table)
            config.apply_options(options, pt2)
            pt2.show()
            

            
            matching_fragments['ppm error'] = matching_fragments['ppm error'].astype(float)
            matching_fragments['Raw Intensity'] = matching_fragments['Raw Intensity'].astype(float)
            matching_fragments['TIC adjusted Intensity'] = matching_fragments['TIC adjusted Intensity'].astype(float)

        def graph_all_ions():

            # Fragment Abundance Plot with all ion types
            frag_plot_fig = fragment_abundance_plot(
                matching_fragments, protein_sequence)
            plot(frag_plot_fig)

        def graph_NC_terminal():

            # Fragment Abundance Plot with only N/C terminal
            NC_plot_fig = fragment_abundance_plot_NC_term(
                matching_fragments, protein_sequence)
            plot(NC_plot_fig)

        def graph_seq_coverage():

            # Sequence Coverage Plot
            
            seq_fig = sequence_coverage_plot(
                protein_sequence, matching_fragments)
            plot(seq_fig)

            
        # Go button
        go_button = tk.Button(self, text='Go', command=execute,
                              bg='brown', fg='white', font=('helvetica', 12, 'bold'))
        self.rowconfigure(12, weight=2)
        go_button.grid(column=1, row=12, columnspan=4, sticky = tk.NSEW, 
                       padx = 200, pady = 50)
        

        # Graph all ions button
        graph_button = tk.Button(self, text='Graph All Ion Types', 
                                 command=graph_all_ions, bg='brown', 
                                 fg='white', font=('helvetica', 9, 'bold'))
        graph_button.grid(column=0, row=22)

        # Graph N/C terminal button
        graph2_button = tk.Button(self, text='Graph N/C-Terminal', 
                                  command=graph_NC_terminal, bg='brown', 
                                  fg='white', font=('helvetica', 9, 'bold'))
        graph2_button.grid(column=2, row=22, columnspan=2)

        # Graph Seq Coverage button

        graph3_button = tk.Button(self, text='Graph Sequence Coverage',
                                  command=graph_seq_coverage, bg='brown',
                                  fg='white', font=('helvetica', 9, 'bold'))
        graph3_button.grid(column=5, row=22)
        

        # Next Page button
        next_page_button = tk.Button(self, text='Charge Site Analysis', 
                                     command=lambda: controller.show_frame(ResultsPage), bg='brown', 
                                     fg='white', font=('helvetica', 9, 'bold'))
        next_page_button.grid(column=5, row=26)


# Charge Site Analysis Page
class ResultsPage(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.controller = controller
        canvas2 = tk.Canvas(self, width=1200, height=1000,  relief = 'raised')
        canvas2.pack()

        # Page Title
        title_label = tk.Label(self, text='Charge State Analysis')
        title_label.config(font=('helvetica', 14))
        canvas2.create_window(400, 25, window=title_label)

        # Open Deconvoluted Mass + Intensity List
        csv_label = tk.Label(
            self, text='Upload Deconvoluted Mass + Intensity List:')
        csv_label.config(font=('helvetica', 10, 'bold'))
        canvas2.create_window(73, 100, window=csv_label, anchor='w')

        def clicked():
            global filename
            filename = askopenfilename(filetypes=(("CSV Files", "*.csv"),))
            file_label['text'] = filename
        open_csv_button = tk.Button(
            self, text="Upload", command=clicked, bg='brown', fg='white', font=('helvetica', 9, 'bold'))
        canvas2.create_window(400, 100, window=open_csv_button)
        file_label = tk.Label(self, text='')
        canvas2.create_window(73, 140, window=file_label, anchor='w')

        # Open Deconvoluted Mass + Charge List
        xls_label = tk.Label(
            self, text='Upload Deconvoluted Mass + Charge Information:')
        xls_label.config(font=('helvetica', 10, 'bold'))
        canvas2.create_window(73, 200, window=xls_label, anchor='w')

        def clicked2():
            global charge_state_filename
            charge_state_filename = askopenfilename(
                filetypes=(("XLS Files", "*.xls"),))
            charge_state_file_label['text'] = charge_state_filename
        open_xls_button = tk.Button(
            self, text="Upload", command=clicked2, bg='brown', fg='white', font=('helvetica', 9, 'bold'))
        canvas2.create_window(450, 200, window=open_xls_button)
        charge_state_file_label = tk.Label(self, text='')
        canvas2.create_window(
            73, 240, window=charge_state_file_label, anchor='w')

        # Protein Sequence Entry
        seq_label = tk.Label(self, text='Protein Sequence: ')
        seq_label.config(font=('helvetica', 10, 'bold'))
        canvas2.create_window(73, 280, window=seq_label,  anchor='w')
        seq_entry = tk.Text (self, height=10, width=50, wrap = tk.WORD)
        seq_entry.configure(font=('helvetica', 10))
        canvas2.create_window(220, 270, window=seq_entry,  anchor='nw')

        # PTM position Entry
        PTM_position_var = tk.StringVar()
        PTM_position_label = tk.Label(
            self, text='PTM position(s): \n (comma separated)')
        PTM_position_label.config(font=('helvetica', 10, 'bold'))
        canvas2.create_window(64, 460, window=PTM_position_label,  anchor='w')
        PTM_position_entry = tk.Entry (self, textvariable=PTM_position_var, width = 58)
        canvas2.create_window(220, 462, window=PTM_position_entry,  anchor='w')

        # PTM mass Entry
        PTM_mass_var = tk.StringVar()
        PTM_mass_label = tk.Label(
            self, text='PTM mass(es): \n (comma separated)')
        PTM_mass_label.config(font=('helvetica', 10, 'bold'))
        canvas2.create_window(58, 520, window=PTM_mass_label,  anchor='w')
        PTM_mass_entry = tk.Entry (self, textvariable=PTM_mass_var, width= 58)
        canvas2.create_window(220, 522, window=PTM_mass_entry, anchor='w')

        # PPM Error Entry
        PPM_var = tk.StringVar(value='10')
        PPM_label = tk.Label(self, text='PPM error:')
        PPM_label.config(font=('helvetica', 10, 'bold'))
        canvas2.create_window(73, 580, window=PPM_label,  anchor='w')
        PPM_entry = tk.Entry (self, textvariable=PPM_var, width=58)
        canvas2.create_window(220, 582, window=PPM_entry,  anchor='w')

        # Functions (same as previous page)

        def user_PTMs(x, y):

            if x == '':
                x = [0]

            else:
                x = [int(x) for x in x.split(',')]

            if y == '':
                y = [0]

            else:
                y = [float(x) for x in y.split(',')]

            PTM_dict = {x[i]: y[i] for i in range(len(x))}

            return (PTM_dict)

        def theoretical_mass(seq, ptm_dict):
            AA_masses = {"A": 71.03711, "C": 103.00919, "D": 115.02694, "E": 129.04259,
                         "F": 147.06841, "G": 57.02146, "H": 137.05891, "I" : 113.08406,
                         "K": 128.09496, "L": 113.08406, "M": 131.04049, "N" : 114.04293,
                         "P": 97.05276, "Q": 128.05858, "R": 156.10111, "S" : 87.03203,
                         "T": 101.04768, "V": 99.06841, "W": 186.07931, "Y" : 163.06333}

            protein_aa_list = [aa for aa in seq]  # Split sequence into list

            protein_monoisotopic_mass = 0

            for i in protein_aa_list:
                if i in AA_masses:
                    protein_monoisotopic_mass += AA_masses[i]

            for value in ptm_dict.values():
                protein_monoisotopic_mass += value

            protein_monoisotopic_mass += 18.01

            return round(protein_monoisotopic_mass, 2)

        def theoretical_fragment_calc(seq, ptm_dict):

            AA_masses = {"A": 71.03711, "C": 103.00919, "D": 115.02694, "E": 129.04259,
                         "F": 147.06841, "G": 57.02146, "H": 137.05891, "I" : 113.08406,
                         "K": 128.09496, "L": 113.08406, "M": 131.04049, "N" : 114.04293,
                         "P": 97.05276, "Q": 128.05858, "R": 156.10111, "S" : 87.03203,
                         "T": 101.04768, "V": 99.06841, "W": 186.07931, "Y" : 163.06333}

            protein_aa_list = [aa for aa in seq]  # Split sequence into list

            mass_list = []

            for i in protein_aa_list:
                if i in AA_masses:
                    mass_list.append(AA_masses[i])

            for index, value in enumerate(mass_list):
                for k, v in ptm_dict.items():
                    if index == (k-1):
                        mass_list[index] += v

            fragment_dict = {'a': [], 'a+1': [], 'b': [], 'c': [], 'x':[], 'x+1':[], 'y':[], 'y-1':[], 'y-2':[], 'z':[]}   ###Dictionary of lists for fragment ion types

            fragment_mass = 0
            N_term_fragment_list = []  # List for N-term fragments

            for index, value in enumerate(mass_list):
                N_term_fragment_list.append(fragment_mass)
                fragment_mass += mass_list[index]

            fragment_mass = 0
            C_term_fragment_list = []  # List for C-term fragments

            for index, value in reversed(list(enumerate(mass_list))):
                C_term_fragment_list.append(fragment_mass)
                fragment_mass += mass_list[index]

            for i in N_term_fragment_list:  # Add fragments to dictionary with appropriate mass shifts
                fragment_dict['a'].append(i-27.9949)
                fragment_dict['a+1'].append(i-26.987075)
                fragment_dict['b'].append(i)
                fragment_dict['c'].append(i+17.026549)

            for i in C_term_fragment_list:  # Add fragments to dictionary with appropriate mass shifts
                fragment_dict['x'].append(i+43.989830)
                fragment_dict['x+1'].append(i+44.997655)
                fragment_dict['y'].append(i+18.010565)
                fragment_dict['y-1'].append(i+17.00274)
                fragment_dict['y-2'].append(i+15.994915)
                fragment_dict['z'].append(i+1.991841)

            return fragment_dict

        def ppm_difference(observed, theoretical):

            ppm_difference = (((observed-theoretical)/theoretical)*1000000)

            return (ppm_difference)

        def match_fragments(deconvoluted_masses, theoretical_fragments, ion_types):


            fragment_matches = {'ion type': [], 'position': [], 'Raw Intensity': [], 'TIC adjusted': [], 'Observed mass':[], 'theoretical mass': [], 'ppm error':[], 'Relative Ion Number':[], 'N/C-term':[]}
            k = 0

            for i in deconvoluted_masses.itertuples():
                for key, value in theoretical_fragments.items():
                    for j in value:
                        if j > 0:
                            if key in ion_types:
                                if abs(ppm_difference(i.Mass, j)) <= int(PPM_error):
                                    fragment_matches['Observed mass'].append(
                                        i.Mass)
                                    fragment_matches['ion type'].append(key)
                                    fragment_matches['position'].append(
                                        theoretical_fragments[key].index(j))
                                    fragment_matches['ppm error'].append(
                                        ppm_difference(i.Mass, j))
                                    fragment_matches['theoretical mass'].append(
                                        j)
                                    fragment_matches['Raw Intensity'].append(
                                        i.Intensity)
                                    fragment_matches['TIC adjusted'].append(
                                        (i.Intensity)/1)
                                    if fragment_matches['ion type'][k] == 'y' or fragment_matches['ion type'][k] == 'y-1' or fragment_matches['ion type'][k] == 'y-2' \
                                            or fragment_matches['ion type'][k] == 'x' or fragment_matches['ion type'][k] == 'x+1' or fragment_matches['ion type'][k] == 'z':
                                        fragment_matches['Relative Ion Number'].append(
                                            len(protein_sequence)-theoretical_fragments[key].index(j))
                                    else:
                                        fragment_matches['Relative Ion Number'].append(
                                            theoretical_fragments[key].index(j))
                                    if key == 'a' or key == 'a+1' or key == 'b' or key == 'c':
                                        fragment_matches['N/C-term'].append(
                                            'N-terminal')
                                    else:
                                        fragment_matches['N/C-term'].append(
                                            'C-terminal')
                                    k += 1

            return (fragment_matches)

        def execute():
            # Get variables from entry boxes
            seq = seq_entry.get("1.0", "end-1c")
            PTM_position = PTM_position_var.get()
            PTM_mass = PTM_mass_var.get()
            global PPM_error
            PPM_error = PPM_var.get()
            file = open(filename, 'r')
            global protein_sequence
            protein_sequence = seq

            # Read in deconvoluted mass list, remove duplicates
            with file as f:
                for l in f:
                    if l.startswith("Mass"):
                        columns = l.split(',')
                        new_columns = []
                        for x in columns:
                             new_columns.append(x.strip())
                        break
                deconvoluted_masses = pd.read_csv(
                    file, float_precision=None, names=new_columns)
            deconvoluted_masses = deconvoluted_masses.sort_values(
                'Intensity', ascending=False).drop_duplicates('Mass').sort_index()

            # Define ion types in search
            ion_types = ['a', 'a+1', 'b', 'c', 'x','x+1','y','y-1','y-2','z']

            # Adjust protein mass and fragments according to PTM entry
            ptm_dict = user_PTMs(PTM_position, PTM_mass)
            protein_monoisotopic_mass = theoretical_mass(
                protein_sequence, ptm_dict)
            theoretical_fragments = theoretical_fragment_calc(
                protein_sequence, ptm_dict)

            # Match fragments
            global matching_fragments
            matching_fragments = pd.DataFrame.from_dict(match_fragments(
                deconvoluted_masses, theoretical_fragments, ion_types))

            ##################################################################

            charge_state_file = pd.read_excel(charge_state_filename)

            charge_state_file['No.'] = charge_state_file['No.'].fillna(
                method='ffill')

            # extract deconvoluted monoisotopic masses
            mass_only = charge_state_file.dropna().iloc[:, [0, 1]]

            # tidy charge state data
            charge_site = charge_state_file.loc[charge_state_file['RT Range'].isnull()].iloc[:, [0, 1, 4, 5]].dropna()
            charge_site.columns = charge_site.iloc[0]
            charge_site = charge_site.drop(
                charge_site.loc[charge_site['Charge State'] == 'Charge State'].index).rename(columns={1.0: 'No.'})

            # merge tables
            merge_table = pd.merge(charge_site, mass_only, on='No.', how= 'outer').astype(float)
            merge_table = merge_table[[
                'No.', 'Monoisotopic Mass', 'Mostabund m/z', 'Charge State', 'Charge Normalized Intensity']]
            merge_table['Corrected Charge Normalized Intensity'] = merge_table['Charge Normalized Intensity'] / \
                merge_table['Charge State']
            merge_table['Charge State'] = merge_table['Charge State'].astype(
                int)

            # Filter so only contains a,a+1,x,x+1 ions
            a_x = ['a', 'a+1', 'x', 'x+1']
            combined_table = pd.merge(matching_fragments, merge_table, left_on='Observed mass', right_on= 'Monoisotopic Mass', how = 'outer'
                                      ).iloc[:, [0, 1, 7, 10,11,12,13,14]]
            combined_table2 = combined_table.loc[combined_table['ion type'].isin(
                a_x)]
            combined_table2['Charge State'] = pd.Categorical(
                combined_table2['Charge State'])

            # Make plot with a ions
            color_map = {'1': '#e6194B' ,'2': '#3f51b5','3': '#3cb44b','4': '#4363d8','5': '#f58231',
                         '6': '#911eb4','7': '#42d4f4','8': '#f032e6','9': '#000075','10': '#9A6324',
                         '11': '#e6194B','12': '#3f51b5','13': '#3cb44b','14': '#4363d8','15': '#f58231',
                         '16': '#911eb4','17': '#42d4f4','18': '#f032e6','19': '#000075','20': '#9A6324',
                         '21': '#e6194B','22': '#3f51b5','23': '#3cb44b','24': '#4363d8','25': '#f58231',
                         '26': '#911eb4','27': '#42d4f4','28': '#f032e6','29': '#000075','30': '#9A6324'}
            
            a = ['a', 'a+1']
            a_table = combined_table2[combined_table2['ion type'].isin(a)]
            a_table['Percent Intensity'] = a_table['Corrected Charge Normalized Intensity'] / \
                a_table.groupby('Relative Ion Number')[
                'Corrected Charge Normalized Intensity'].transform('sum')
            a_table = a_table.sort_values(by = 'Charge State')
            fig1 = px.bar(a_table, x='Relative Ion Number', y='Percent Intensity',
                          color='Charge State', barmode = 'stack', title = 'Charge Site Analysis with a and a+1 Ions',
                          labels={'Relative Ion Number': 'Backbone Position'}, color_discrete_map=color_map)
            
            fig1.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True, title_font=dict(size=25), tickfont=dict(size=25))
            fig1.update_xaxes(range=[0, len(protein_sequence)])
            fig1.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True, title_font=dict(size=25), tickfont=dict(size=25), tickformat = '.1f')
            fig1.update_layout(template='seaborn', plot_bgcolor = 'white')
            fig1.update_traces(marker=dict(pattern_solidity=0.4, pattern_fgcolor = 'black'))
            fig1.update_layout(title_font_family='Arial', title_x=0.5, title_font = dict(size = 35), legend_font = dict(size = 25))

            plot(fig1)

            # Make plot with x ions
            x = ['x', 'x+1']
            x_table = combined_table2[combined_table2['ion type'].isin(x)]
            x_table['Percent Intensity'] = x_table['Corrected Charge Normalized Intensity'] / \
                x_table.groupby('Relative Ion Number')[
                'Corrected Charge Normalized Intensity'].transform('sum')
            x_table = x_table.sort_values(by = 'Charge State')
            x_table['Charge State'] = x_table['Charge State'].astype(str)
            fig2 = px.bar(x_table, x='Relative Ion Number', y='Percent Intensity',
                          color='Charge State', barmode = 'stack', title = 'Charge Site Analysis with x and x+1 Ions',
                          labels={'Relative Ion Number': 'Backbone Position'}, color_discrete_map=color_map)

            fig2.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True, title_font=dict(size=25), tickfont=dict(size=25))
            fig2.update_xaxes(range=[0, len(protein_sequence)])
            fig2.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True, title_font=dict(size=25), tickfont=dict(size=25), tickformat = '.1f')
            fig2.update_layout(template='seaborn', plot_bgcolor = 'white')
            fig2.update_traces(marker=dict(pattern_solidity=0.4, pattern_fgcolor = 'black'))
            fig2.update_layout(title_font_family='Arial', title_x=0.5, title_font = dict(size = 35), legend_font = dict(size = 25))

            plot(fig2)

        # Go button
        go_button = tk.Button(self, text='Go', command=execute,
                              bg='brown', fg='white', font=('helvetica', 9, 'bold'))
        canvas2.create_window(400, 750, window=go_button)

        # Next Page button
        last_page_button = tk.Button(self, text='Fragment Intensity', command=lambda: controller.show_frame(StartPage), bg='brown', fg='white', font=('helvetica', 9, 'bold'))
        canvas2.create_window(800, 950, window=last_page_button)


# Driver Code
app = tkinterApp()
app.mainloop()
