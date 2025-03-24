# -*- coding: utf-8 -*-

class mineral_id():
    def __init__(self, path = 'blank', analytical_tool = 'EDS', output = 'recent'):
        from pandas import read_excel
        self.__path = path
        self.data = read_excel(io = path)
        self.periodic_table = read_excel(io="data_sheets/periodic_table.xlsx")
        self.working = 'blank'
        if analytical_tool in ['SEM - EDS', 'SEM-EDS', 'EDS', 'eds', 'SEM', 'sem', 'sem-eds']:
            self.formulas = read_excel(io="data_sheets/eds_formulas.xlsx")
        self.output = output

    def find_headers(self, frame):
        '''
        Produces a list of non element column names, which must be removed before
        many arithmentic manipulations are done. Called by to_mol_percent and min_id

        Returns
        -------
        string_cols: list
            A list of the non element column headers

        '''
        headers = list(frame.columns)
        elements = list(self.periodic_table['Chemical sym'])
        string_cols = []
        for header in headers:
            if header not in elements:
                string_cols.append(header)   
        return string_cols
        
    def remove_headers(self, df):
        '''
        Removes non element headers of dataframes

        Parameters
        ----------
        df : pandas DataFrame
            A dataframe containing element column names (e.g. 'H', 'Ag', 'Li', Hf', etc.),
            and non element column names (e.g. 'total', 'sample id', etc.)

        Returns
        -------
        frame : pandas DataFrame
            The input dataframe with all indexes and collumns in their origonal location, 
        except for all columns with names that are not elemental symbols (e.g. He, Xe, U, etc.)

        '''
        drop = self.find_headers(df)
        frame = df.copy()
        for i in drop:
            frame.pop(i)
        return frame
        
    
    def to_mol_percent(self):
        '''
        Renormalizes wt% chemistry to mol% chemistry

        Returns
        -------
        self.working: A dataframe containing only mol% elemental abundance, ordered as
        self.data is, i.e. with the same index order and the same columns in the same order.

        '''
        #removing non-element headers
        drop = self.find_headers(self.data)
        wrk = self.data.copy()
        for i in drop:
            wrk.pop(i)
        
        #the bulk of the method
        for column in wrk.columns:
            #defining the index (as an int) of self.__periodic_table that contains an element in wrk.columns
            e_index = self.periodic_table.index[self.periodic_table['Chemical sym']==column].tolist()
            for l in e_index:
                e_index = l
            
            #Multiplying each column of wrk by the molar mass of its column header element
            e_mass = self.periodic_table.at[e_index, 'Atomic mass (g/mol)']
            wrk[column] = wrk[column]/e_mass
            
        #normalizing the molar makeup each sample to a percent
        total = wrk.sum(axis=1)/1
        
        #returning a dataframe containing only mol% elemental abundance
        self.working = wrk.divide(total,axis =0)
    
    
    def ident(self):
        '''
        Perfroms the arithmetic discrimination of sample identification. Measured mol% chemistry
        of each sample is subtracted from ideal mol% chemistry. The absolute value of those differences
        are summed (and multiplied by 50) producing a mismatch percentage as the discrimination
        value. For each unknown, the ideal formula with the lowest mismatch percentage is presented
        as the correct mineral match, along with its missmatch percentage. Below 10% is usually safe
        to call a correct match. Please note that non-ideal chemistry (i.e. phases exhibiting solid
        solution) will have a higher mismatch percentage, perhaps higher than the 10% threshold.

        '''
        unknown = self.working
        targets = self.remove_headers(self.formulas)
        minerals = []
        miss = []
        for row in list(unknown.index):
            line = unknown.iloc[row]
            difference = abs(targets.subtract(line))
            total_diff = difference.sum(axis = 1)*50
            miss.append(total_diff.min())
            best_fit = total_diff.idxmin()
            identity = self.formulas['name'].iloc[best_fit]
            minerals.append(identity)
        self.data['mineral estimate'] = minerals
        self.data['mismatch percent'] = miss

    def identify(self):
        self.to_mol_percent()
        self.ident()
        self.data.to_csv(path_or_buf=self.output+'.csv', index=False, na_rep=0)

##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################        
    
class mineral_formula():
    def __init__(self, name, formula):
        from pandas import read_excel
        self.formula = formula
        self.form = formula
        self.name = name
        self.periodic_table = read_excel(io="data_sheets/periodic_table.xlsx")
        self.new_line = 'blank'
        
    def check_min(self):
        from pandas import read_excel
        ref_data = read_excel(io ="data_sheets/trunk_formulas.xlsx")
        name = False
        formula = False
        match = False
        line = 'blank'
        if self.name in list(ref_data['name']):
            name = True
        if self.name in list(ref_data['formula']):
            formula = True
        if name == True:
            match = True
            i = ref_data.index[ref_data['name']==self.name].tolist() # a list of 
            for l in i:
                i = l  
            line = ref_data.loc[i]
        elif formula == True:
            match = True
            i = ref_data.index[ref_data['name']==self.formula].tolist() # a list of 
            for l in i:
                i = l  
            line = ref_data.loc[i]
        return match, line
        
    def to_chemistry(self):    
        '''
        DEV NOTE - This docstring has been wholesale adopted from a previous standalone
        function and may not apply perfectly. PLEASE CHANGE
        
        This function is built to turn chemical formulas into ideal weight percent
        chemistry for the formula. *CURRENTLY* the formula must be input as elements
        and stoichiometry (e.g. H2O,CaCO3, etc.) only. It cannot manage any forms 
        of brackets, for repeated chemical units (e.g. 
        K2(Mg6)(Si6Al2)O20(OH)4), geochemical "solid solution" style 
        replacement (e.g. [Fe,Mn][Ta,Nb]2O6), or unit valence (e.g.
        Fe{2+}Fe{3+}2O4). It also cannot yet manage repeated elements in a formula
        (e.g. FeOOH). All of these should be implemented later in order to make a
        robust, reusable function.
    
        Parameters
        ----------
        formula : string 
            The formula in question. Elements must be in elemental symbol format 
            (e.g. H, He, Al, Fe) with any stoichiometric variation following as 
            whole integers (e.g. H2, CaCO4), fractional stoichiometry cannot be 
            understood by this function, so any fractions must be multiplied to
            whole integers first.
            
        name : string
            The name of the formula. (e.g. if formula is H2O name should be 'water')

        Returns
        -------
        df : pandas dataframe
            A single row of a pandas dataframe. The index of the produced row is 0.
            The collumns are 'name', 'formula', and then the elements from the 
            input formula in the order they appear. For example, if the input 
            formula was CaCO4, df.columns would return ['name', 'formula', 'Ca', 
                                                    'C', 'O']. Values in the element columns are 'weight fractions' of the
            element in the entire formula. True weight percent would be these units 
            * 100.
        
        '''
        import pandas as pd
        #setting the first two collumns of the product dataframe
        
        data = [[self.name, self.formula]]
        df = pd.DataFrame(data = data, columns = ['name', 'formula'])
    
        #Defining needed lists
        element_caps = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                        'M', 'N', 'O', 'P', 'R', 'S', 'T', 'U', 'V', 'W', 'X',
                        'Y', 'Z']
        numbers = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']
   
        for letter in element_caps:
            divider = '-'+letter
            self.formula = self.formula.replace(letter, divider)
            self.formula = self.formula.lstrip('-')
        #turning the formula into a list of single element chunks
        self.formula = self.formula.split('-')
        forms = []
    #separating single element chunks into element names and multipliers
        for element_set in self.formula:
            counter = 0
            for number in numbers:
                if number in element_set:
                    divider='-'+number
                    element_set = element_set.replace(number, divider)
                    forms.append(element_set)
                    break
                else:
                    counter +=1
            if counter == 10:
                element_set = element_set+'-1'
                forms.append(element_set)
    
        multipliers = []
        headers = []
        #moving the elements into a spreadsheet
        for element_set in forms:
            element, value = element_set.split('-')
            multipliers.append(value)
            headers.append(element)
        multipliers = [multipliers]
        wrk = pd.DataFrame(data = multipliers, columns = headers)

        #multiplies the working columns by their molar mass
        for c in wrk.columns:
            i = self.periodic_table.index[self.periodic_table['Chemical sym']==c].tolist() # a list of 
            for l in i:
                i = l            
            k=float(wrk.at[0,c])
            j=self.periodic_table.at[i, 'Atomic mass (g/mol)']
            j = j * k
            wrk[c] = j
    
        #defining a divisor
        div = wrk.sum(axis = 1)
    
        #Finalizing the arithmetic expression wt% = (atomic units * atomic mass)/ sum(atomic units*atomic masses)
        wrk = wrk.divide(div, axis = 0)
    
        #combining the working table with the leading collumns of the final dataframe
        df = pd.concat([df, wrk], axis=1)
        self.new_line = df
    
    def renorm(self, excl):
        '''
        This removes elements from the ideal mineral chemistry spreadsheet, and 
        it renormalizes each formula to 100% without the removed elements. It
        will be called when creating ideal mineral chemistry for different instruments
        that cannot detect certain elements (e.g. for SEM EDS; Li, H, He, etc.)

        Parameters
        ----------
        excl : type: list
            The elements to be removed from the new formula. Must be a list.

        Returns
        -------
        renormalized : pandas DataFram
            The formula with elements from excl removed, and the remaining elements
            renormalized to a 1.0 total

        '''
        from pandas import DataFrame
        from pandas import concat
        renormalized = DataFrame(data = [[self.name, self.form]], columns = ['name', 'formula'])
        renormalized = concat([renormalized, self.new_line], axis = 1)
        elms = []
        for element in excl:
            if element in self.new_line.columns:
                elms += element
        elements = self.new_line.drop(labels = elms, axis = 1)
        renorm_factor = DataFrame(data = [0], columns = [1])
        
        for element in elms:
            temp = renormalized.pop(element)
            renorm_factor[1]+= temp
        for column in elements:
            renormalized[column] = renormalized[column]/(1-renorm_factor[1])
        return renormalized
        
    
    def find_headers(self, frame):
        '''
        Produces a list of non element column names, which must be removed before
        many arithmentic manipulations are done. Called by to_mol_percent and min_id

        Returns
        -------
        string_cols: list
            A list of the non element column headers

        '''
        headers = list(frame.columns)
        elements = list(self.periodic_table['Chemical sym'])
        string_cols = []
        for header in headers:
            if header not in elements:
                string_cols.append(header)   
        return string_cols
        
    def remove_headers(self, df):
        '''
        Removes non element headers of dataframes

        Parameters
        ----------
        df : pandas DataFrame
            A dataframe containing element column names (e.g. 'H', 'Ag', 'Li', Hf', etc.),
            and non element column names (e.g. 'total', 'sample id', etc.)

        Returns
        -------
        frame : pandas DataFrame
            The input dataframe with all indexes and collumns in their origonal location, 
        except for all columns with names that are not elemental symbols (e.g. He, Xe, U, etc.)

        '''
        drop = self.find_headers(df)
        frame = df.copy()
        for i in drop:
            frame.pop(i)
        return frame
        
    
    def to_mol_percent(self):
        '''
        Renormalizes wt% chemistry to mol% chemistry

        Returns
        -------
        self.working: A dataframe containing only mol% elemental abundance, ordered as
        self.data is, i.e. with the same index order and the same columns in the same order.

        '''
        #removing non-element headers
        drop = self.find_headers(self.new_line)
        wrk = self.new_line.copy()
        for i in drop:
            wrk.pop(i)
        
        #the bulk of the method
        for column in wrk.columns:
            #defining the index (as an int) of self.__periodic_table that contains an element in wrk.columns
            e_index = self.periodic_table.index[self.periodic_table['Chemical sym']==column].tolist()
            for l in e_index:
                e_index = l
            
            #Multiplying each column of wrk by the molar mass of its column header element
            e_mass = self.periodic_table.at[e_index, 'Atomic mass (g/mol)']
            wrk[column] = wrk[column]/e_mass
            
        #normalizing the molar makeup each sample to a percent
        total = wrk.sum(axis=1)/1
        
        #returning a dataframe containing only mol% elemental abundance
        self.new_line = wrk.divide(total,axis =0)

    
    def append_spreadsheets(self):
        from pandas import read_excel
        from pandas import concat
        self.to_chemistry()
        self.to_mol_percent()
        
        #apending the trunk datasheet
        trunk = self.renorm([])
        trunk_data = read_excel(io ="data_sheets/trunk_formulas.xlsx")
        trunk_data = concat([trunk_data, trunk], ignore_index = True)
        trunk_data = trunk_data.fillna(0)
        trunk_data.to_excel(excel_writer="data_sheets/trunk_formulas.xlsx", index=False, na_rep=0, engine="openpyxl")

        #apending the EDS datasheet
        eds = self.renorm(excl = ['H', 'He', 'Li', 'Be', 'B', 'Ne', 'Ar', 'Kr', 'Xe', 'Rn', 'Og', 'Tc'])
        eds_data = read_excel(io ="data_sheets/eds_formulas.xlsx")
        eds_data = concat([eds_data, eds], ignore_index = True)
        eds_data = eds_data.fillna(0)
        eds_data.to_excel(excel_writer="data_sheets/eds_formulas.xlsx", index=False, na_rep=0, engine="openpyxl")
        
        
        
        
def example():
    """
    An example of how to add a mineral formula to the mineral formula spreadsheet
    """
    x = 'native silver'
    y = 'Ag'
    run = mineral_formula(name = x, formula = y)
    run.append_spreadsheets()
    
    
if __name__ == "__main__":
    run = mineral_id(path = "data_sheets\example_data.xlsx", analytical_tool = "SEM-EDS", output = 'Output/example.xlsx')
    run.identify()