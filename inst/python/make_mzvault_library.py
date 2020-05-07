import sqlite3
import pandas as pd
import numpy as np
import re
from rdkit import Chem
from rdkit.Chem import PandasTools
from io import StringIO
import datetime
import struct
import os
import sys
from tqdm import tqdm

def read_cfm(spec):
    """ Read CFM spectra into a dictionary"""
    if not isinstance(spec, bytes):
        return(None)
        #sys.exit('Spectrum object is not of type bytes')
    if pd.isna(spec):
        return(None)
    if len(spec) == 0:
        return(None)
    d = pd.read_table(StringIO(spec.decode('utf-8')), header = None)
    xx = {}
    for i in [0,1,2]:
        if i < 2:
            idx = range(np.where(d.values == 'energy'+str(i))[0][0]+1,np.where(d.values == 'energy'+str(i+1))[0][0])
        else:
            idx = range(np.where(d.values == 'energy'+str(i))[0][0]+1,d.shape[0])
            idx = range(idx[0],np.where([xx.split(' ')[0] == '0' for xx in d[0].values])[0][0])
        s = []
        for a in d.loc[idx][0].values:
            s.append({k:v for k,v in zip(['mz','intensity'],a.split(' ')[:2])})
        xx['energy'+str(i)] = pd.DataFrame(s).to_dict()
    return(xx)

def PackBinaryBuffer(x):
    x_b = bytearray(8*len(x))
    offset = 0
    for v in x:
        struct.pack_into('<d', x_b, offset, float(v))
        offset = offset + 8 
    buff = memoryview(x_b)
    return(buff)

def clean_spec(x,n_peaks,massAccuracy):
    if not isinstance(x, dict):
        sys.exit('Input should be a dictionary')
    if not [k for k in x] == ['mz','intensity']:
        sys.exit('Dictionary should have keys mz and intensity')
    df = pd.DataFrame(data = x, dtype = float)
    df = df.sort_values(by = "intensity", ascending = False)
    df = df.head(n = n_peaks)
    df = df.sort_values(by = "mz")
    df['accuracy'] = (df.mz/1e6)*massAccuracy
    return(df)

def make_empty_mzvault(mzvault_path):
    """Create an empty mzVault library"""
    conn = sqlite3.connect(mzvault_path)
    curs = conn.cursor()
    curs.execute("""CREATE TABLE CompoundTable([CompoundId] INTEGER PRIMARY KEY, 
                        [Formula] TEXT, [Name] TEXT, [Synonyms] BLOB_TEXT, [Tag] TEXT, 
                        [Sequence] TEXT, [CASId] TEXT, [ChemSpiderId] TEXT, [HMDBId] TEXT, 
                        [KEGGId] TEXT, [PubChemId] TEXT, [Structure] BLOB_TEXT, 
                        [mzCloudId] INTEGER, [CompoundClass] TEXT)""")
    curs.execute("""CREATE TABLE HeaderTable(version INTEGER NOT NULL DEFAULT 0, 
                        [CreationDate] TEXT,[LastModifiedDate] TEXT,
                        [Description] TEXT,[Company] TEXT)""")
    curs.execute("""CREATE TABLE InChIKeyTable([CompoundId] INTEGER PRIMARY KEY, 
                        [InChIKey] TEXT, [InChIKey14] TEXT)""")
    curs.execute("""CREATE TABLE 'SpectrumTable' ([SpectrumId] INTEGER PRIMARY KEY,
                        [CompoundId] INTEGER REFERENCES [CompoundTable] ([CompoundId]),
                        [mzCloudURL] TEXT, [ScanFilter] TEXT, [RetentionTime] DOUBLE, 
                        [ScanNumber] INTEGER, [PrecursorMass] DOUBLE, [NeutralMass] DOUBLE,
                        [CollisionEnergy] TEXT, [Polarity] TEXT, [FragmentationMode] TEXT,
                        [IonizationMode] TEXT, [MassAnalyzer] TEXT, [InstrumentName] TEXT, 
                        [InstrumentOperator] TEXT, [RawFileURL] TEXT, [blobMass] BLOB, 
                        [blobIntensity] BLOB, [blobAccuracy] BLOB, [blobResolution] BLOB, 
                        [blobFlags] BLOB, [blobTopPeaks] BLOB, [Version] INTEGER,
                        [CreationDate] TEXT, [Curator] TEXT, PrecursorType TEXT, 
                        blobRecalMass BLOB, blobRecalIntensity BLOB, blobExplanations BLOB, 
                        MassAccuracy DOUBLE)""")
    curs.execute("""CREATE INDEX idx_inchikey ON InChIKeyTable(InChIKey14)""")
    curs.execute("""CREATE INDEX idx_spectrumid ON SpectrumTable(SpectrumId)""")
    headerData = pd.DataFrame(data = [[4,datetime.datetime.today().strftime('%Y-%m-%d'),
                                          datetime.datetime.today().strftime('%Y-%m-%d'),
                                          "CFM MS/MS Library","Duke University"]], 
                                  columns = ["version", "CreationDate", "LastModifiedDate", "Description", "Company"])
    headerData.to_sql(name = "HeaderTable", con = conn, if_exists = "append", index = False)
    conn.commit()
    return(conn)

def make_mzvault(input_db, mol_tab, mol_idx, identifier_tab, cfmpred_pos, cfmpred_neg, ppm = 5):
    """ 
    Create an mz vault library from a database file.
    
    Parameters:
    -----------
    input_db: Full file path to the database containing cfm spectra 
    mol_tab: Name of the table containing molecules as SMILES strings 
    mol_idx: Name of the variable used for indexing the molecules and predicted
       spectra tables. 
    identifier_tab: Name of table in input database containing identifiers. 
      One of the identifiers must be an InChIKey. Others recognized (optional) 
      values include CASRN, Name, or InChI. Unrecognized values are concatenated 
      into a single string and added to a synonyms column. 
    cfmpred_pos: Name of table containing the positive ion predicted spectra. 
    cfmpred_neg: Name of table containing the negative ion predicted spectra. 
    ppm: Mass error for constructing the the mass accuracy binary buffer. 
    """
    conn = sqlite3.connect(input_db)
    df = pd.read_sql('SELECT * FROM '+mol_tab+';', conn,index_col=mol_idx)
    if identifier_tab is not None:
        identifiers = pd.read_sql('SELECT * FROM '+identifier_tab+';', conn)
    else:
        identifiers = None 
    cfm_neg = pd.read_sql('SELECT * FROM '+cfmpred_neg+';', con=conn,index_col=mol_idx)
    cfm_pos = pd.read_sql('SELECT * FROM '+cfmpred_pos+';', con=conn,index_col=mol_idx)
    cfm = pd.merge(cfm_neg, cfm_pos, on = mol_idx, how = 'outer',suffixes = ('_neg', '_pos'))
    df = pd.merge(df, cfm, on = mol_idx, how = 'left')
    smiles_col = [x for x in df.columns.values if re.search('smiles', x, re.I)][0]
    PandasTools.AddMoleculeColumnToFrame(df, smiles_col, 'ROMol',False)
    output_db = re.sub('.db', '_mzvault.db', input_db)
    conn = make_empty_mzvault(output_db)
    curs = conn.cursor()
    accession = 0
    pbar = tqdm(total = len(df))
    for x,mol,spec_neg,spec_pos in zip(df.index.values,df['ROMol'],df['spectrum_neg'],df['spectrum_pos']):
        pbar.update(1)
        if mol is None: continue
        neutralMass = Chem.rdMolDescriptors.CalcExactMolWt(mol)
        precmz = {'+':neutralMass + 1.00727647, '-':neutralMass - 1.00727647}
        precType = {'+':'[M+H]+', '-':'[M-H]-'}
        ionMode = {'+':'P', '-':'N'}
        formula = Chem.rdMolDescriptors.CalcMolFormula(mol)   
        molblock = Chem.rdmolfiles.MolToMolBlock(mol)
        inchikey = Chem.MolToInchiKey(mol)
        if identifiers is not None:
            name_col = [x for x in identifiers.columns.values if re.search('name',x,re.I)][0]
            inchikey_col = [x for x in identifiers.columns.values if re.search('^inchikey',x,re.I)][0]
            inchi_col = [x for x in identifiers.columns.values if re.search('inchi(?!key)',x,re.I)][0]
            casrn_col = [x for x in identifiers.columns.values if re.search('casrn',x,re.I)][0]
            synonym_col = [x for x in identifiers.columns.values if x not in [name_col,inchikey_col,casrn_col,inchi_col]]
            try:
                name = identifiers[identifiers[inchikey_col] == inchikey][name_col]
                name = ','.join([x for x in name.values.tolist()])
            except IndexError:
                name = " "
        
            try:
                casrn = identifiers[identifiers[inchikey_col] == inchikey][casrn_col]
                casrn = ','.join([x for x in casrn.values.tolist()])
            except IndexError:
                casrn = " "
        
            try:
                synonyms = identifiers[identifiers[inchikey_col] == inchikey][synonym_col].transpose().values.tolist()
                synonyms = ','.join([y for x in synonyms for y in x])
            except IndexError:
                synonyms = " "
        else:
            name = " "
            casrn = " "
            synonyms = " "
            
        spec = {'+':read_cfm(spec_pos), '-':read_cfm(spec_neg)}
        spec = {k:spec[k] for k in spec if spec[k] is not None}
        for pol in spec:
            for energy in spec[pol]:
                # Format scan header
                scanHeader = " ".join(filter(None,["CFM", str(accession), 'CFM', 'prediction', 'EE', pol, "CE=", re.sub("energy","", energy)]))
                # Get spectrum and keep only top 200 most intense peaks
                spectrum = clean_spec(spec[pol][energy], 200, ppm)
                
                # Create binary buffers for mz and intensity for SQLite BLOB
                mzBuf = PackBinaryBuffer(spectrum.mz)  
                intensityBuf = PackBinaryBuffer(spectrum.intensity)
                accuracyBuf = PackBinaryBuffer(spectrum.accuracy)
                
                # Write Compound Data to database
                curs.execute("SELECT EXISTS(SELECT 1 FROM InChIKeyTable WHERE InChIKey14= ? LIMIT 1)", (inchikey[:14],))
                cmpd_exists = curs.fetchall()[0][0]>0
                if not cmpd_exists:
                    # Set compound ID
                    cmpdID = pd.read_sql('SELECT MAX(CompoundId) FROM CompoundTable', conn).values[0][0]
                    if cmpdID is None:
                        cmpdID = 1
                    else:
                        cmpdID = cmpdID + 1
                    # Append compound data to database
                    cmpdData = pd.DataFrame(data = [[cmpdID,formula,name,synonyms,casrn,None,molblock]], 
                                  columns = ["CompoundId", "Formula", "Name", "Synonyms", "CASId", "PubChemId", "Structure"])
                    cmpdData.to_sql(name = "CompoundTable", con = conn, if_exists = "append", index = False)
                    conn.commit()
                    # Append InChIKey data to database
                    inchikeyData = pd.DataFrame(data = [[cmpdID,inchikey,inchikey[:14]]],
                                    columns = ["CompoundId", "InChIKey", "InChIKey14"])
                    inchikeyData.to_sql(name = "InChIKeyTable", con = conn, if_exists = "append", index = False)
                    conn.commit()
                else:
                    # Get compound ID from InChIKey table
                    curs.execute("SELECT CompoundId FROM InChIKeyTable WHERE InChIKey14 = ?", (inchikey[:14],))
                    cmpdID = curs.fetchall()[0][0]
                # Write Spectrum Data to database
                spectrumData = pd.DataFrame(data = [[accession,cmpdID,precmz[pol],neutralMass,re.sub("energy","", energy),pol,'PREDICTION',
                                         precType[pol],ionMode[pol],'PREDICTION','CFM-ID',mzBuf,
                                         intensityBuf,accuracyBuf,None,None,
                                         None,ppm,scanHeader]],
                                    columns = ["SpectrumId","CompoundId","PrecursorMass",
                                               "NeutralMass","CollisionEnergy","Polarity",
                                               "FragmentationMode","PrecursorType","IonizationMode",
                                               "MassAnalyzer","InstrumentName","blobMass",
                                               "blobIntensity","blobAccuracy","blobRecalMass",
                                               "blobRecalIntensity","blobExplanations",
                                               "MassAccuracy","ScanFilter"])
                spectrumData.to_sql(name = "SpectrumTable", con = conn, if_exists = "append", index = False)
                conn.commit()
                accession+=1
    pbar.close()
    return(output_db)
