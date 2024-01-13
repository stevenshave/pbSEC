MSCAT_PATH="""\"c:\\Program Files\\ProteoWizard\\ProteoWizard 3.0.10462\\mscat.exe\" """
PICKER_PATH="pickerppm.exe"


import tkinter as tk
import ntpath
import subprocess
import time
import tkinter.messagebox as messagebox
from tkinter.filedialog import askdirectory, askopenfilename
import csv
import os
import glob
import os



def ltqmassesmain():
    tk.Tk().withdraw()
    pools={}
    messagebox.showinfo("MassFinder", "Please select the pool definition file.")
    pooldef=askopenfilename()
    resultsfile="results-"+ntpath.basename(pooldef)
        
    with open(pooldef, "r") as csvFile:
        reader=csv.reader(csvFile)
        for row in reader:
            if row[0] in pools:
                pools[row[0]].append(row[1])
            else:
                pools[row[0]]=[]
                pools[row[0]].append(row[1])
    
    print("Found definitions + masses for", len(pools), "pools")

    messagebox.showinfo("MassFinder", "Please specify the data directory containing mass spec files")
    datadir = askdirectory()
    
    os.chdir(datadir)
    try:
        os.remove(resultsfile)
    except OSError:
        pass

    f = open(resultsfile, 'w')
    f.write('Well File Name,Mass,IonCount,IonCount-Minute3, IonCount-Minute4,...\n')

    for pool in pools:
        dirlisting = []
        dirlisting = glob.glob("*"+pool+"*.RAW")
        print(pool, "-----DIRLISTING", dirlisting)
        if(len(dirlisting)!=1):
            continue
        massesstring=""
        for mass in pools[pool]:
            massesstring=massesstring+str(mass)+" "
        commandstring = MSCAT_PATH +  dirlisting[0] + "| "+PICKER_PATH+ " " + resultsfile + " " +dirlisting[0]+" 10.0 "+massesstring
        print(commandstring)
        os.system(commandstring)
    messagebox.showinfo("MassFinder completed", "Results written to "+datadir+"/"+resultsfile)
    return


if __name__ == '__main__':
        try:
            ltqmassesmain()
            input()
        except:
            import sys
            print(sys.exc_info()[0])
            import traceback
            print(traceback.format_exc())
            print("Press enter to exit ...")
            input()
