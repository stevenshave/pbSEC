# Introduction
Plate based size exclusion chromatography (pbSEC) is a primary screening method used within the Auer lab.  Analysis of pools of compounds by hand, manually looking at spectra is time consuming an error prone.  For these reasons, an automated approach was required for analysis of this mass spectrometery data.  The program 'MassFinder' fulfils this purpose although the pipeline will likely need heavy modification for use on different mass spectrometer instruments.  It is however written to be as general as possible, relying upon the ProteoWizard package for reading a variety of mass spec formats.

# Requirements
The software pipeline consists of two programs
- massfinder.py
    - the main control/interface program written in python.
- pickerppm
    - a fast C++ text parser for proteowizard output which is called from the massfinder.py script.

Jointly, these two programs assume the computer on which it is running also has the following software installed:
- Python3 - https://www.python.org/ (>=3.6)
- Proteowizard - http://proteowizard.sourceforge.net/
- Visual C++ Redistributable for Visual Studio 2015 (x86 version) - https://www.microsoft.com/en-us/download/details.aspx?id=48145

With clang, pickerppm may be compiled for the target system with the command:

```bash
clang++ -o pickerppm pickerppm.cpp -O3
```

# Configuration

The MassFinder pipeline requires some configuration before being run by setting some paths in the massfinder.py script. They should be absolute paths and be surrounded by escaped quotes as shown in the defaults below. Newer versions of Python supporting f-strings would alieviate this clunky requirement, however, we required this to run on locked down instrumentation equipment and so had to compromise. The following paths must be set in the massfinder.py file:

- MSCAT_PATH
    - Points at the mscat.exe application which comes with the install of ProteoWizzard
    - By default """\"c:\\Program Files\\ProteoWizard\\ProteoWizard 3.0.10462\\mscat.exe\" """
- PICKER_PATH
    - Points at the pickerppm program 
    - By default """\"c:\\Program Files\\pickerppm.exe\""""
        - (Assuming compilation for Windows)
       
# Usage
Upon execution of massfinder.py, the user will be asked to specify a platemap defining which masses are to be looked for in each well. An extract from a platemap CSV file (comma separated, no header) follows:

|      |        |
|------|--------|
|pool01|270.0608|
|pool01|280.1121|
|pool01|323.1874|
|pool01|335.1543|
|pool02|275.0608|
|pool02|283.1121|
|pool02|281.0186|
|pool02|293.1227|
|pool02|336.1217|

Whilst the above names pools pool01 and pool02, any unique identifier may be used which matches mass spec output for that pool/well.

With the platemap and data directories identified, the MassFinder program will begin calling mscat.exe to decode vendor mass spec data files and pass output to Picker.exe.  This is a long process and best left running overnight.  The time consuming part of this is conversion from the vendor’s proprietary mass spec data format to something easily readable.  This conversion is achieved using the ProteoWizzard tool mscat.exe, and as such, there is nothing that we can do to speed the process up.

Results for each well (ion counts for each mass) are written to a results CSV file taking the form “results-[platemapname].csv”, where [platemapname] is the name of the platemap.  This file will be written to the data directory where the mass spec data is located. 

Results may then be inspected with excel, or similar spreadsheet analysis software, sorting each pool by ion counts, or the overall plate by ion count to see the priority compounds.  Once compounds are prioritised by sorting by overall ion count, a manual check must be undertaken, looking at the columns following “IonCount”.  These additional columns give the ion count per minute (ignoring, as the total count does, the first 2 minutes of the run).  If a compound comes out consistently over the course of the run, it should be deprioritised as it is not a true hit – simply going through the column on its own and not eluting with the protein (indicating binding).
 