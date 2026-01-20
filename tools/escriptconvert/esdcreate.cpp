
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* https://github.com/LutzGross/esys-escript.github.io
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#include <iostream>
#include <fstream>

using namespace std;

int main(int argc, char** argv)
{
    string dsName;

    if (argc > 1) {
        cerr << "esdcreate - escript datafile creator version 1.0" << endl;
        cerr << "This program takes no arguments." << endl;
        return -1;
    }

    for (;;) {
        cout << "Please enter a name for this dataset (e.g. simulation): ";
        cin >> dsName;
        dsName.append(".esd");
        ifstream f;
        f.open(dsName.c_str(), ifstream::in);
        // File exists
        if (f.good()) {
            cout << "A file by the name " << dsName << " already exists! Overwrite? (y/n) ";
            f.close();
            char a = ' ';
            while (a != 'y' && a != 'n')
                cin >> a;
            if (a == 'y')
                break;

        } else
            break;
    }

    ofstream esdFile;
    esdFile.open(dsName.c_str(), ifstream::out);
    esdFile << "#escript datafile V1.0" << endl;

    int numTS = 0;
    while (numTS <= 0) {
        cout << "Number of timesteps: ";
        cin >> numTS;
        if (!cin.good()) {
            cin.clear();
            cin.ignore(100, '\n');
        }
        if (numTS <= 0)
            cout << "Please enter a value > 0!" << endl;
    }
    esdFile << "T=" << numTS << endl;

    string meshName;
    int numParts = 0;
    for (;;) {
        cout << "Filename of the mesh (without .nc): ";
        cin >> meshName;
        string fname = meshName;
        fname.append(".nc");
        ifstream f;
        f.open(fname.c_str(), ifstream::in);
        if (f.good()) {
            f.close();
            numParts = 1;
        } else {
            char partnum[10];
            for (;; numParts++) {
                snprintf(partnum, 10, ".%04d", numParts);
                string pname = fname;
                pname.append(partnum);
                f.open(pname.c_str(), ifstream::in);
                if (f.good())
                    f.close();
                else
                    break;
            }
        }
        if (numParts == 0)
            cout << "The mesh was not found! Please try again." << endl;
        else
            break;
    }
    esdFile << "M=" << meshName << endl;
    esdFile << "N=" << numParts << endl;

    cout << endl;
    cout << "For each variable you would like to add please enter the filename" << endl;
    cout << "of the variable (without .nc) followed by a colon and a name for" << endl;
    cout << "the variable. For example: temp:temperature" << endl;
    if (numTS > 1) {
        cout << "Since you indicated that there are multiple timesteps, you have to put" << endl;
        cout << "a format string at the right place of the filename which will be replaced" << endl;
        cout << "by the actual timestep when loading the files." << endl;
        cout << "Examples:" << endl;
        cout << "  temp.%d:temperature would be resolved to temp.0.nc, temp.1.nc etc." << endl;
        cout << "  temp.%04d:temperature would be resolved to temp.0000.nc, temp.0001.nc etc." << endl;
    }
    cout << endl;
    cout << "When you are finished please enter: done" << endl;

    for (;;) {
        string varString;
        cout << "Filename and variable name: ";
        cin >> varString;
        if (varString == string("done"))
            break;
        if (varString.find(':') == varString.npos) {
            cout << "Your input is invalid! Please enter filename:varname or done." << endl;
        } else
            esdFile << "V=" << varString << endl;
    }

    if (numTS > 1) {
        int tsIncrement = 0;
        cout << endl;
        cout << "Please supply the timestep increment to find your files." << endl;
        cout << "For example, if your files are named temp.0, temp.10, temp.20, etc." << endl;
        cout << "the increment is 10, whereas temp.0, temp.1, temp.2 have an increment of 1." << endl;
        while (tsIncrement <= 0) {
            cout << "Timestep increment for filenames: ";
            cin >> tsIncrement;
            if (!cin.good()) {
                cin.clear();
                cin.ignore(100, '\n');
            }
            if (tsIncrement <= 0)
                cout << "Please enter a value > 0!" << endl;
        }
        esdFile << "DT=" << tsIncrement << endl;
    }

    esdFile.close();

    cout << endl;
    cout << "Your escript datafile has been written to " << dsName << "." << endl;
    cout << "You can now open this file from within VisIt to visualise your data." << endl;

    return 0;
}

