#include "Data.h"

string myText;

void Data::getParametersTest(ifstream& MyReadFile) {
    int nr = 0;
    int coma = 0;
    string number = "";
    try {
        while (getline(MyReadFile, myText) && myText[0] >= 'A' && myText[0] <= 'z') {
            for (unsigned int i = 0; i < myText.size(); i++)
                if (myText[i] >= '0' && myText[i] <= '9')
                    params[nr] = params[nr] * 10 + myText[i] - 48;
            nr++;
        }
        nr = 0;
        d_nodes = new double* [params[8]];
        for (int i = 0; i < params[8]; i++) {
            d_nodes[i] = new double[2];
        }
        d_elems = new int* [params[9]];
        for (int i = 0; i < params[9]; i++) {
            d_elems[i] = new int[4];
        }
        for (int i = 0; i < params[8]; i++) {
            for (int j = 0; j < 2; j++) {
                d_nodes[i][j] = 0.;
            }
        }
        for (int i = 0; i < params[9]; i++) {
            for (int j = 0; j < 4; j++) {
                d_elems[i][j] = 0.;
            }
        }
        nr = 0;

        while (getline(MyReadFile, myText) && myText[1] != 'E') {
            for (int i = 0; i < myText.size(); i++) {
                if ((coma == 1 || coma == 2) && (myText[i] >= '0' && myText[i] <= '9' || myText[i] <= '.')) {
                    number += myText[i];
                }

                if (myText[i] == ',') {
                    if (coma == 1) {
                        d_nodes[nr][0] = stod(number);
                        number = "";
                    }
                    coma++;
                }
                if (i == myText.size() - 1) {
                    d_nodes[nr][1] = stod(number);
                    number = "";
                }
            }
            nr++;
            coma = 0;
        }

        nr = 0;
        number = "";
        coma = 0;

        while (getline(MyReadFile, myText) && myText[1] != 'B') {
            for (int i = 0; i < myText.size(); i++) {
                if ((coma != 0) && (myText[i] >= '0' && myText[i] <= '9')) {
                    number += myText[i];
                }

                if (myText[i] == ',') {
                    if (coma != 0) {
                        d_elems[nr][coma - 1] = stoi(number);
                        number = "";
                    }
                    coma++;
                }
                if (i == myText.size() - 1) {
                    d_elems[nr][3] = stoi(number);
                    number = "";
                }
            }
            nr++;
            coma = 0;
        }

        number = "";
        coma = 0;

        while (getline(MyReadFile, myText)) {
            for (int i = 0; i < myText.size(); i++) {
                if (myText[i] >= '0' && myText[i] <= '9') {
                    number += myText[i];
                }

                if (myText[i] == ',') {
                    number = "";
                    coma++;
                    bc_counter++;

                }
                if (i == myText.size() - 1) {
                    bc_counter++;
                    number = "";
                }

            }
            coma = 0;
            //cout << "from while loop number counter:" << bc_counter << '\n';
            bcs = new int[bc_counter];

            for (int i = 0; i < myText.size(); i++) {
                if (myText[i] >= '0' && myText[i] <= '9') {
                    number += myText[i];
                }

                if (myText[i] == ',') {
                    bcs[coma] = stoi(number);
                    //cout << "from for loop number:" << bcs[coma] << '\n';
                    number = "";
                    coma++;
                }
                if (i == myText.size() - 1) {
                    bcs[coma] = stoi(number);
                    //cout << "from for loop number:" << bcs[coma] << '\n';
                    number = "";
                }

            }
            coma = 0;
        }

        // Assign data values
        simTime = params[0];
        dTau = params[1];
        k = params[2];
        alpha = params[3];
        t_env = params[4];
        T0 = params[5];
        ro = params[6];
        c = params[7];
        nN = params[8];
        nE = params[9];
    }
    catch (...) {
        std::cout << "Data import failed." << std::endl;
    }
}

void Data::setParameters(Grid &G){

    G.nodes = new Node[G.nN];
    G.elements = new Element[G.nE];

    for (int i = 0; i < G.nN; i++) {
        G.nodes[i].x = d_nodes[i][0];
        G.nodes[i].y = d_nodes[i][1];
        G.nodes[i].t0 = T0;
    }
    for (int i = 0; i < G.nE; i++) {
        G.elements[i].ID[0] = d_elems[i][0];
        G.elements[i].ID[1] = d_elems[i][1];
        G.elements[i].ID[2] = d_elems[i][2];
        G.elements[i].ID[3] = d_elems[i][3];
    }
    for (int i = 0; i < bc_counter; i++) {
        G.nodes[bcs[i]-1].BC = 1;
    }
}
