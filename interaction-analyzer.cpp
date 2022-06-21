/* 
 * File:   interaction-analyzer.cpp
 * Author: Martin Smiesko, University of Basel, All rights reserved
 * 
 * Dissemination and reuse only with explicit permission from the author
 *
 * Created on February 7, 2020, 1:33 PM
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "definitions.h"
#include "points_vectors.h"
#include "structures.h"

// main routine to read the structural data in the PDB format
void    readProteinPDB(char * filename, Residue * &lastResidue, Molecule * &lastMolecule, Water * &lastWater, RawAtom * &lastRawAtom) {
	int             atomCounter = 0;
	int     	residueCounter = 0;
        int             waterCounter = 0;
        int             moleculeCounter = 0;
	int             resNum;
	int             lastNum  = -1;
	int             atomNum;
        unsigned char   N;
	char		letter = '-';	// amino acid one letter code + extended ones
	char		lineBuffer[255];
	char		resName[5];
	char		lastName[5];
	char		atomType[5];
        char		ts[10];         // temporary string
        char            element[3];
        char            chain;
        char            lastChain;
        char            altLoc;
        double          occupancy;
        double          Bfactor;
	Point		pP;		// pP holding protein atom's coordinates
        int             charge;
        Residue *       curResidue = NULL;
        Molecule *      curMolecule = NULL;
        Water *         curWater = NULL;
        char            flag;
        
        ProteinAtom *   pa;
        LigandAtom *    la;
        RawAtom *       ra;
        RawAtom *       ra2;
        
        lastResidue = NULL; 
        lastMolecule = NULL;
        lastWater = NULL;
        
        
        FILE *	input = fopen(filename, "rt");
	
	// open file for reading
	if (input == NULL) { printf("Error! Problem opening file %s.\n", filename); exit(1); }

        // clear resName
	strcpy(lastName, "");
        lastChain = '-';
	
        // START OF READING COORDINATES
        // READ ATOM & HETATM
        while (fgets(lineBuffer, sizeof(lineBuffer), input) != NULL  &&  strstr(lineBuffer, "ENDMDL") != lineBuffer)
	{
                // STRUCTURAL DATA LINES
                     if(strstr(lineBuffer, "ATOM") == lineBuffer) flag = ATOM_LINE;
                else if(strstr(lineBuffer, "HETATM") == lineBuffer) flag = HETATM_LINE;
                else flag = NO;
                
                if(flag != NO) {
                        // atom number
                        if(flag == ATOM_LINE) { strncpy(ts, lineBuffer+4, 7); ts[7] = 0; }
                        else if(flag == HETATM_LINE) { strncpy(ts, lineBuffer+6, 5); ts[5] = 0; }
                        
                        if(sscanf(ts, "%d", &atomNum) != 1) { printf("Error1 (atom num) reading PDB line : %s\n", lineBuffer); exit(1); };
                        if(atomNum < 0) { printf("Error1 atom number smaller than zero at line : %s\n", lineBuffer); exit(1); };
                        
                        
                        // OVERRIDE TO READ HETATM LIKE ATOM
                        flag = ATOM_LINE;
                        
                        // atom type
                        strncpy(atomType, lineBuffer+12, 4);
                        atomType[4] = 0;

                        // alternative orientation
                        altLoc = lineBuffer[16];
                        
                        // res name
                        strncpy(resName, lineBuffer+17, 3);
                        resName[3] = 0;
                        
                        // res num
                        strncpy(ts, lineBuffer+22, 4);
                        ts[4] = 0;
                        if(sscanf(ts, "%d", &resNum) != 1) { 
                            if(strcmp(resName, "T3P") == 0  ||  strcmp(resName, "WAT") == 0  ||  strcmp(resName, "HOH") == 0  ||  strcmp(resName, "DOD") == 0) {
                                resNum = 0;
                                
                                if(ts[0] >= 'A') resNum = resNum + (ts[0] - 'A' + 10) * 1000;  
                                else resNum = resNum + (ts[0] - '0') * 1000;

                                if(ts[1] >= 'A') resNum = resNum + (ts[1] - 'A' + 10) * 100;  
                                else resNum = resNum + (ts[1] - '0') * 100;

                                if(ts[2] >= 'A') resNum = resNum + (ts[2] - 'A' + 10) * 10;  
                                else resNum = resNum + (ts[2] - '0') * 10;
                                
                                if(ts[3] >= 'A') resNum = resNum + (ts[3] - 'A' + 10);  
                                else resNum = resNum + (ts[3] - '0');
                                
                            } else {
                                printf("Error4 (residue number) reading PDB line : %s\n", lineBuffer);
                                exit(1);
                            }
                        }
                        
                        // chain ID : 1 char
                        chain = lineBuffer[21];

                        // X
                        strncpy(ts, lineBuffer+30, 8);
                        ts[8] = 0;
                        if(sscanf(ts, "%lf", &pP.x) != 1) { printf("Error5 (x-coord) reading PDB line : %s\n", lineBuffer); exit(1); }
                        // Y
                        strncpy(ts, lineBuffer+38, 8);
                        ts[8] = 0;
                        if(sscanf(ts, "%lf", &pP.y) != 1) { printf("Error6 (y-coord) reading PDB line : %s\n", lineBuffer); exit(1); }
                        // Z
                        strncpy(ts, lineBuffer+46, 8);
                        ts[8] = 0;
                        if(sscanf(ts, "%lf", &pP.z) != 1) { printf("Error7 (z-coord) reading PDB line : %s\n", lineBuffer); exit(1); }
                        
                        // occupancy
                        strncpy(ts, lineBuffer+56, 4);
                        ts[4] = 0;
                        if(sscanf(ts, "%lf", &occupancy) != 1) { printf("Error8 (occupancy) reading PDB line : %s\n", lineBuffer); exit(1); }
                        if(occupancy < 0.0  ||  occupancy > 1.0) { printf("Warning8 (occupancy) reading PDB line : %s\n", lineBuffer); occupancy = 1.0; }
                        
                        // B-factor
                        strncpy(ts, lineBuffer+60, 5);
                        ts[5] = 0;
                        if(sscanf(ts, "%lf", &Bfactor) != 1) { printf("Error9 (B-factor) reading PDB line : %s\n", lineBuffer); exit(1); }
                        
                        // element
                        strncpy(element, lineBuffer + 76, 2);
                        element[2] = '\0';
                        toUpperCase(element);
                        N = giveProtonNumber(element);
                        
                        // charge   e.g. " N1+"
                        strncpy(ts, lineBuffer + 78, 2);
                        ts[1] = '\0';
                        if(sscanf(ts, "%d", &charge) != 1) charge = 0;
                        else {
                                if(lineBuffer[79] == '-') charge = -charge;
                        }
			

                        // if the name or number is different, make a new prot.residue
                        if(strcmp(lastName, resName) != 0  ||  lastNum != resNum  ||  lastChain != chain ){
                                // set new name and number
                                strcpy(lastName, resName);
                                lastNum = resNum;
                                lastChain = chain;

                                if(flag == ATOM_LINE) {
                                        // encode into one letter amino acid code
                                        letter = residueName_C2L(resName);
                                        // create a new residue
                                        lastResidue = new Residue(resNum, letter, resName, chain, giveRotBonds(letter), lastResidue);
                                        // increase counter as there is a new residue
                                        residueCounter++;
                                }
                                else if(flag == HETATM_LINE) {
                                        if(strcmp(resName, "HOH") == 0  ||  strcmp(resName, "WAT") == 0  ||  strcmp(resName, "DOD") == 0  ||  strcmp(resName, "T3P") == 0) {
                                                lastWater = new Water(resNum, '~', resName, chain, lastWater);
                                                waterCounter++;
                                                letter = '~';
                                        }
                                        else if(strcmp(resName, "MSE") == 0  ||  strcmp(resName, "SEC") == 0) {
                                                // encode into one letter amino acid code
                                                letter = residueName_C2L(resName);
                                                // create a new residue
                                                lastResidue = new Residue(resNum, letter, resName, chain, giveRotBonds(letter), lastResidue);
                                                // increase counter as there is a new residue
                                                residueCounter++;
                                        }
                                        else {
                                                lastMolecule = new Molecule(resNum, 'X', resName, chain, lastMolecule);
                                                moleculeCounter++;
                                                letter = 'X';
                                                printf("new HET residue: %c%s %c %d\n", altLoc, resName, chain, resNum);
                                        }
                                }
                        }

                        // add atom
                        if(flag == ATOM_LINE) {
                                // add atom to the residue
                                lastResidue->addAtom(atomNum, atomType, altLoc, pP, occupancy, Bfactor, N, charge);
                        }
                        else if(flag == HETATM_LINE) {
                                // add water
                                if(letter == '~') {
                                        // add atom to the water
                                        lastWater->addAtom(atomNum, atomType, altLoc, pP, occupancy, Bfactor, N, charge);
                                }
                                // add molecule
                                else if (letter == 'X') {
                                        // add atom to the molecule
                                        lastMolecule->addAtom(atomNum, atomType, altLoc, pP, occupancy, Bfactor, N, charge);
                                }
                                // add Seleno residues
                                else if (letter == 'o'  ||  letter == 'u') {
                                        // add atom to the molecule
                                        lastResidue->addAtom(atomNum, atomType, altLoc, pP, occupancy, Bfactor, N, charge);
                                }
                                else {
                                        printf("Unknown structural element.\n");
                                        exit(1);
                                }
                        }
                        // increase the counter
                        atomCounter++;
		}
        }
        // END OF READING COORDINATES
        
        if(atomCounter == 0) {
                printf("No useful atoms were read... Strange.\n");
                exit(1);
        }
        if(residueCounter == 0) {
                printf("No protein atoms were read... Exit.\n");
                exit(1);
        }

        
        ProteinAtom *   pa2;
        LigandAtom *    la2;
        
        // initialize last pointer
        lastRawAtom = lastResidue->getLastAtom();
        
        // chain all atoms together
        curResidue = lastResidue;
        ra = NULL;
        while(curResidue) {
                pa = curResidue->getLastAtom();
                while(pa) {
                        // remember the last one for linking to other elements
                        ra = (RawAtom *) pa;
                        // if pointer to the previous is null, link to the previous one
                        if(pa->getPrevious() == NULL  &&  curResidue->getPrevious() != NULL) {
                                pa2 = curResidue->getPrevious()->getLastAtom();
                                ra2 = (RawAtom *) pa2;
                                ra->setPrevious(ra2);
                        }
                        // set species
                        pa->setSpecies(PROTEIN);
                        
                        ra->setDAC(atom_giveDonAcc(pa));
                        
                        pa = pa->getPrevious();
                }
                curResidue = curResidue->getPrevious();
        }

        curMolecule = lastMolecule;
        if(curMolecule != NULL) ra->setPrevious((RawAtom *) curMolecule->getLastAtom());
        while(curMolecule) {
                la = curMolecule->getLastAtom();
                while(la) {
                        // if pointer to the previous is null, link to the previous one
                        ra = (RawAtom *) la;
                        if(la->getPrevious() == NULL  &&  curMolecule->getPrevious() != NULL) {
                                la2 = curMolecule->getPrevious()->getLastAtom();
                                ra2 = (RawAtom *) la2;
                                ra->setPrevious(ra2);
                        }

                        // remember the last one for linking to other elements
                        la = la->getPrevious();
                }
                curMolecule = curMolecule->getPrevious();
        }
        
        curWater = lastWater;
        if(curWater != NULL) ra->setPrevious((RawAtom *) curWater->getLastAtom());
        while(curWater) {
                pa = curWater->getLastAtom();
                while(pa) {
                        // remember the last one for linking to other elements
                        ra = (RawAtom *) pa;
                        // if pointer to the previous is null, link to the previous one
                        if(pa->getPrevious() == NULL  &&  curWater->getPrevious() != NULL) {
                                pa2 = curWater->getPrevious()->getLastAtom();
                                ra2 = (RawAtom *) pa2;
                                ra->setPrevious(ra2);
                        }
                        // set species
                        pa->setSpecies(WATER);
                        
                        pa = pa->getPrevious();
                }
                curWater = curWater->getPrevious();
        }        
        	
	// SUMMARY
	printf("Read PDB file %s : %d atoms, %d residues, %d molecules, %d waters\n", filename, atomCounter, residueCounter, moleculeCounter, waterCounter);
        
        // close the template input file
	if (fclose(input) != 0) { printf("Error! Closing file %s failed.\n", filename); exit(1); }
}


// simple protein connectivity
void    protein_doFullConnectivity(Residue * rLast) {
    Residue *       r;
    Residue *       r2;
    ProteinAtom *   pa;
    ProteinAtom *   pa2;
    double          dist;

    r = rLast;
    while(r) {
        r2 = rLast;
        while(r2) {
            pa = r->getLastAtom();
            while(pa) {
                pa2 = r2->getLastAtom();
                while(pa2) {
                    if(pa != pa2) {
                        dist = giveDistance(pa, pa2);
                        // 0.528 is a scaling factor based onyears of experience and work very reliably for connectivity algorithms
                        if(dist  <  0.528 * (pa->getRadius() + pa2->getRadius())) {
                            pa->setNeighbor(pa2);
                        }
                        if(pa->getProtonNum() == SULPHUR  &&  pa2->getProtonNum() == SULPHUR  &&  dist  <  0.6 * (pa->getRadius() + pa2->getRadius())) {
                            pa->setNeighbor(pa2);
                        }
                    }
                    pa2 = pa2->getPrevious();
                }
                pa = pa->getPrevious();
            }
            r2 = r2->getPrevious();
        }
        r = r->getPrevious();
    }
}

// map intramolecular H-bonds
void    protein_mapIntramolecularHbond(Residue * r) {
        ProteinAtom *   pa;
        ProteinAtom *   pa2;
        ProteinAtom *   X;
        Residue *       rBkp = r;
        Residue *       r2 = NULL;
        double          d, angle;
        double          eHB;
        char            a1, a2;
        
        printf("Finding all H-bonds with E better than %4.1f kca;/mol\n", HB_MIN_ENERGY);
        // all vs all residues
        while(r) {
            pa = r->getLastAtom();
            while(pa) {
                r2 = rBkp;
                while(r2) {
                    pa2 = r2->getLastAtom();
                    while(pa2) {
                        if(pa->getProtonNum() == HYDROGEN  &&  pa->getNeighbor(0)->getProtonNum() != CARBON  &&  pa2->getDAC() == ACCEPTOR) {
                            
                            // prevent C=O...H-N within residue
                            if(r == r2  &&  strcmp(pa->getAtomType(), " H  ") == 0  &&  strcmp(pa2->getAtomType(), " O  ") == 0) {
                                // do nothing
                            }
                            else{
                                X = pa->getNeighbor(0);
                                eHB = giveHbondEnergy(X->getCoord(), pa->getCoord(), pa2->getCoord(), X->getProtonNum(), pa2->getProtonNum(), YES);

                                // if(d <= HB_MAX_DIST  &&  angle >= HB_THRESHOLD_ANGLE) {
                                if(eHB < HB_MIN_ENERGY) {
                                    if(strcmp(pa->getAtomType(), " H  ") == 0) a1 = 'b';
                                    else a1 = 's';
                                    
                                    if(strcmp(pa2->getAtomType(), " O  ") == 0  ||  strcmp(pa2->getAtomType(), " OXT") == 0) a2 = 'b';
                                    else a2 = 's';

                                    d = giveDistance(pa2->getCoord(), pa->getCoord());
                                    angle = giveAngle(pa2->getCoord(), pa->getCoord(), X->getCoord());

                                    printf("iHB:  %c%c  %c-H...%c  r=%4.2f  a=%5.1f  e=%5.2f  !  %4s %3s %3d ... %4s %3s %3d  !  %3s %3d ... %3s %3d\n", a1, a2,
                                            giveElement(X->getProtonNum()), giveElement(pa2->getProtonNum()), d, I80_DIV_PI * angle, eHB,
                                            pa->getAtomType(), r->getResidueName(), r->getUnitNumber(), pa2->getAtomType(), r2->getResidueName(), r2->getUnitNumber(),
                                            r->getResidueName(), r->getUnitNumber(), r2->getResidueName(), r2->getUnitNumber());
                                }
                            }
                        }
                        
                        pa2 = pa2->getPrevious();
                    }
                    r2 = r2->getPrevious();
                }
                pa = pa->getPrevious();
            }
            r = r->getPrevious();
        }

}

void    protein_mapIntramolecularHbond_LIGAND(Residue * r) {
        ProteinAtom *   pa;
        ProteinAtom *   pa2;
        ProteinAtom *   X;
        Residue *       rBkp = r;
        Residue *       r2 = NULL;
        double          d, angle;
        double          eHB;
        char            a1, a2;
        
        printf("Finding all L-L H-bonds with E better than %4.1f kcal/mol\n", HB_MIN_ENERGY);

        // clear Ehb
        while(r) {
            r->clearEhb();
            r = r->getPrevious();
        }
        
        r = rBkp;
        // all vs all residues
        while(r) {
            if(r->getActiveSite() == NO) {
                r2 = rBkp;
                while(r2) {
                    if(r2->getActiveSite() == NO) {
                        pa = r->getLastAtom();
                        while(pa) {
                            pa2 = r2->getLastAtom();
                            while(pa2) {
                                if(pa->getProtonNum() == HYDROGEN  &&  pa->getNeighbor(0)->getProtonNum() != CARBON  &&  pa2->getDAC() == ACCEPTOR) {

                                    // prevent C=O...H-N within residue
                                    if(r == r2  &&  strcmp(pa->getAtomType(), " H  ") == 0  &&  strcmp(pa2->getAtomType(), " O  ") == 0) {
                                        // do nothing
                                    } else {
                                        X = pa->getNeighbor(0);
                                        eHB = giveHbondEnergy(X->getCoord(), pa->getCoord(), pa2->getCoord(), X->getProtonNum(), pa2->getProtonNum(), YES);

                                        if(eHB < HB_MIN_ENERGY) {
                                            if(strcmp(pa->getAtomType(), " H  ") == 0) a1 = 'b';
                                            else a1 = 's';

                                            if(strcmp(pa2->getAtomType(), " O  ") == 0  ||  strcmp(pa2->getAtomType(), " OXT") == 0) a2 = 'b';
                                            else a2 = 's';

                                            d = giveDistance(pa2->getCoord(), pa->getCoord());
                                            angle = giveAngle(pa2->getCoord(), pa->getCoord(), X->getCoord());

                                            printf("iHB_L-L:  %c%c  %c-H...%c  r=%4.2f  a=%5.1f  e=%5.2f  !  %4s %3s %3d ... %4s %3s %3d  !  %3s %3d ... %3s %3d\n", a1, a2,
                                                    giveElement(X->getProtonNum()), giveElement(pa2->getProtonNum()), d, I80_DIV_PI * angle, eHB,
                                                    pa->getAtomType(), r->getResidueName(), r->getUnitNumber(), pa2->getAtomType(), r2->getResidueName(), r2->getUnitNumber(),
                                                    r->getResidueName(), r->getUnitNumber(), r2->getResidueName(), r2->getUnitNumber());
                                        }
                                        
                                        if(r->getID() != '~'  &&  r2->getID() != '~') {
                                            r->addEhb(eHB);
                                            r2->addEhb(eHB);
                                        }
                                    }
                                }

                                pa2 = pa2->getPrevious();
                            }
                            pa = pa->getPrevious();
                        }
                    }    
                    r2 = r2->getPrevious();
                }
            }
            r = r->getPrevious();
        }
        
        printf("Final EhbSum_L-L scores per residue:\n");
        r = rBkp;
        while(r) {
            if(r->getActiveSite() == NO) printf("ehb_L-L %3s %3d : %5.2f\n", r->getResidueName(), r->getUnitNumber(), r->getEhbSum()); 
            r = r->getPrevious();
        }

}

int     isWater(Residue * r) {
        if(strcmp(r->getResidueName(), "T3P") == 0  || strcmp(r->getResidueName(), "WAT") == 0  || strcmp(r->getResidueName(), "HOH") == 0  || strcmp(r->getResidueName(), "DOD") == 0) {
            return YES;
        }
        return NO;
}

void    protein_mapIntermolecularHbond(Residue * r) {
        ProteinAtom *   pa;
        ProteinAtom *   pa2;
        ProteinAtom *   X;
        Residue *       rBkp = r;
        Residue *       r2 = NULL;
        double          d, angle;
        double          eHB;
        char            a1, a2;
        
        printf("Finding all L-P H-bonds with E better than %4.1f kcal/mol\n", HB_MIN_ENERGY);
        
        // clear Ehb
        while(r) {
            r->clearEhb();
            r = r->getPrevious();
        }
        
        r = rBkp;
        // all vs all residues
        while(r) {
            r2 = rBkp;
            while(r2) {
                // if from different worlds
                if(r->getActiveSite() + r2->getActiveSite() == YES) {
                    pa = r->getLastAtom();
                    while(pa) {
                        pa2 = r2->getLastAtom();
                        while(pa2) {
                            if(pa->getProtonNum() == HYDROGEN  &&  pa->getNeighbor(0)->getProtonNum() != CARBON  &&  pa2->getDAC() == ACCEPTOR) {

                                // prevent C=O...H-N within residue
                                if(r == r2  &&  strcmp(pa->getAtomType(), " H  ") == 0  &&  strcmp(pa2->getAtomType(), " O  ") == 0) {}
                                // ignore Wat - Wat bridges
                                else if(isWater(r) == YES  &&  isWater(r2) == YES) {}
                                else {
                                    X = pa->getNeighbor(0);
                                    eHB = giveHbondEnergy(X->getCoord(), pa->getCoord(), pa2->getCoord(), X->getProtonNum(), pa2->getProtonNum(), YES);

                                    if(eHB < HB_MIN_ENERGY) {
                                        if(strcmp(pa->getAtomType(), " H  ") == 0) a1 = 'b';
                                        else a1 = 's';

                                        if(strcmp(pa2->getAtomType(), " O  ") == 0  ||  strcmp(pa2->getAtomType(), " OXT") == 0) a2 = 'b';
                                        else a2 = 's';

                                        d = giveDistance(pa2->getCoord(), pa->getCoord());
                                        angle = giveAngle(pa2->getCoord(), pa->getCoord(), X->getCoord());

                                        printf("HB_L-P:  %c%c  %c-H...%c  r=%4.2f  a=%5.1f  e=%5.2f  !  %4s %3s %3d ... %4s %3s %3d  !  %3s %3d ... %3s %3d\n", a1, a2,
                                                giveElement(X->getProtonNum()), giveElement(pa2->getProtonNum()), d, I80_DIV_PI * angle, eHB,
                                                pa->getAtomType(), r->getResidueName(), r->getUnitNumber(), pa2->getAtomType(), r2->getResidueName(), r2->getUnitNumber(),
                                                r->getResidueName(), r->getUnitNumber(), r2->getResidueName(), r2->getUnitNumber());
                                    }
                                    
                                    if(r->getID() != '~'  &&  r2->getID() != '~') {
                                        r->addEhb(eHB);
                                        r2->addEhb(eHB);
                                    }
                                }
                            }

                            pa2 = pa2->getPrevious();
                        }
                        pa = pa->getPrevious();
                    }
                }   
                r2 = r2->getPrevious();
            }
            r = r->getPrevious();
        }
        
        printf("Final hbSum_L-P scores per residue LIGAND to PROTEIN:\n");
        r = rBkp;
        while(r) {
            if(r->getActiveSite() == NO) printf("ehb_L-P %3s %3d : %5.2f\n", r->getResidueName(), r->getUnitNumber(), r->getEhbSum()); 
            r = r->getPrevious();
        }
}

double  give_e_hphob(ProteinAtom * a1, ProteinAtom * a2, double d) {
    double  scale, e = 0;

    scale = 2.0 * (d - a1->getRadius() - a2->getRadius() - 2.0) / 3.0;
    
          if(scale <= -1.0) e = EHPHOB;
    else if (scale <   1.0) e = EHPHOB * (0.25 * scale * scale * scale  -  0.75 * scale  +  0.5);
    else                    e = 0;
    
    return e;
}

int     isLipophilic(ProteinAtom  * pa){
    if(pa->getProtonNum() == CARBON) {
        // omit backbone
        if(strcmp(pa->getAtomType(), " C  ") == 0) return NO;
        return YES;
    }
    else if(pa->getProtonNum() == SULPHUR) return YES;
    else return NO;
}

int  check_1_4(ProteinAtom * pa, ProteinAtom * pa2) {
    int             i, j, k;
    ProteinAtom *   n1;
    ProteinAtom *   n2;
    ProteinAtom *   n3;
    
    // loop over all neighbors of pa
    for (i = 0; i < pa->getNeighborCount(); i++) {
        // 1 - 2
        n1 = pa->getNeighbor(i);
        

        // check if n1 = pa2
        if (n1->getNumber() == pa2->getNumber()) return YES;
        
        // loop over all neighbors of n1
        for (j = 0; j < n1->getNeighborCount(); j++) {
            // 1 - 3
            n2 = n1->getNeighbor(j);
            if (n2 == pa2) return YES;

            // loop over all neighbors of n2
            for (k = 0; k < n2->getNeighborCount(); k++) {
                // 1 - 4
                n3 = n2->getNeighbor(k);
                if (n3 == pa2) return YES;
            }
        }
    }
    return NO;
}

void    protein_mapIntermolecularLipo(Residue * r) {
        ProteinAtom *   pa;
        ProteinAtom *   pa2;
        Residue *       rBkp = r;
        Residue *       r2 = NULL;
        double          d;
        double          eLipo;
        double          eLipoSum;
        
        // clear Lipo
        while(r) {
            r->clearLipo();
            r = r->getPrevious();
        }
        
        printf("Finding all L-P residue pairs with E better than %4.1f kca;/mol\n", LIPO_MIN_ENERGY);
        // all vs all residues
        r = rBkp;
        while(r) {
            // avoid duplicating
            r2 = r->getPrevious();
            while(r2) {
                eLipoSum = 0;

                // in one activesite = YES (1) and one activesite = NO (0)
                // 0 + 0 = 0
                // 1 + 0 = 1
                // 1 + 1 = 2
                if(r->getActiveSite() + r2->getActiveSite() == YES) {
                    pa = r->getLastAtom();
                    while(pa) {
                        pa2 = r2->getLastAtom();
                        while(pa2) {
                            if(isLipophilic(pa) == YES  &&  isLipophilic(pa2) == YES) {
                                d = giveDistance(pa2->getCoord(), pa->getCoord());
                                eLipo = give_e_hphob(pa, pa2, d);
                                eLipoSum = eLipoSum + eLipo;
                            }

                            pa2 = pa2->getPrevious();
                        }
                        pa = pa->getPrevious();
                    }
                    if (eLipoSum < LIPO_MIN_ENERGY) {
                        printf("Lipo:  e=%5.2f  !  %3s %3d ... %3s %3d\n",
                                eLipoSum,
                                r->getResidueName(), r->getUnitNumber(), r2->getResidueName(), r2->getUnitNumber());

                    }
                    r->addElipo(eLipoSum);
                    r2->addElipo(eLipoSum);
                }

                r2 = r2->getPrevious();
            }
            r = r->getPrevious();
        }
        
        printf("Final lipoSum_L-P scores per residue LIGAND to PROTEIN:\n");
        r = rBkp;
        while(r) {
            // print all because of automated processing
            if(r->getActiveSite() == NO) printf("lipo_L-P %3s %3d : %5.2f\n", r->getResidueName(), r->getUnitNumber(), r->getElipoSum()); 
            r = r->getPrevious();
        }
}

void    protein_mapIntramolecularLipo_LIGAND(Residue * r) {
        ProteinAtom *   pa;
        ProteinAtom *   pa2;
        Residue *       rBkp = r;
        Residue *       r2 = NULL;
        double          d;
        double          eLipo;
        double          eLipoSum;
        
        // clear Lipo
        while(r) {
            r->clearLipo();
            r = r->getPrevious();
        }
        
        printf("Finding all L-L residue pairs with E better than %4.1f kcal/mol\n", LIPO_MIN_ENERGY);
        // all vs all residues
        r = rBkp;
        while(r) {
            if(r->getActiveSite() == NO) {
                // avoid duplicating
                r2 = r->getPrevious();
                while(r2) {
                    if(r2->getActiveSite() == NO) {
                        eLipoSum = 0;

                        pa = r->getLastAtom();
                        while(pa) {
                            pa2 = r2->getLastAtom();
                            while(pa2) {
                                if(isLipophilic(pa) == YES  &&  isLipophilic(pa2) == YES  &&  check_1_4(pa, pa2) == NO) {
                                    // if consecutive residues the switch off backbone
                                    if(abs((r->getUnitNumber() - r2->getUnitNumber()) == 1)  &&  strcmp(pa->getAtomType(), " CA ") == 0  &&  strcmp(pa2->getAtomType(), " CA ") == 0) {
                                         // printf("Skipping %d and %d\n", r->getUnitNumber(), r2->getUnitNumber());
                                    }
                                    else {
                                        d = giveDistance(pa2->getCoord(), pa->getCoord());
                                        eLipo = give_e_hphob(pa, pa2, d);
                                        eLipoSum = eLipoSum + eLipo;
                                    }
                                }

                                pa2 = pa2->getPrevious();
                            }
                            pa = pa->getPrevious();
                        }
                        if (eLipoSum < LIPO_MIN_ENERGY) {
                            printf("lipo_E:  e=%5.2f  !  %3s %3d ... %3s %3d\n",
                                    eLipoSum,
                                    r->getResidueName(), r->getUnitNumber(), r2->getResidueName(), r2->getUnitNumber());

                        }
                        r->addElipo(eLipoSum);
                        r2->addElipo(eLipoSum);
                    }
                    r2 = r2->getPrevious();
                }
            }
            r = r->getPrevious();
        }
        
        printf("Final ElipoSum_L-L scores per residue:\n");
        r = rBkp;
        while(r) {
            // if(r->getElipoSum() <  LIPO_MIN_ENERGY) printf("lipo_L-L %3s %3d : %5.2f\n", r->getResidueName(), r->getUnitNumber(), r->getElipoSum()); 
            if(r->getActiveSite() == NO) printf("lipo_L-L %3s %3d : %5.2f\n", r->getResidueName(), r->getUnitNumber(), r->getElipoSum()); 
            r = r->getPrevious();
        }
}

void    protein_mark_all(Residue * r) {
    while(r) {
        r->setActiveSite();
        r = r->getPrevious();
    }
}

void    protein_clear_ligand(Residue * r) {
    Residue *   rBkp = r;
    Residue *   r2 = r;
    ProteinAtom *   pa;
    // ProteinAtom *   pa2;
    int         cleared = YES;
    int         i;
    
    while(r) {
        // find the residue belonging to the ligand
        if(strcmp(r->getResidueName(), "TRP") == 0  &&  (r->getUnitNumber() == 7  ||  r->getUnitNumber() == 8)) {
            // 
            printf("Ligand marking residue found: %s %d\n", r->getResidueName(), r->getUnitNumber());
            r->clearActiveSite();
        }
        r = r->getPrevious();
    }
    
    printf("Residue clearing\n");
    while(cleared == YES) {
        cleared = NO;
        r = rBkp;
        while(r) {
            if(r->getActiveSite() == NO) {
                // clear all of its atoms
                pa = r->getLastAtom();
                while(pa){
                    // look all neighbors of all cleared residues
                    for(i = 0; i < pa->getNeighborCount(); i++) {
                        r2 = (Residue *)pa->getNeighbor(i)->getResidue();
                        if(r2->getActiveSite() == YES) {
                            r2->clearActiveSite();
                            printf("%s %3d set as ligand\n", r2->getResidueName(), r2->getUnitNumber());
                            cleared = YES;
                        }
                    }
                    pa = pa->getPrevious();
                }
            }
            r = r->getPrevious();
        } 
    }
    printf("Done\n");
}

// main prougram routine
int main(int argc, char *argv[]) {
	printf("interaction-analyzer ver. %4.2f by Martin Smiesko\n", VERSION);

        Molecule * 	lastRawLigand = NULL;
        Water *         water = NULL;
        Residue *       residue = NULL;
        RawAtom *       lastRawAtom = NULL;
        
        // if a datafile supplied on the command line
        if(argc == 2  &&  strcmp(argv[1], "")) {
            
                // R E A D   P D B
                // read protein-ligand complex - MD frame in the PDB format
                readProteinPDB(argv[1], residue, lastRawLigand, water, lastRawAtom);
                
                
                
                // P R E P A R E   T H E   D A T A
                // do protein & ligand connectivity
                protein_doFullConnectivity(residue);
                
                // mark all as active site
                protein_mark_all(residue);
                
                // Clear ligand flag based on residue TRP 7 or 8 - compstatin or CP01 for the study of Nature
                protein_clear_ligand(residue);
                
                
                
                // A N A L Y Z E   I N T E R A C T I O N S
                // detect intramolecular H-bonds - LIGAND
                protein_mapIntramolecularHbond_LIGAND(residue);
                
                // detect intermolecular H-bonds  L...P
                protein_mapIntermolecularHbond(residue);
                
                // detect intramolecular lipophilic interactions - LIGAND
                protein_mapIntramolecularLipo_LIGAND(residue);
                                
                // detect intermolecular lipophilic interactions L...P
                protein_mapIntermolecularLipo(residue);

                
                // all fine :)
                printf("interaction-analyzer finished gracefully.\n");
	}
	else { 
		// if not formated properly on the command line, print this short help
                printf("\nProgram to perform H-bonding and Lipophilic interaction analyses on the PDB snapshots from MDs.\n");
		printf("\nUsage: ./interaction-analyzer.x  PDBfile\n");
	}

	return 0;
}
