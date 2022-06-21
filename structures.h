/* 
 * File:   structures.h
 * Author: Martin Smiesko, University of Basel
 *
 * Created on June 19, 2014, 4:57 PM
 */

#ifndef STRUCTURES_H
#define	STRUCTURES_H

class   StructuralUnit;

class RawAtom {
protected:
        RawAtom *               previous;
	Point			coord;		// coordinates
	unsigned char		protonNum;	// proton number of the atom
	char			atomType[5];	// atom label = PDBX atom type
	double			radius;		// VdW radius in Angstroms
	double			charge;		// electrostatic charge
	int                     number;		// atom number as appears in the source PDB file
	char                    hybridization;	// hybridization state
        char                    dac;            // donor = 1, acceptor = 2, nothing = 0
	unsigned char           neighborCount;	// number of neighbors of the atom
        unsigned char           metalNeighborCount;	// number of neighbors of the atom
        StructuralUnit *        residue;        // residue pointer to which it belongs
        
        double                  occupancy;
        double                  Bfactor;
        char                    altLocation;
        char                    species;        // protein, water, ligand, metal
public:
	void			setCoord(Point p) { coord = p; }
	void			setNumber(int n) { number = n; }
	void			setCharge(double q) { charge = q; }
        void                    setHybridization(char value) { hybridization = value; }
        void                    setDAC(char value) { dac = value; }
        void                    setResidue(StructuralUnit * r) { residue = r; }
        void                    setRadius(double r) { radius = r; }
        void                    setOccupancy(double value) { occupancy = value; }
        void                    setBfactor(double value) { Bfactor = value; }
        void                    setAltLocation(char value) { altLocation = value; }
        void                    setPrevious(RawAtom * value) { previous = value;  }
        void                    setSpecies(char value) { species = value; }

	int             	getNumber() { return number; }
	Point			getCoord() { return coord; }
	unsigned char		getProtonNum() { return protonNum; }
	double			getRadius() { return radius; }
	double			getCharge() { return charge; }
	char *			getAtomType() { return atomType; }
	char                    getHybridization() { return hybridization; }
        RawAtom *               getPrevious() { return previous; }
        char                    getDAC() { return dac; }
	unsigned char           getNeighborCount() { return neighborCount; }
        unsigned char           getMetalNeighborCount() { return metalNeighborCount; }
        StructuralUnit *        getResidue() { return residue; }
        double                  getOccupancy() { return occupancy; }
        double                  getBfactor() { return Bfactor; }
        char                    getAltLocation() { return altLocation; }
        char                    getSpecies() { return species; }
        
        void                    incNeighborCount() { neighborCount++; };

        // number  atom_type  alternate_location  x  y  z  occupancy  B-factor  proton_num  charge  index
        RawAtom(int n, char * at, char aL, Point p, double occ, double Bf, unsigned char pN, double Q) {
 		number = n;
		strcpy(atomType, at);
		coord = p;
		charge = Q;
                previous = NULL;

                protonNum = pN;
                radius = giveVdWRadius(pN);

                hybridization = UNCLEAR;
                dac = NO;
                neighborCount = 0;
                metalNeighborCount = 0;
                residue = NULL;
                
                altLocation = aL;
                occupancy = occ;
                Bfactor = Bf;
        };
        
        RawAtom() {};
};

class	Ring;

class	LigandAtom : public RawAtom {
	LigandAtom *	previous;	// pointer to the previous atom in the list
	LigandAtom *	next;		// pointer to the next atom in the list
	unsigned char	inRing;		// flag on, if atom already in a ring
	unsigned char	written;	// flag if the atom has been written in a ring searching algorithm
	LigandAtom *	neighbor[ATOM_MAX_NEIGHBORS];	// pointers to neighboring atoms
	char		bond[ATOM_MAX_NEIGHBORS];	// bond order, bonds in the same order like neighbors
	char		formalQ;	// formal charge
	Ring *		ring[ATOM_MAX_NEIGHBORS];	// pointers to rings, where this atom belong to
        unsigned char   nonPolarH;
        unsigned char   rotatableH;
        RawAtom *       metalNeighbor[METAL_MAX_NEIGHBORS];
        int             molNum;
        double          surfaceArea;
        double          polarArea;
        LigandAtom *    original;
        LigandAtom *    linked;         // atom linked to another residue
        char            neighT;         // count of neighbors as appeared in the template file
        int             fga;
        
	public:
	void		setNext(LigandAtom * a) { next = a; }
	void		setRingFlag(Ring * r) { ring[inRing] = r; inRing++; };
	void		setWritten() { written = YES; }
	void		clearWritten() { written = NO; }
	void		setNeighbor(LigandAtom * ptr) {
                                // security routine
                                int     i;
                                for(i = 0; i < neighborCount; i++)
                                        if(neighbor[i] == ptr) { 
                                                return;
                                        }
                
				neighbor[neighborCount] = ptr;
                                bond[neighborCount] = UNCLEAR;
				neighborCount++;
				if(neighborCount > ATOM_MAX_NEIGHBORS) { printf("Error! Too many atoms connected to atom %3d\n",  neighborCount); exit(1); }
			}
	void		setNeighborAndBond(LigandAtom * ptr, char bondOrder) {
                                // security routine
                                int     i;
                                for(i = 0; i < neighborCount; i++) 
                                        if(neighbor[i] == ptr) { 
                                                return;
                                        }
                
				neighbor[neighborCount] = ptr;
                                bond[neighborCount] = bondOrder;
				neighborCount++;
				if(neighborCount > ATOM_MAX_NEIGHBORS) { printf("Error! Too many atoms connected to atom %3d\n",  neighborCount); exit(1); }
			}
	void		setMetalNeighbor(RawAtom * ptr) {
                                // security routine
                                int     i;
                                for(i = 0; i < neighborCount; i++) 
                                        if(neighbor[i] == ptr) { 
                                                return;
                                        }
                
				metalNeighbor[metalNeighborCount] = ptr;
				metalNeighborCount++;
				if(metalNeighborCount > METAL_MAX_NEIGHBORS) { printf("Error! Too many atoms connected to metal %3d\n",  metalNeighborCount); exit(1); }
			}        
	void		setBond(unsigned char n, char value) { bond[n] = value; }
	void		incFormalQ() { formalQ++; }
	void		decFormalQ() { formalQ--; }
        void            setVdW(double value) { radius = value; }
        void            setNonPolarH(unsigned char v) { nonPolarH = v; }
        void            setRotatableH(unsigned char v) { rotatableH = v; }
        void            setMolNum(int value) { molNum = value; }
        void            setOriginal(LigandAtom * value) { original = value; } 
        void            setLinked(LigandAtom * v) { linked = v; }
        void            setAtomType(char * v) { strcpy(atomType, v); }
        void            setNeighT(char v) { neighT = v; }
        void            setFGA(int v) { fga = v; }
        
	unsigned char	getRingFlag() { return inRing; };	
	unsigned char	getWritten() { return written; }
	LigandAtom *	getPrevious() { return previous; }
	LigandAtom *	getNext() { return next; }
	LigandAtom *	getNeighbor(unsigned char i) { return neighbor[i]; }
        RawAtom *	getMetalNeighbor(unsigned char i) { return metalNeighbor[i]; }
	char		getBond(unsigned char n) { return bond[n]; }
	char		getFormalQ() { return formalQ; }
	Ring *		getNthRing(unsigned char n) { return ring[n]; }
        int             getMolNum() { return molNum; }
        double          getPolarArea() { return polarArea; }
        double          getSurfaceArea() { return surfaceArea; }
        LigandAtom *	getOriginal() { return original; }
        LigandAtom *    getLinked() { return linked; }
        char            getNeighT() { return neighT; }
        unsigned char	getRotatableH() { return rotatableH; };

        unsigned char	isNonPolarH() { return nonPolarH; };
        void            addSurfaceArea(double value) { surfaceArea += value; }
        void            addPolarArea(double value) { polarArea += value; }
        int             getFGA() { return fga; }
        

        LigandAtom(int n, char * at, char aL, Point p, double occ, double Bf, unsigned char pN, double Q, LigandAtom * prev) : RawAtom(n, at, aL, p, occ, Bf, pN, Q) {
                for(neighborCount = 0; neighborCount < ATOM_MAX_NEIGHBORS; neighborCount++) {
                        neighbor[neighborCount] = NULL;
                        bond[neighborCount] = UNCLEAR;
                        ring[neighborCount] = NULL;
                }
                for(neighborCount = 0; neighborCount < METAL_MAX_NEIGHBORS; neighborCount++) {
                        metalNeighbor[neighborCount] = NULL;
                }
                previous = prev;
		if(prev != NULL) prev->setNext(this);
		next = NULL;
        	neighborCount = 0;
		inRing = NO;
		written = NO;
		formalQ = 0;
                nonPolarH = YES;
                this->setPrevious((RawAtom *) prev);
                
                if(isMetal(pN)) {
                        species = METAL;
                }
                else species = LIGAND;
                dac = NO;
                molNum = -1;
                surfaceArea = polarArea = 0;
                original = NULL;
                linked = NULL;
                neighT = -1;
                rotatableH = NO;
                fga = 0;
	}
        
        LigandAtom() : RawAtom() { fga = 0; }
};

class	ProteinAtom : public RawAtom {
	ProteinAtom *	previous;	// pointer to the previous atom in the list
	ProteinAtom *	next;		// pointer to the next atom in the list
        unsigned char   Hcount;
	ProteinAtom *	neighbor[ATOM_MAX_NEIGHBORS];	// pointers to neighboring atoms
        int             flag;           // multipurpose flag variable - clear before use :)   
	
	public:
	void		setNext(ProteinAtom * a) { next = a; }
        void            setHcount(unsigned char value) { Hcount = value; }
        void            setVdW(double value) { radius = value; }
        void            setFlag(int value) { flag = value; }

        unsigned char   getHcount() { return Hcount; }
	ProteinAtom *	getPrevious() { return previous; }
	ProteinAtom *	getNext() { return next; }
        ProteinAtom *	getNeighbor(unsigned char i) { return neighbor[i]; }
        int             getFlag() { return flag; }

	void		setNeighbor(ProteinAtom * ptr) {
                                // security routine
                                int     i;
                                for(i = 0; i < neighborCount; i++)
                                        if(neighbor[i] == ptr) {
                                                return;
                                        }
                
				neighbor[neighborCount] = ptr;
				neighborCount++;
				if(neighborCount > ATOM_MAX_NEIGHBORS) { printf("Error! Too many atoms (%3d) connected to protein atom %s %d\n",  neighborCount, this->atomType, this->getNumber()); exit(1); }
			}
        // number  atom_type  alternate_location  x  y  z  occupancy  B-factor  proton_num  charge  index
        ProteinAtom(int n, char * at, char aL, Point p, double occ, double Bf, unsigned char pN, double Q, ProteinAtom * prev) : RawAtom(n, at, aL, p, occ, Bf, pN, Q) {
		previous = prev;
		if(prev != NULL) prev->setNext(this);
		next = NULL;
                neighborCount = 0;
                Hcount = -1;
                this->setPrevious((RawAtom *) prev);
                flag = 0;
                radius = giveVdWRadius(pN);
	}
};

class	StructuralUnit {
protected:
	char			ID;			// one letter identifier of an amino acid or a residue, t for template, l for ligand
	int             	unitNumber;		// the unit number as read from the file
	unsigned char		activesite;		// active site flag
	ProteinAtom *           lastAtom;		// pointer to the last atom of the molecule
	ProteinAtom *		firstAtom;              // pointer to the first atom of the molecule
        char                    residueName[4];         // residue name long
        char                    chain;
public:
	void			setActiveSite() { activesite = YES; }
        void			clearActiveSite() { activesite = NO; }
	ProteinAtom *           getLastAtom() { return lastAtom; }
	ProteinAtom *           getFirstAtom() { return firstAtom; }
	int                     getUnitNumber() { return unitNumber; }
	char			getID() { return ID; }
        char *			getResidueName() { return residueName; }
	unsigned char		getActiveSite() { return activesite; }
	char                    getChain() { return chain; }
};

class	Residue : public StructuralUnit {
	Residue *		previous;	// pointer to the previous residue
	Residue *		next;		// pointer to the next residue
	Residue *		lastConf;	// pointer to the last conformation of this residue
	unsigned char		rotBonds;	// number of rotatable bonds
	int                     cloneN;		// number of clones of this residue
	int			bumpsToLigand;	// number of bumps a residue has with the ligand
	int			bumpsToProtein;	// number of bumps a residue has with all other residues including all atoms; -1 if not yet calculated
        Residue *               master;         // pointer to the master residue, NULL for the master itself
        int                     confID;         // conformation number
        unsigned char           bridged;        // flag indicating that the residue is H-bonded to a residue in the active site
        unsigned char           missAtom;       // missing atoms flag
        double                  ElipoSum;       // sum of lipophilic interaction energies to others
        double                  EhbSum;         // sum of H-bonding interaction energies to others
	
	
	public:
	void			setPrevious(Residue * r) { previous = r; }
	void			setNext(Residue * r) { next = r; }
	void			setBumpsToLigand(int b) { bumpsToLigand = b; }
	void			setRotBonds(unsigned char n) { rotBonds = n; }
        void                    setMaster(Residue * r) { master = r; }
        void                    setConfID(int n) { confID = n; }
        void                    setBridged() { bridged = YES; }
        void                    setLastconf(Residue * r) { lastConf = r; }
        void                    setMissAtom() { missAtom++; }
        void                    setID(char c) { ID = c; }
        
	Residue *		getPrevious() { return previous; }
	Residue *		getNext() { return next; }
	Residue * 		getLastConf() { return lastConf; }
	unsigned char		getRotBonds() { return rotBonds; }
	int	 		getBumpsToLigand() { return bumpsToLigand; }
	unsigned int		getCloneN() { return cloneN; }
        Residue *               getMaster() { return master; }
        int                     getConfID() { return confID; }
        unsigned char           getBridged() { return bridged; }
        unsigned char           getMissAtom() { return missAtom; }
        double			getElipoSum() { return ElipoSum; }
        double			getEhbSum() { return EhbSum; }

        void			decreaseCloneN(){ cloneN--; }
        void                    addElipo(double el) { ElipoSum = ElipoSum + el; }
        void                    addEhb(double ehb) { EhbSum = EhbSum + ehb; }
        void                    clearLipo() { ElipoSum = 0; } 
        void                    clearEhb() { EhbSum = 0; } 
        void			addAtom(int num, char * at, char aL, Point p, double occ, double Bf, unsigned char protNum, double Q){
                                        lastAtom = new ProteinAtom(num, at, aL, p, occ, Bf, protNum, Q, lastAtom);
                                        lastAtom->setResidue(this);
                                        lastAtom->setSpecies(PROTEIN);
					// if it is the first atom
					if(lastAtom->getPrevious() == NULL) firstAtom = lastAtom;
                                }
     
        // constructor		
	Residue(int n, char l, char * rname, char ch, unsigned char rotBs, Residue * p) {
		unitNumber = n;
		ID = l;
		previous = p;
		if(p != NULL) p->setNext(this);
		next = NULL;
		rotBonds = rotBs;
		activesite = NO;
		lastConf = NULL;
		cloneN = 0;
		lastAtom = NULL;
		bumpsToProtein = -1;
		bumpsToLigand = 0;
                master = NULL;
                confID = 0;
                bridged = NO;
                strcpy(residueName, rname);
                residueName[3] = 0;
                chain = ch;
                missAtom = 0;
                ElipoSum = 0;
                EhbSum = 0;
	}
};

class   ResidueList {
protected:
        ResidueList *           previous;
        Residue *		residue;
public:
	ResidueList *		getPrevious() { return previous; }
	Residue *		getResidue() { return residue; }

        ResidueList(ResidueList * p, Residue * r) { previous = p; residue = r; }
        ResidueList(Residue * r) { residue = r; }
};

class Water : public StructuralUnit {
	Water *		previous;
	Water *		next;
	unsigned char	structural;
        unsigned char	bridging;
        unsigned char   solvate;        // water which has 2 and more H-bonds to the same ligand
        unsigned char   status;         // YES = on, NO = off
        unsigned char   bond2W;         // bound to other water
        unsigned char   bond2L;         // bound to ligand
        unsigned char   bond2P;         // bound to protein
        unsigned char   bond2M;         // bound to metal
        unsigned char   flag;

	public:
                
        // number  atom_type  alternate_location  x  y  z  occupancy  B-factor  proton_num  charge  index
        void		addAtom(int num, char * at, char aL, Point p, double occ, double Bf, unsigned char protNum, double Q){                
                                lastAtom = new ProteinAtom(num, at, aL, p, occ, Bf, protNum, Q, lastAtom);
                                lastAtom->setResidue(this);
                                lastAtom->setSpecies(WATER);
                                lastAtom->setDAC(DON_ACC);
                        }

	void		setStructural() { structural = YES; }
        void		setSolvate(unsigned char v) { solvate = v; }
        void		setBridging() { bridging = YES; }
	void		setNext(Water * w) { next = w; }
        void            setStatus(unsigned char s) { status = s; }
        void            setFlag(unsigned char v) { flag = v; }
        
        void            addBondToWater() { bond2W++; }
        void            addBondToLigand() { bond2L++; }
        void            addBondToProtein() { bond2P++; }
        void            addBondToMetal() { bond2M++; }
        
	unsigned char	getStructural() { return structural; }
        unsigned char	getSolvate() { return solvate; }
        unsigned char	getBridging() { return bridging; }
	Water * 	getPrevious() { return previous; }
	Water * 	getNext() { return next; }
        unsigned char   getStatus() { return status; }
        unsigned char   getBonds2Water() { return bond2W; }
        unsigned char   getBonds2Ligand() { return bond2L; }
        unsigned char   getBonds2Protein() { return bond2P; }
        unsigned char   getBonds2Metal() { return bond2M; }
        unsigned char   getFlag() { return flag; }
	
        Water(int n, char l, char * rname, char ch, Water * p) {
		unitNumber = n;
                structural = NO;
		previous = p;
		lastAtom = NULL;
		next = NULL;
		if(p != NULL) p->setNext(this);
		ID = l;
                status = YES;
                strcpy(residueName, rname);
                residueName[3] = 0;
                chain = ch;
                bond2W = bond2L = bond2P = bond2M = 0;
                activesite = 0;
                bridging = NO;
                flag = 0;
                solvate = NO;
	}
};

class	Ring {
	unsigned char	size;
	unsigned char	lipophilic;	// set to 1 if all atoms are carbons or sulphur, and no substitutents
	unsigned char	linkerRing;	// set if the ring has more than 1 chain connected
	Point		centroid;
	Point		norm;
	LigandAtom *	member[MAX_RING_SIZE];
	Ring *		previous;
	Ring *		fusedto[MAX_RING_SIZE];
	unsigned char	aromatic;
	unsigned char	bondsdone;
	unsigned char	fusedCount;
	unsigned char	number;		// ring number as found by the ring finding routine
        unsigned char   dd;             // delocalized and double count
        unsigned char   ddouts;         // delocalized and double count
        
	public:
	void		setSize(unsigned char value) { size = value; }
	void		setCentroid(Point value) { centroid = value; }
	void		setNormal(Point v1) { norm = v1; }
	void		setMember(LigandAtom * value) { 
				if(size < MAX_RING_SIZE) {
                                        member[size] = value;
                                        size++;
                                }
                                else {
					printf("Error: Maximum ring size exceeded (n = %d).\n", MAX_RING_SIZE);
					exit(1);
				}
			}
	void		setFused(Ring * value) { 
				fusedto[fusedCount] = value;
				fusedCount++;
				if(fusedCount > MAX_RING_SIZE) {
					printf("Error: Maximum number of fused rings exceeded (n = %d).\n", MAX_RING_SIZE);
					exit(1);
				}
			}
	void		setLipophilic(unsigned char value) { lipophilic = value; };
	void		setLinkerRing(unsigned char value) { linkerRing = value; };
	void		setAromatic() { aromatic = YES; }
	void		setBondsDone() { bondsdone = YES; }
        void            setDelocAndDouble(unsigned char value) { dd = value; }
        void            setDelocAndDoubleOuts(unsigned char value) { ddouts = value; }
		
	unsigned char	getSize() { return size; }
	Point		getCentroid() { return centroid; }
	LigandAtom *	getMember(unsigned char value) { return member[value]; }
	unsigned char	getFusedCount() { return fusedCount; }
	Ring *		getFused(unsigned char value) { return fusedto[value]; }
	Ring *		getPrevious() { return previous; }
	unsigned char	getLipophilic() { return lipophilic; }
	unsigned char	getLinkerRing() { return linkerRing; }
	Point		getNorm() { return norm; }
	unsigned char	getAromatic() { return aromatic; }
	unsigned char	getBondsDone() { return bondsdone; }
	unsigned char	getNumber() { return number; }
        unsigned char	getDelocAndDouble() { return dd; }
        unsigned char	getDelocAndDoubleOuts() { return ddouts; }

	Ring(Ring * p) { 
		previous = p;
		size = 0;
		lipophilic = YES;
		linkerRing = NO;
		aromatic = NO;
		bondsdone = NO;
		fusedCount = 0;
                dd = ddouts = 0;
		if(p == NULL) number = 0;
		else number = p->getNumber() + 1;
	}
};

double	giveDistance(RawAtom * a1, RawAtom * a2) { return giveDistance(a1->getCoord(), a2->getCoord()); }


// function, that returns, how many neighbors of certain element an atom has

template <class AtomType>
unsigned char	giveNeighborSum(AtomType * a, unsigned char protonNum) {
	unsigned char	i;
	unsigned char	sum = 0;
	
	for(i = 0; i < a->getNeighborCount(); i++) if(a->getNeighbor(i)->getProtonNum() == protonNum) sum++;
	return sum;
}

template <class AtomType>
unsigned char	giveHalogenSum(AtomType * a) {
	unsigned char	i;
	unsigned char	sum = 0;

	for(i = 0; i < a->getNeighborCount(); i++)
		if(a->getNeighbor(i)->getProtonNum() == FLUORINE  ||
		   a->getNeighbor(i)->getProtonNum() == CHLORINE ||
		   a->getNeighbor(i)->getProtonNum() == BROMINE  ||
		   a->getNeighbor(i)->getProtonNum() == IODINE) sum++;
	return sum;
}

class	Centroid {
	private:
	Point		cPoint;
	Centroid *	previous;
	Point		ringNormal;
	Ring *		ring;
	
	public:
	void	setCPoint(Point pt) { cPoint = pt; }
	void	setNormal(Point pt) { ringNormal = pt; }
	void	setPrevious(Centroid * value) { previous = value; }
	
	Point		getCPoint() { return cPoint; }
	Centroid *	getPrevious() { return previous; }
	Point		getNormal() { return ringNormal; }
	Ring *		getRing() { return ring; }
	
	Centroid(Point p, Point rn, Ring * r, Centroid * c) { cPoint = p; ringNormal = rn; ring = r; previous = c; }
};

class   Pharmacophore;

class	Molecule : public StructuralUnit  {
	private:
	LigandAtom *	lastAtom;
	LigandAtom *	help;
	Centroid *	centroid;
	Molecule *	previous;
	Ring *		ring;
	char		totalQ;
        unsigned char   flexible;
        int             rotBonds;
        double          selfIndex;
        double          MW;
        int             n_carbohyd_atoms;
        int             n_atoms;
        int             n_atoms_in_template;
        
        unsigned char   proteinBound;
        unsigned char   metalloligand;
        unsigned char   property;       // heme, modified aa, etc.
        unsigned char   carbohydrate;
        unsigned char   useful;
        unsigned char   wBridge;
        unsigned char   wClose;
        unsigned char   wBound;
        unsigned char   missAtom;
        
        char		dbrefStr[MAX_DBREF_LENGTH];
        char            fullResname[MAX_RES_NAME_LENGTH];
        char            molRootName[128];

	Pharmacophore *	pharmacophore;
        
	public:
        void            addAtom(int num, char * at, char aL, Point p, double occ, double Bf, unsigned char protNum, double Q){
                                help = lastAtom;
                                lastAtom = new LigandAtom(num, at, aL, p, occ, Bf, protNum, Q, lastAtom);
                                lastAtom->setResidue(this);

                                if(help == NULL) firstAtom = (ProteinAtom *) lastAtom;
                                else help->setNext(lastAtom);
                                
                                if(protNum == OXYGEN  ||  protNum == OXYGEN  ||  protNum == SULPHUR) lastAtom->setDAC(DON_ACC_UNSURE);
                        }
	void	removeAllAtoms() {
			while(lastAtom){
				help = lastAtom;
				lastAtom = lastAtom->getPrevious();
				delete help;
			}
		}
	void	setCentroid(Centroid * value) { centroid = value; }
	void	setRing(Ring * value) { ring = value; }
	void	setPrevious(Molecule * value) { previous = value; }
	void	setTotalQ(char value) { totalQ = value; }
	void    setSelfIndex(double value) { selfIndex = value; }
        void    setRotBonds(int value) { rotBonds = value; }
        void    setProperty(unsigned char value) { property = value; }
        void    setProteinBound(unsigned char value) { proteinBound = value; }
        void    setMetalloligand(unsigned char value) { metalloligand = value; }
        void    setCarbohydrate() { carbohydrate = YES; }
        void    setUseful(unsigned char value) { useful = value; }
        void    setCloseWatersCount(unsigned char value) { wClose = value; }
        void    setBridgeWatersCount(unsigned char value) { wBridge = value; }
        void    setBoundWatersCount(unsigned char value) { wBound = value; }
        void    setAtomNumber(int value) { n_atoms = value; }
        void    setAtomNumberInTemplate(int value) { n_atoms_in_template = value; }
        void    setCBHAtomNumber(int value) { n_carbohyd_atoms = value; }
        void    setMissAtom(unsigned char value) { missAtom = value; }
        void    setDBrefStr(char * v) { strcpy(dbrefStr, v); }
        void    setFullResname(char * v) { strcpy(fullResname, v); }
        void    setResname(char * v) { strcpy(residueName, v); }
        void    setMolRootName(char * v) { strcpy(molRootName, v); }
        
        void    addResname(char * v) { 
                        strcat(fullResname, v);
                        if(strlen(fullResname) >= sizeof(fullResname)) {
                                printf("Maximum string length (%d) for the full residue name exceeded\n", MAX_RES_NAME_LENGTH);
                                exit(1);
                        }
                }
        
        void	incTotalQ() { totalQ++; }
        void	decTotalQ() { totalQ--; }
        void    addFlexible() { flexible++; }
        void    addWeight(double w) { MW = MW + w; }

	LigandAtom *	getLastAtom() { return lastAtom; }
	Centroid *	getCentroid() { return centroid; }
	Ring *		getRing() { return ring; }
	Molecule *	getPrevious() { return previous; }
	char		getTotalQ() { return totalQ; }
        LigandAtom *	getFirstAtom() { return (LigandAtom *) firstAtom; }
        unsigned char   getFlexible() { return flexible; }
        double          getSelfIndex() { return selfIndex; }
        int             getRotBonds() { return rotBonds; }
        unsigned char   getProperty() { return property; }
        unsigned char   getProteinBound() { return proteinBound; }
        unsigned char   getMetalloligand() { return metalloligand; }
        unsigned char   getCarbohydrate() { return carbohydrate; }
        unsigned char   getUseful() { return useful; }
        unsigned char   getCloseWatersCount() { return wClose; }
        unsigned char   getBridgeWatersCount() { return wBridge; }
        unsigned char   getBoundWatersCount() { return wBound; }
        int             getAtomNumber() { return n_atoms; }
        int             getAtomNumberInTemplate() { return n_atoms_in_template; }
        int             getCBHAtomNumber() { return n_carbohyd_atoms; }
        double          getMolecularWeight() { return MW; }
        char *		getDBrefStr() { return dbrefStr; }
        unsigned char   getMissAtom() { return missAtom; }
        char *		getFullResname() { return fullResname; }
        char *		getMolRootName() { return molRootName; }
	Pharmacophore *	getPharmacophore() { return pharmacophore; }
        
        Molecule(int n, char l, char * rname, char ch, Molecule * p) {
		unitNumber = n;
		lastAtom = NULL;
		strcpy(residueName, rname);
                residueName[3] = 0;
		previous = p;
                ID = l;
		totalQ = 0;
                flexible = 0;
                selfIndex = rotBonds = 0;
                property = 0;
                proteinBound = NO;
                chain = ch;
                metalloligand = NO;
                carbohydrate = NO;
                useful = YES;
                wBridge = wClose = wBound = 0;
                MW = n_carbohyd_atoms = n_atoms = 0;
                dbrefStr[0] = 0;
                missAtom = 0;
                n_atoms_in_template = 0;
                fullResname[0] = 0;
                strcpy(molRootName, "noname.dat");
	}
};

// coefficients taken from ATOMTYP from schrodinger
double	giveRadiusSpecialized(ProteinAtom * a) {
        unsigned char   flag;     // multipurpose variable
	switch (a->getProtonNum()) {
		case   HYDROGEN   :
                        flag = a->getNeighbor(0)->getProtonNum();
                             if(flag == CARBON)   return 1.09;
                        else if(flag == NITROGEN) return 1.05;
                        else if(flag == OXYGEN)   return 1.00;
                        else if(flag == SULPHUR)  return 1.09;
                        else { printf("Problem with VdW radius -> atom %s %d\n", a->getAtomType(), a->getNumber()); exit(1); }
                        break;
                        
		case   CARBON     :
                        // if not heme, apply united atom rules
                        if(a->getResidue()->getID() != 'X') {
                                // get H count
                                flag = a->getHcount();
                                // SP3
                                if(a->getHybridization() == SP3) {
                                             if(flag == 0) return 1.7;
                                        else if(flag == 1) return 1.75;
                                        else if(flag == 2) return 1.8;
                                        else if(flag == 3) return 1.8;
                                        else { printf("Problem with VdW radius -> atom %s %d\n", a->getAtomType(), a->getNumber()); exit(1); }
                                }
                                // SP2
                                else if(a->getHybridization() == SP2) {
                                             if(flag == 0) return 1.72;
                                        else if(flag == 1) return 1.75;
                                        else if(flag == 2) return 1.8;
                                        else { printf("Problem with VdW radius -> atom %s %d\n", a->getAtomType(), a->getNumber()); exit(1); }
                                }
                                // SP1
                                else if(a->getHybridization() == SP) {
                                             if(flag == 0) return 1.78;
                                        else if(flag == 1) return 1.78;
                                        else { printf("Problem with VdW radius -> atom %s %d\n", a->getAtomType(), a->getNumber()); exit(1); }     
                                }
                                else { printf("Problem with VdW radius -> atom %s %d\n", a->getAtomType(), a->getNumber()); exit(1); }
                        }
                        else {
                                flag = a->getNeighborCount();
                                // SP3
                                     if(flag == 4) return 1.70;
                                // SP2
                                else if(flag == 3) return 1.72;
                                // SP1
                                else if(flag == 2) return 1.78;
                        }
        
                        break;
                        
		case   NITROGEN   :
                             if (a->getHybridization() == SP3) return 1.6;
                        else if (a->getHybridization() == SP2) return 1.55;
                        else if (a->getHybridization() == SP)  return 1.55;
                        else { printf("Problem with VdW radius -> atom %s %d\n", a->getAtomType(), a->getNumber()); exit(1); }
                        break;
                        
		case   OXYGEN     :
                        
                             if (a->getHybridization() == SP3) return 1.52;
                        else if (a->getHybridization() == SP2) return 1.5;
                        else if (strcmp(a->getAtomType(), "O2-") == 0) return 1.5;     
                        else { printf("Problem with VdW radius -> atom %s %d\n", a->getAtomType(), a->getNumber()); exit(1); }     
                        break;
                        
		case   FLUORINE   : return  1.47   ;
		case   SILICON    : return  2.1    ;
		case   PHOSPHORUS : return  1.8    ;
		case   SULPHUR    : return  1.8    ;
		case   CHLORINE   : return  1.75   ;
		case   IRON       : return  2.0    ;
		case   ZINC       : return  1.39   ;
		case   BROMINE    : return  1.85   ;
		case   IODINE     : return  1.98   ;
		default : printf("VdW routine: Unknown atom N=%d\n", a->getProtonNum()); exit(1); break;
	}
	return 0;
}

double	giveRadiusSpecialized(LigandAtom * a) {
        unsigned char   flag;     // multipurpose variable
	switch (a->getProtonNum()) {
		case   HYDROGEN   :
                        flag = a->getNeighbor(0)->getProtonNum();
                             if(flag == CARBON)   return 1.09;
                        else if(flag == NITROGEN) return 1.05;
                        else if(flag == OXYGEN)   return 1.00;
                        else if(flag == SULPHUR)  return 1.09;
                        else { printf("Problem with VdW radius -> atom %s %d\n", a->getAtomType(), a->getNumber()); exit(1); }
                        break;
                        
		case   CARBON     :
                        flag = a->getHybridization();
                        // SP3
                             if(flag == SP3) return 1.70;
                        // SP2
                        else if(flag == SP2) return 1.72;
                        // SP1
                        else if(flag == SP ) return 1.78;
                        else { printf("Problem with VdW radius -> atom %s %d\n", a->getAtomType(), a->getNumber()); exit(1); }
                        break;
                        
		case   NITROGEN   :
                        flag = a->getHybridization();
                             if (flag == SP3) return 1.60;
                        else if (flag == SP2) return 1.55;
                        else if (flag == SP ) return 1.55;
                        else { printf("Problem with VdW radius -> atom %s %d\n", a->getAtomType(), a->getNumber()); exit(1); }
                        break;
                        
		case   OXYGEN     :
                        flag = a->getHybridization();
                             if (flag == SP3) return 1.52;
                        else if (flag == SP2) return 1.50;
                        else { printf("Problem with VdW radius -> atom %s %d\n", a->getAtomType(), a->getNumber()); exit(1); }     
                        break;
                        
		case   SULPHUR    : return  1.8    ;
		case   FLUORINE   : return  1.47   ;
		case   CHLORINE   : return  1.75   ;
		case   BROMINE    : return  1.85   ;
		case   IODINE     : return  1.98   ;

                case   SILICON    : return  2.1    ;
		case   PHOSPHORUS : return  1.8    ;
		default : printf("VdW routine: Unknown atom N=%d\n", a->getProtonNum()); exit(1); break;
	}
	return 0;
}

char    atom_giveDonAcc(ProteinAtom * pa2) {
        char    pa[5];
        sscanf(pa2->getAtomType(), "%s", pa);
        pa[4] = 0;
        if (pa2->getProtonNum() != CARBON) {
                // backbone
                if(strcmp(pa, "N") == 0) return DONOR;
                if(strcmp(pa, "O") == 0) return ACCEPTOR;
                if(strcmp(pa, "OXT") == 0) return ACCEPTOR;
                // GLN - GLU
                if(strcmp(pa, "OE1") == 0) return ACCEPTOR;
                if(strcmp(pa, "OE2") == 0) return ACCEPTOR;
                if(strcmp(pa, "NE2") == 0) return DONOR;
                // ASN - ASP
                if(strcmp(pa, "OD1") == 0) return ACCEPTOR;
                if(strcmp(pa, "OD2") == 0) return ACCEPTOR;
                if(strcmp(pa, "ND2") == 0) return DONOR;
                // THR, SER, CYS, TYR
                if(strcmp(pa, "OG1") == 0) return DON_ACC;
                if(strcmp(pa, "OG") == 0) return DON_ACC;
                if(strcmp(pa, "SG") == 0) return DON_ACC;
                if(strcmp(pa, "OH") == 0) return DON_ACC;
                // MET
                if(strcmp(pa, "SD") == 0) return ACCEPTOR;
                // TRP, LYS, ARG
                if(strcmp(pa, "NE1") == 0) return DONOR;
                if(strcmp(pa, "NZ") == 0) return DONOR;
                if(strcmp(pa, "NE") == 0) return DONOR;
                if(strcmp(pa, "NH1") == 0) return DONOR;
                if(strcmp(pa, "NH2") == 0) return DONOR;
                // HIS - always unsure due to protonation
                if(strcmp(pa, "ND1") == 0) return DON_ACC_UNSURE;
                if(strcmp(pa, "NE2") == 0) return DON_ACC_UNSURE;
        }
        return NO;
}

#endif	/* STRUCTURES_H */

