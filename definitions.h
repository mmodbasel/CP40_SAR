/* 
 * File:   definitions.h
 * Author: Martin Smiesko, University of Basel
 *
 * Created on June 19, 2014, 4:29 PM
 */

#ifndef DEFINITIONS_H
#define	DEFINITIONS_H

#define	VERSION			1.00
#define	YES			1
#define NO			0

#define	MAX_NEIGHBOR_LIST_SIZE 	64
#define MAX_RING_SIZE 		64
#define MAX_DBREF_LENGTH        32
#define MAX_RES_NAME_LENGTH    128     // length of the combined residue name from all subunits XXX-A-9999
#define MIN_RESIDUES_AT_BINDINGSITE     10
#define MAXCHAINS               26      // alphabet A...Z
#define MIN_CHAIN_LENGTH        50

#define	SP4			4	// middle atom of azide (-N=N=N), cyanate (-N=C=O) or propadiene (-C=C=CH2)
#define	SP3			3
#define	SP2			2
#define	SP			1
#define	NONE			0
#define	UNCLEAR			-1

#define	SINGLE			1
#define	DOUBLE			2
#define TRIPLE			3
#define	DELOC			4

#define	DONOR                   1
#define	ACCEPTOR		2
#define DON_ACC                 3
#define DON_ACC_UNSURE          4

#define CCDC_CUTOFF             0.4
#define HB_THRESHOLD_ANGLE       RAD_105
#define HB_MAX_DIST             2.5                     // maximum distance H...Y  (in X-H...Y)
#define HB_MIN_ENERGY           -0.5                     // maximum distance H...Y  (in X-H...Y)
#define LIPO_MIN_ENERGY         -0.5                     // maximum distance H...Y  (in X-H...Y)
#define EHPHOB              -0.3
#define CUTOFF_NonBonded    12.0
#define ACTIVESITE_RADIUS       6.0

#define BUMP_CHECK_FACTOR_1_4	0.632                   //  with new VdW radii 0.795 * 0.795
#define BUMP_CHECK_FACTOR_AA	0.716477603             // 5% shorter
#define BUMP_CHECK_FACTOR_HB	0.325			// less stringent HB factor 1.50 / 2.52 = 0.59524 ---> 0.3543
#define BUMP_CHECK_FACTOR_DA	0.65                    // second order factor i.e. for O....O in O-H..O  2.5 / 3.04 = 0.8224 ----> 0.67629

#define FAMILY_MAX_RESIDUES     128
// #define ATOM_MAX_NEIGHBORS      4
#define ATOM_MAX_NEIGHBORS      6
#define METAL_MAX_NEIGHBORS     8

#define ATOM_LINE               1
#define HETATM_LINE             2

#define MW_LOW_LIMIT    100.0
#define MW_HIGH_LIMIT   800.0
#define MIN_CBH_FRACTION        0.333
#define ResNoLimit      20
#define RESNUM_LOW_LIMIT        6
#define RESNUM_HIGH_LIMIT       50
#define MIN_ATOMS       5
#define MAX_ATOMS       250

#define PROTEIN         0
#define WATER           1
#define LIGAND          2
#define METAL           3
#define METALLOLIGAND   4
#define COVALENT_LIGAND 5
#define COFACTOR        6       // e.g. heme -> needed to be include as a protein for the surface areas

#define WATER_CLOSE     1
#define WATER_BOUND     2
#define WATER_BRIDGE    3

#define OUT_OF_PLANE_LIMIT    PI * 15.0 / 180.0
#define LINEARITY_LIMIT         RAD_165
#define WATER_HB_MAX_DIST       3.5
#define METAL_LIGAND_MAX_DIST   3.0

#define 	ACTINIUM	89
#define 	ALUMINIUM	13
#define 	AMERICIUM	95
#define 	ANTIMONY	51
#define 	ARGON           18
#define 	ARSENIC         33
#define 	ASTATINE	85
#define 	BARIUM          56
#define 	BERKELIUM	97
#define 	BERYLLIUM	4
#define 	BISMUTH         83
#define 	BOHRIUM         107
#define 	BORON           5
#define 	BROMINE         35
#define 	CADMIUM         48
#define 	CAESIUM         55
#define 	CALCIUM         20
#define 	CALIFORNIUM	98
#define 	CARBON          6
#define 	CERIUM          58
#define 	CHLORINE	17
#define 	CHROMIUM	24
#define 	COBALT          27
#define 	COPPER          29
#define 	CURIUM          96
#define 	DARMSTADTIUM	110
#define 	DUBNIUM         105
#define 	DYSPROSIUM	66
#define 	EINSTEINIUM	99
#define 	ERBIUM          68
#define 	EUROPIUM	63
#define 	FERMIUM         100
#define 	FLUORINE	9
#define 	FRANCIUM	87
#define 	GADOLINIUM	64
#define 	GALLIUM         31
#define 	GERMANIUM	32
#define 	GOLD            79
#define 	HAFNIUM         72
#define 	HASSIUM         108
#define 	HELIUM          2
#define 	HOLMIUM         67
#define 	HYDROGEN	1
#define 	INDIUM          49
#define 	IODINE          53
#define 	IRIDIUM         77
#define 	IRON            26
#define 	KRYPTON         36
#define 	LANTHANUM	57
#define 	LAWRENCIUM	103
#define 	LEAD            82
#define 	LITHIUM         3
#define 	LUTETIUM	71
#define 	MAGNESIUM	12
#define 	MANGANESE	25
#define 	MEITNERIUM	109
#define 	MENDELEVIUM	101
#define 	MERCURY         80
#define 	MOLYBDENUM	42
#define 	NEODYMIUM	60
#define 	NEON            10
#define 	NEPTUNIUM	93
#define 	NICKEL          28
#define 	NIOBIUM         41
#define 	NITROGEN	7
#define 	NOBELIUM	102
#define 	OSMIUM          76
#define 	OXYGEN          8
#define 	PALLADIUM	46
#define 	PHOSPHORUS	15
#define 	PLATINUM	78
#define 	PLUTONIUM	94
#define 	POLONIUM	84
#define 	POTASSIUM	19
#define 	PRASEODYMIUM	59
#define 	PROMETHIUM	61
#define 	PROTACTINIUM	91
#define 	RADIUM          88
#define 	RADON           86
#define 	RHENIUM         75
#define 	RHODIUM         45
#define 	RUBIDIUM	37
#define 	RUTHENIUM	44
#define 	RUTHERFORDIUM	104
#define 	SAMARIUM	62
#define 	SCANDIUM	21
#define 	SEABORGIUM	106
#define 	SELENIUM	34
#define 	SILICON         14
#define 	SILVER          47
#define 	SODIUM          11
#define 	STRONTIUM	38
#define 	SULPHUR         16
#define 	TANTALUM	73
#define 	TECHNETIUM	43
#define 	TELLURIUM	52
#define 	TERBIUM         65
#define 	THALLIUM	81
#define 	THORIUM         90
#define 	THULIUM         69
#define 	TIN             50
#define 	TITANIUM	22
#define 	TUNGSTEN	74
#define 	URANIUM         92
#define 	VANADIUM	23
#define 	XENON           54
#define 	YTTERBIUM	70
#define 	YTTRIUM         39
#define 	ZINC            30
#define 	ZIRCONIUM	40
#define         DEUTERIUM       1

#define N_AMINOACID     22
#define GLY     0
#define ALA     1
#define VAL     2
#define LEU     3
#define ILE     4
#define PRO     5
#define PHE     6
#define SER     7
#define THR     8
#define TYR     9
#define TRP     10
#define HIS     11
#define MET     12
#define CYS     13
#define ASP     14
#define ASN     15
#define GLU     16
#define GLN     17
#define ARG     18
#define LYS     19
// selenocysteine
#define SEC     20
#define MSE     21

int     residueName_getIndex(char * resname) {
        if(strcmp(resname, "GLY") == 0) return GLY;
        if(strcmp(resname, "ALA") == 0) return ALA;
        if(strcmp(resname, "VAL") == 0) return VAL;
        if(strcmp(resname, "LEU") == 0) return LEU;
        if(strcmp(resname, "ILE") == 0) return ILE;
        if(strcmp(resname, "PRO") == 0) return PRO;
        if(strcmp(resname, "PHE") == 0) return PHE;
        if(strcmp(resname, "SER") == 0) return SER;
        if(strcmp(resname, "THR") == 0) return THR;
        if(strcmp(resname, "TYR") == 0) return TYR;
        if(strcmp(resname, "TRP") == 0) return TRP;
        if(strcmp(resname, "HIS") == 0) return HIS;
        if(strcmp(resname, "MET") == 0) return MET;
        if(strcmp(resname, "CYS") == 0) return CYS;
        if(strcmp(resname, "ASP") == 0) return ASP;
        if(strcmp(resname, "ASN") == 0) return ASN;
        if(strcmp(resname, "GLU") == 0) return GLU;
        if(strcmp(resname, "GLN") == 0) return GLN;
        if(strcmp(resname, "ARG") == 0) return ARG;
        if(strcmp(resname, "LYS") == 0) return LYS;
        // selenocysteine
        if(strcmp(resname, "SEC") == 0) return SEC;
        if(strcmp(resname, "MSE") == 0) return MSE;
        
        printf("Error: No index for residue %s\n", resname);
        exit(1);
}

void    residueIndex_getName(int i, char * resname) {
             if(i == GLY) strcpy(resname, "GLY");
        else if(i == ALA) strcpy(resname, "ALA");
        else if(i == VAL) strcpy(resname, "VAL");
        else if(i == LEU) strcpy(resname, "LEU");
        else if(i == ILE) strcpy(resname, "ILE");
        else if(i == PRO) strcpy(resname, "PRO");
        else if(i == PHE) strcpy(resname, "PHE");
        else if(i == SER) strcpy(resname, "SER");
        else if(i == THR) strcpy(resname, "THR");
        else if(i == TYR) strcpy(resname, "TYR");
        else if(i == TRP) strcpy(resname, "TRP");
        else if(i == HIS) strcpy(resname, "HIS");
        else if(i == MET) strcpy(resname, "MET");
        else if(i == CYS) strcpy(resname, "CYS");
        else if(i == ASP) strcpy(resname, "ASP");
        else if(i == ASN) strcpy(resname, "ASN");
        else if(i == GLU) strcpy(resname, "GLU");
        else if(i == GLN) strcpy(resname, "GLN");
        else if(i == ARG) strcpy(resname, "ARG");
        else if(i == LYS) strcpy(resname, "LYS");
        // selenocysteine
        else if(i == SEC) strcpy(resname, "SEC");
        else if(i == MSE) strcpy(resname, "MSE");
        else {
                printf("Error: No residue name for index %d\n", i);
                exit(1);
        }
}

int	giveProtonNumber(char * el) {
	     if(strcmp(el, " H") == 0) return HYDROGEN;
        else if(strcmp(el, " D") == 0) return DEUTERIUM;
	else if(strcmp(el, " C") == 0) return CARBON;
	else if(strcmp(el, " N") == 0) return NITROGEN;
	else if(strcmp(el, " O") == 0) return OXYGEN;
	else if(strcmp(el, " S") == 0) return SULPHUR;
	else if(strcmp(el, " P") == 0) return PHOSPHORUS;
	else if(strcmp(el, " F") == 0) return FLUORINE;
	else if(strcmp(el, " I") == 0) return IODINE;
	else if(strcmp(el, "CL") == 0) return CHLORINE;
	else if(strcmp(el, "BR") == 0) return BROMINE;
        else if(strcmp(el, "ZN") == 0) return ZINC;
        else if(strcmp(el, "FE") == 0) return IRON;
        else if(strcmp(el, "AC") == 0) return ACTINIUM;
        else if(strcmp(el, "AL") == 0) return ALUMINIUM;
        else if(strcmp(el, "AM") == 0) return AMERICIUM;
        else if(strcmp(el, "SB") == 0) return ANTIMONY;
        else if(strcmp(el, "AR") == 0) return ARGON;
        else if(strcmp(el, "AS") == 0) return ARSENIC;
        else if(strcmp(el, "AT") == 0) return ASTATINE;
        else if(strcmp(el, "BA") == 0) return BARIUM;
        else if(strcmp(el, "BK") == 0) return BERKELIUM;
        else if(strcmp(el, "BE") == 0) return BERYLLIUM;
        else if(strcmp(el, "BI") == 0) return BISMUTH;
        else if(strcmp(el, "BH") == 0) return BOHRIUM;
        else if(strcmp(el, " B") == 0) return BORON;
        else if(strcmp(el, "CD") == 0) return CADMIUM;
        else if(strcmp(el, "CS") == 0) return CAESIUM;
        else if(strcmp(el, "CA") == 0) return CALCIUM;
        else if(strcmp(el, "CF") == 0) return CALIFORNIUM;
        else if(strcmp(el, "CE") == 0) return CERIUM;
        else if(strcmp(el, "CR") == 0) return CHROMIUM;
        else if(strcmp(el, "CO") == 0) return COBALT;
        else if(strcmp(el, "CU") == 0) return COPPER;
        else if(strcmp(el, "CM") == 0) return CURIUM;
        else if(strcmp(el, "DS") == 0) return DARMSTADTIUM;
        else if(strcmp(el, "DB") == 0) return DUBNIUM;
        else if(strcmp(el, "DY") == 0) return DYSPROSIUM;
        else if(strcmp(el, "ES") == 0) return EINSTEINIUM;
        else if(strcmp(el, "ER") == 0) return ERBIUM;
        else if(strcmp(el, "EU") == 0) return EUROPIUM;
        else if(strcmp(el, "FM") == 0) return FERMIUM;
        else if(strcmp(el, "FR") == 0) return FRANCIUM;
        else if(strcmp(el, "GD") == 0) return GADOLINIUM;
        else if(strcmp(el, "GA") == 0) return GALLIUM;
        else if(strcmp(el, "GE") == 0) return GERMANIUM;
        else if(strcmp(el, "AU") == 0) return GOLD;
        else if(strcmp(el, "HF") == 0) return HAFNIUM;
        else if(strcmp(el, "HS") == 0) return HASSIUM;
        else if(strcmp(el, "HE") == 0) return HELIUM;
        else if(strcmp(el, "HO") == 0) return HOLMIUM;
        else if(strcmp(el, "IN") == 0) return INDIUM;
        else if(strcmp(el, "IR") == 0) return IRIDIUM;
        else if(strcmp(el, "KR") == 0) return KRYPTON;
        else if(strcmp(el, "LA") == 0) return LANTHANUM;
        else if(strcmp(el, "LR") == 0) return LAWRENCIUM;
        else if(strcmp(el, "PB") == 0) return LEAD;
        else if(strcmp(el, "LI") == 0) return LITHIUM;
        else if(strcmp(el, "LU") == 0) return LUTETIUM;
        else if(strcmp(el, "MG") == 0) return MAGNESIUM;
        else if(strcmp(el, "MN") == 0) return MANGANESE;
        else if(strcmp(el, "MT") == 0) return MEITNERIUM;
        else if(strcmp(el, "MD") == 0) return MENDELEVIUM;
        else if(strcmp(el, "HG") == 0) return MERCURY;
        else if(strcmp(el, "MO") == 0) return MOLYBDENUM;
        else if(strcmp(el, "ND") == 0) return NEODYMIUM;
        else if(strcmp(el, "NE") == 0) return NEON;
        else if(strcmp(el, "NP") == 0) return NEPTUNIUM;
        else if(strcmp(el, "NI") == 0) return NICKEL;
        else if(strcmp(el, "NB") == 0) return NIOBIUM;
        else if(strcmp(el, "NO") == 0) return NOBELIUM;
        else if(strcmp(el, "OS") == 0) return OSMIUM;
        else if(strcmp(el, "PD") == 0) return PALLADIUM;
        else if(strcmp(el, "PT") == 0) return PLATINUM;
        else if(strcmp(el, "PU") == 0) return PLUTONIUM;
        else if(strcmp(el, "PO") == 0) return POLONIUM;
        else if(strcmp(el, " K") == 0) return POTASSIUM;
        else if(strcmp(el, "PR") == 0) return PRASEODYMIUM;
        else if(strcmp(el, "PM") == 0) return PROMETHIUM;
        else if(strcmp(el, "PA") == 0) return PROTACTINIUM;
        else if(strcmp(el, "RA") == 0) return RADIUM;
        else if(strcmp(el, "RN") == 0) return RADON;
        else if(strcmp(el, "RE") == 0) return RHENIUM;
        else if(strcmp(el, "RH") == 0) return RHODIUM;
        else if(strcmp(el, "RB") == 0) return RUBIDIUM;
        else if(strcmp(el, "RU") == 0) return RUTHENIUM;
        else if(strcmp(el, "RF") == 0) return RUTHERFORDIUM;
        else if(strcmp(el, "SM") == 0) return SAMARIUM;
        else if(strcmp(el, "SC") == 0) return SCANDIUM;
        else if(strcmp(el, "SG") == 0) return SEABORGIUM;
        else if(strcmp(el, "SE") == 0) return SELENIUM;
        else if(strcmp(el, "SI") == 0) return SILICON;
        else if(strcmp(el, "AG") == 0) return SILVER;
        else if(strcmp(el, "NA") == 0) return SODIUM;
        else if(strcmp(el, "SR") == 0) return STRONTIUM;
        else if(strcmp(el, "TA") == 0) return TANTALUM;
        else if(strcmp(el, "TC") == 0) return TECHNETIUM;
        else if(strcmp(el, "TE") == 0) return TELLURIUM;
        else if(strcmp(el, "TB") == 0) return TERBIUM;
        else if(strcmp(el, "TL") == 0) return THALLIUM;
        else if(strcmp(el, "TH") == 0) return THORIUM;
        else if(strcmp(el, "TM") == 0) return THULIUM;
        else if(strcmp(el, "SN") == 0) return TIN;
        else if(strcmp(el, "TI") == 0) return TITANIUM;
        else if(strcmp(el, " W") == 0) return TUNGSTEN;
        else if(strcmp(el, " U") == 0) return URANIUM;
        else if(strcmp(el, " V") == 0) return VANADIUM;
        else if(strcmp(el, "XE") == 0) return XENON;
        else if(strcmp(el, "YB") == 0) return YTTERBIUM;
        else if(strcmp(el, " Y") == 0) return YTTRIUM;
        else if(strcmp(el, "ZR") == 0) return ZIRCONIUM;
	else {
                printf("Cannot find proton number for atom %s\n", el);
                exit(1);
	}
}

double	giveVdWRadius(unsigned char N) {
        switch(N) {
            case	ACTINIUM	: return 2.00;	break;
            case	ALUMINIUM	: return 2.00;	break;
            case	AMERICIUM	: return 2.00;	break;
            case	ANTIMONY	: return 2.00;	break;
            case	ARGON		: return 1.88;	break;
            case	ARSENIC		: return 1.85;	break;
            case	ASTATINE	: return 2.00;	break;
            case	BARIUM		: return 2.00;	break;
            case	BERKELIUM	: return 2.00;	break;
            case	BERYLLIUM	: return 2.00;	break;
            case	BISMUTH		: return 2.00;	break;
            case	BOHRIUM		: return 2.00;	break;
            case	BORON		: return 2.00;	break;
            case	BROMINE		: return 1.85;	break;
            case	CADMIUM		: return 1.58;	break;
            case	CAESIUM		: return 2.00;	break;
            case	CALCIUM		: return 2.00;	break;
            case	CALIFORNIUM	: return 2.00;	break;
            case	CARBON		: return 1.70;	break;
            case	CERIUM		: return 2.00;	break;
            case	CHLORINE	: return 1.75;	break;
            case	CHROMIUM	: return 2.00;	break;
            case	COBALT		: return 2.00;	break;
            case	COPPER		: return 1.40;	break;
            case	CURIUM		: return 2.00;	break;
            case	DARMSTADTIUM	: return 2.00;	break;
            case	DUBNIUM		: return 2.00;	break;
            case	DYSPROSIUM	: return 2.00;	break;
            case	EINSTEINIUM	: return 2.00;	break;
            case	ERBIUM		: return 2.00;	break;
            case	EUROPIUM	: return 2.00;	break;
            case	FERMIUM		: return 2.00;	break;
            case	FLUORINE	: return 1.47;	break;
            case	FRANCIUM	: return 2.00;	break;
            case	GADOLINIUM	: return 2.00;	break;
            case	GALLIUM		: return 1.87;	break;
            case	GERMANIUM	: return 2.00;	break;
            case	GOLD		: return 1.66;	break;
            case	HAFNIUM		: return 2.00;	break;
            case	HASSIUM		: return 2.00;	break;
            case	HELIUM		: return 1.40;	break;
            case	HOLMIUM		: return 2.00;	break;
            case	HYDROGEN	: return 1.09;	break;
            case	INDIUM		: return 1.93;	break;
            case	IODINE		: return 1.98;	break;
            case	IRIDIUM		: return 2.00;	break;
            case	IRON		: return 2.00;	break;
            case	KRYPTON		: return 2.02;	break;
            case	LANTHANUM	: return 2.00;	break;
            case	LAWRENCIUM	: return 2.00;	break;
            case	LEAD		: return 2.02;	break;
            case	LITHIUM		: return 1.82;	break;
            case	LUTETIUM	: return 2.00;	break;
            case	MAGNESIUM	: return 1.73;	break;
            case	MANGANESE	: return 2.00;	break;
            case	MEITNERIUM	: return 2.00;	break;
            case	MENDELEVIUM	: return 2.00;	break;
            case	MERCURY		: return 1.55;	break;
            case	MOLYBDENUM	: return 2.00;	break;
            case	NEODYMIUM	: return 2.00;	break;
            case	NEON		: return 1.54;	break;
            case	NEPTUNIUM	: return 2.00;	break;
            case	NICKEL		: return 1.63;	break;
            case	NIOBIUM		: return 2.00;	break;
            case	NITROGEN	: return 1.55;	break;
            case	NOBELIUM	: return 2.00;	break;
            case	OSMIUM		: return 2.00;	break;
            case	OXYGEN		: return 1.52;	break;
            case	PALLADIUM	: return 1.63;	break;
            case	PHOSPHORUS	: return 1.80;	break;
            case	PLATINUM	: return 1.72;	break;
            case	PLUTONIUM	: return 2.00;	break;
            case	POLONIUM	: return 2.00;	break;
            case	POTASSIUM	: return 2.75;	break;
            case	PRASEODYMIUM	: return 2.00;	break;
            case	PROMETHIUM	: return 2.00;	break;
            case	PROTACTINIUM	: return 2.00;	break;
            case	RADIUM		: return 2.00;	break;
            case	RADON		: return 2.00;	break;
            case	RHENIUM		: return 2.00;	break;
            case	RHODIUM		: return 2.00;	break;
            case	RUBIDIUM	: return 2.00;	break;
            case	RUTHENIUM	: return 2.00;	break;
            case	RUTHERFORDIUM	: return 2.00;	break;
            case	SAMARIUM	: return 2.00;	break;
            case	SCANDIUM	: return 2.00;	break;
            case	SEABORGIUM	: return 2.00;	break;
            case	SELENIUM	: return 1.90;	break;
            case	SILICON		: return 2.10;	break;
            case	SILVER		: return 1.72;	break;
            case	SODIUM		: return 2.27;	break;
            case	STRONTIUM	: return 2.00;	break;
            case	SULPHUR		: return 1.80;	break;
            case	TANTALUM	: return 2.00;	break;
            case	TECHNETIUM	: return 2.00;	break;
            case	TELLURIUM	: return 2.06;	break;
            case	TERBIUM		: return 2.00;	break;
            case	THALLIUM	: return 1.96;	break;
            case	THORIUM		: return 2.00;	break;
            case	THULIUM		: return 2.00;	break;
            case	TIN		: return 2.17;	break;
            case	TITANIUM	: return 2.00;	break;
            case	TUNGSTEN	: return 2.00;	break;
            case	URANIUM		: return 1.86;	break;
            case	VANADIUM	: return 2.00;	break;
            case	XENON		: return 2.16;	break;
            case	YTTERBIUM	: return 2.00;	break;
            case	YTTRIUM		: return 2.00;	break;
            case	ZINC		: return 1.39;	break;
            case	ZIRCONIUM	: return 2.00;	break;
            default: break;
        }
        
        printf("VdW: Unsupported element proton number %d\n", N);
        exit(1);
        
	return 0;
}

char	giveElement(unsigned char N) {
        switch(N) {
            case	HYDROGEN	: return 'H'; break;
            case	CARBON		: return 'C'; break;
            case	NITROGEN        : return 'N'; break;
            case	OXYGEN		: return 'O'; break;
            case	FLUORINE	: return 'F'; break;
            case	PHOSPHORUS	: return 'P'; break;
            case	SULPHUR		: return 'S'; break;
            case	IODINE		: return 'I'; break;
            
            default: break;
       }
        
        printf("Cannot convert element with proton number %d\n", N);
        exit(1);
        
	return 0;
}

double	giveCovalentRadius(unsigned char N) {
        switch(N) {
            case	ACTINIUM	: return 2.15; break;
            case	ALUMINIUM	: return 1.21; break;
            case	AMERICIUM	: return 1.80; break;
            case	ANTIMONY	: return 1.39; break;
            case	ARGON		: return 1.51; break;
            case	ARSENIC		: return 1.21; break;
            case	ASTATINE	: return 1.21; break;
            case	BARIUM		: return 2.15; break;
            case	BERKELIUM	: return 1.54; break;
            case	BERYLLIUM	: return 0.96; break;
            case	BISMUTH		: return 1.48; break;
            case	BOHRIUM		: return 1.50; break;
            case	BORON		: return 0.83; break;
            case	BROMINE		: return 1.21; break;
            case	CADMIUM		: return 1.54; break;
            case	CAESIUM		: return 2.44; break;
            case	CALCIUM		: return 1.76; break;
            case	CALIFORNIUM	: return 1.83; break;
            case	CARBON		: return 0.68; break;
            case	CERIUM		: return 2.04; break;
            case	CHLORINE	: return 0.99; break;
            case	CHROMIUM	: return 1.39; break;
            case	COBALT		: return 1.26; break;
            case	COPPER		: return 1.32; break;
            case	CURIUM		: return 1.69; break;
            case	DARMSTADTIUM	: return 1.50; break;
            case	DUBNIUM		: return 1.50; break;
            case	DYSPROSIUM	: return 1.92; break;
            case	EINSTEINIUM	: return 1.50; break;
            case	ERBIUM		: return 1.89; break;
            case	EUROPIUM	: return 1.98; break;
            case	FERMIUM		: return 1.50; break;
            case	FLUORINE	: return 0.64; break;
            case	FRANCIUM	: return 2.60; break;
            case	GADOLINIUM	: return 1.96; break;
            case	GALLIUM		: return 1.22; break;
            case	GERMANIUM	: return 1.17; break;
            case	GOLD		: return 1.36; break;
            case	HAFNIUM		: return 1.75; break;
            case	HASSIUM		: return 1.50; break;
            case	HELIUM		: return 1.50; break;
            case	HOLMIUM		: return 1.92; break;
            case	HYDROGEN	: return 0.23; break;
            case	INDIUM		: return 1.42; break;
            case	IODINE		: return 1.40; break;
            case	IRIDIUM		: return 1.41; break;
            case	IRON		: return 1.52; break;
            case	KRYPTON		: return 1.50; break;
            case	LANTHANUM	: return 2.07; break;
            case	LAWRENCIUM	: return 1.50; break;
            case	LEAD		: return 1.46; break;
            case	LITHIUM		: return 1.28; break;
            case	LUTETIUM	: return 1.87; break;
            case	MAGNESIUM	: return 1.41; break;
            case	MANGANESE	: return 1.61; break;
            case	MEITNERIUM	: return 1.50; break;
            case	MENDELEVIUM	: return 1.50; break;
            case	MERCURY		: return 1.32; break;
            case	MOLYBDENUM	: return 1.54; break;
            case	NEODYMIUM	: return 2.01; break;
            case	NEON		: return 1.50; break;
            case	NEPTUNIUM	: return 1.90; break;
            case	NICKEL		: return 1.24; break;
            case	NIOBIUM		: return 1.64; break;
            case	NITROGEN	: return 0.68; break;
            case	NOBELIUM	: return 1.50; break;
            case	OSMIUM		: return 1.44; break;
            case	OXYGEN		: return 0.68; break;
            case	PALLADIUM	: return 1.39; break;
            case	PHOSPHORUS	: return 1.05; break;
            case	PLATINUM	: return 1.36; break;
            case	PLUTONIUM	: return 1.87; break;
            case	POLONIUM	: return 1.40; break;
            case	POTASSIUM	: return 2.03; break;
            case	PRASEODYMIUM	: return 2.03; break;
            case	PROMETHIUM	: return 1.99; break;
            case	PROTACTINIUM	: return 2.00; break;
            case	RADIUM		: return 2.21; break;
            case	RADON		: return 1.50; break;
            case	RHENIUM		: return 1.51; break;
            case	RHODIUM		: return 1.42; break;
            case	RUBIDIUM	: return 2.20; break;
            case	RUTHENIUM	: return 1.46; break;
            case	RUTHERFORDIUM	: return 1.50; break;
            case	SAMARIUM	: return 1.98; break;
            case	SCANDIUM	: return 1.70; break;
            case	SEABORGIUM	: return 1.50; break;
            case	SELENIUM	: return 1.22; break;
            case	SILICON		: return 1.20; break;
            case	SILVER		: return 1.45; break;
            case	SODIUM		: return 1.66; break;
            case	STRONTIUM	: return 1.95; break;
            case	SULPHUR		: return 1.02; break;
            case	TANTALUM	: return 1.70; break;
            case	TECHNETIUM	: return 1.47; break;
            case	TELLURIUM	: return 1.47; break;
            case	TERBIUM		: return 1.94; break;
            case	THALLIUM	: return 1.45; break;
            case	THORIUM		: return 2.06; break;
            case	THULIUM		: return 1.90; break;
            case	TIN		: return 1.39; break;
            case	TITANIUM	: return 1.60; break;
            case	TUNGSTEN	: return 1.62; break;
            case	URANIUM		: return 1.96; break;
            case	VANADIUM	: return 1.53; break;
            case	XENON		: return 1.50; break;
            case	YTTERBIUM	: return 1.87; break;
            case	YTTRIUM		: return 1.90; break;
            case	ZINC		: return 1.22; break;
            case	ZIRCONIUM	: return 1.75; break;
            default: break;
       }
        
        printf("Covalent parameters: Unsupported element proton number %d\n", N);
        exit(1);
        
	return 0;
}

double	giveAtomicWeight(unsigned char N) {
        switch(N) {
            case	ACTINIUM	: return 227	; break;
            case	ALUMINIUM	: return 26.982	; break;
            case	AMERICIUM	: return 243	; break;
            case	ANTIMONY	: return 121.76	; break;
            case	ARGON		: return 39.948	; break;
            case	ARSENIC		: return 74.922	; break;
            case	ASTATINE	: return 210	; break;
            case	BARIUM		: return 137.327; break;
            case	BERKELIUM	: return 247	; break;
            case	BERYLLIUM	: return 9.012	; break;
            case	BISMUTH		: return 208.98	; break;
            case	BOHRIUM		: return 264	; break;
            case	BORON		: return 10.811	; break;
            case	BROMINE		: return 79.904	; break;
            case	CADMIUM		: return 112.411; break;
            case	CAESIUM		: return 132.905; break;
            case	CALCIUM		: return 40.078	; break;
            case	CALIFORNIUM	: return 251	; break;
            case	CARBON		: return 12.011	; break;
            case	CERIUM		: return 140.116; break;
            case	CHLORINE	: return 35.453	; break;
            case	CHROMIUM	: return 51.996	; break;
            case	COBALT		: return 58.933	; break;
            case	COPPER		: return 63.546	; break;
            case	CURIUM		: return 247	; break;
            case	DARMSTADTIUM	: return 271	; break;
            case	DUBNIUM		: return 262	; break;
            case	DYSPROSIUM	: return 162.5	; break;
            case	EINSTEINIUM	: return 252	; break;
            case	ERBIUM		: return 167.26	; break;
            case	EUROPIUM	: return 151.964; break;
            case	FERMIUM		: return 257	; break;
            case	FLUORINE	: return 18.998	; break;
            case	FRANCIUM	: return 223	; break;
            case	GADOLINIUM	: return 157.25	; break;
            case	GALLIUM		: return 69.723	; break;
            case	GERMANIUM	: return 72.61	; break;
            case	GOLD		: return 196.967; break;
            case	HAFNIUM		: return 178.49	; break;
            case	HASSIUM		: return 269	; break;
            case	HELIUM		: return 4.003	; break;
            case	HOLMIUM		: return 164.93	; break;
            case	HYDROGEN	: return 1.008	; break;
            case	INDIUM		: return 114.818; break;
            case	IODINE		: return 126.904; break;
            case	IRIDIUM		: return 192.217; break;
            case	IRON		: return 55.845	; break;
            case	KRYPTON		: return 83.8	; break;
            case	LANTHANUM	: return 138.906; break;
            case	LAWRENCIUM	: return 262	; break;
            case	LEAD		: return 207.2	; break;
            case	LITHIUM		: return 6.941	; break;
            case	LUTETIUM	: return 174.967; break;
            case	MAGNESIUM	: return 24.305	; break;
            case	MANGANESE	: return 54.938	; break;
            case	MEITNERIUM	: return 268	; break;
            case	MENDELEVIUM	: return 258	; break;
            case	MERCURY		: return 200.59	; break;
            case	MOLYBDENUM	: return 95.94	; break;
            case	NEODYMIUM	: return 144.24	; break;
            case	NEON		: return 20.18	; break;
            case	NEPTUNIUM	: return 237	; break;
            case	NICKEL		: return 58.693	; break;
            case	NIOBIUM		: return 92.906	; break;
            case	NITROGEN	: return 14.007	; break;
            case	NOBELIUM	: return 259	; break;
            case	OSMIUM		: return 190.23	; break;
            case	OXYGEN		: return 15.999	; break;
            case	PALLADIUM	: return 106.42	; break;
            case	PHOSPHORUS	: return 30.974	; break;
            case	PLATINUM	: return 195.078; break;
            case	PLUTONIUM	: return 244	; break;
            case	POLONIUM	: return 210	; break;
            case	POTASSIUM	: return 39.098	; break;
            case	PRASEODYMIUM	: return 140.908; break;
            case	PROMETHIUM	: return 145	; break;
            case	PROTACTINIUM	: return 231.036; break;
            case	RADIUM		: return 226	; break;
            case	RADON		: return 222	; break;
            case	RHENIUM		: return 186.207; break;
            case	RHODIUM		: return 102.906; break;
            case	RUBIDIUM	: return 85.468	; break;
            case	RUTHENIUM	: return 101.07	; break;
            case	RUTHERFORDIUM	: return 261	; break;
            case	SAMARIUM	: return 150.36	; break;
            case	SCANDIUM	: return 44.956	; break;
            case	SEABORGIUM	: return 266	; break;
            case	SELENIUM	: return 78.96	; break;
            case	SILICON		: return 28.086	; break;
            case	SILVER		: return 107.868; break;
            case	SODIUM		: return 22.991	; break;
            case	STRONTIUM	: return 87.62	; break;
            case	SULPHUR		: return 32.066	; break;
            case	TANTALUM	: return 180.948; break;
            case	TECHNETIUM	: return 98	; break;
            case	TELLURIUM	: return 127.6	; break;
            case	TERBIUM		: return 158.925; break;
            case	THALLIUM	: return 204.383; break;
            case	THORIUM		: return 232.038; break;
            case	THULIUM		: return 168.934; break;
            case	TIN		: return 118.71	; break;
            case	TITANIUM	: return 47.867	; break;
            case	TUNGSTEN	: return 183.84	; break;
            case	URANIUM		: return 238.029; break;
            case	VANADIUM	: return 50.942	; break;
            case	XENON		: return 131.29	; break;
            case	YTTERBIUM	: return 173.04	; break;
            case	YTTRIUM		: return 88.906	; break;
            case	ZINC		: return 65.39	; break;
            case	ZIRCONIUM	: return 91.224	; break;
            default: break;
        }

        printf("Atomic Weight: Unsupported element proton number %d\n", N);
        exit(1);
        
	return 0;
}

unsigned char   isMetal(char N) {
             if(N >= LITHIUM    &&  N <= BERYLLIUM) return YES;
        else if(N >= SODIUM     &&  N <= ALUMINIUM) return YES;
        else if(N >= POTASSIUM  &&  N <= GERMANIUM) return YES;
        else if(N >= RUBIDIUM   &&  N <= ANTIMONY)  return YES;
        else if(N >= CAESIUM    &&  N <= POLONIUM)  return YES;
        else if(N >= FRANCIUM)                      return YES;
        else return NO;
}

double  my_dabs(double d) {
        if(d < 0) return -d;
        else return d;
}

void	toUpperCase(char * charPtr) { for (unsigned int i = 0; i < strlen(charPtr); i++) if(charPtr[i] >= 'a'  &&  charPtr[i] <= 'z') charPtr[i] = charPtr[i] - 32; }

unsigned char	giveRotBonds(char l){
	// for rigid residues: P, G, A, c : rotBonds = 0
	     if(l == 'V'  ||  l == 'm') return 1;
	else if(l == 'L'  ||  l == 'I'  ||  l == 'F'  ||  l == 'N'  ||  l == 'D'  ||  l == 'd'  ||  l == 'H'  ||
		l == 'h'  ||  l == 'B'  ||  l == 'b'  ||  l == 'S'  ||  l == 'T'  ||  l == 'W'  ||  l == 'C') return 2;
	else if(l == 'M'  ||  l == 'Y'  ||  l == 'Q'  ||  l == 'E'  ||  l == 'e') return 3;
	else if(l == 'K'  ||  l == 'R') return 4;
        // specials !!!
        else if(l == 'o') return 2;
        else if(l == 'u') return 3;
	else return 0;
}

void    residueName_L2C(char c, char * rname){
        switch(c) {
                case 'G' : strcpy(rname, "GLY"); break;
                case 'A' : strcpy(rname, "ALA"); break;
                case 'V' : strcpy(rname, "VAL"); break;
                case 'L' : strcpy(rname, "LEU"); break;
                case 'I' : strcpy(rname, "ILE"); break;
                case 'F' : strcpy(rname, "PHE"); break;
                case 'M' : strcpy(rname, "MET"); break;
                case 'C' : strcpy(rname, "CYS"); break;
                case 'c' : strcpy(rname, "CYX"); break;
                case 'm' : strcpy(rname, "CYM"); break;
                case 'Y' : strcpy(rname, "TYR"); break;
                case 'W' : strcpy(rname, "TRP"); break;
                case 'T' : strcpy(rname, "THR"); break;
                case 'S' : strcpy(rname, "SER"); break;
                case 'P' : strcpy(rname, "PRO"); break;
                case 'K' : strcpy(rname, "LYS"); break;
                case 'R' : strcpy(rname, "ARG"); break;
                case 'N' : strcpy(rname, "ASN"); break;
                case 'Q' : strcpy(rname, "GLN"); break;
                case 'D' : strcpy(rname, "ASP"); break;
                case 'd' : strcpy(rname, "ASH"); break;
                case 'E' : strcpy(rname, "GLU"); break;
                case 'e' : strcpy(rname, "GLH"); break;
                case 'H' : strcpy(rname, "HID"); break;
                case 'h' : strcpy(rname, "HIE"); break;
                case 'B' : strcpy(rname, "HIM"); break;
                case 'b' : strcpy(rname, "HIP"); break;
                // selenocysteine
                case 'o' : strcpy(rname, "SEC"); break;
                // selenomethionine
                case 'u' : strcpy(rname, "MSE"); break;
                
                default: printf("Unknown residue %c\n", c); exit(1); break;
        }
}

char    residueName_C2L(char * resName){
             if(strcmp(resName, "GLY") == 0) return 'G';
        else if(strcmp(resName, "ALA") == 0) return 'A';
        else if(strcmp(resName, "VAL") == 0) return 'V';
        else if(strcmp(resName, "LEU") == 0) return 'L';
        else if(strcmp(resName, "ILE") == 0) return 'I';
        else if(strcmp(resName, "PHE") == 0) return 'F';
        else if(strcmp(resName, "MET") == 0) return 'M';
        else if(strcmp(resName, "CYS") == 0) return 'C';
        else if(strcmp(resName, "CYX") == 0) return 'c';
        else if(strcmp(resName, "CYM") == 0) return 'm';
        else if(strcmp(resName, "TYR") == 0) return 'Y';
        else if(strcmp(resName, "TRP") == 0) return 'W';
        else if(strcmp(resName, "THR") == 0) return 'T';
        else if(strcmp(resName, "SER") == 0) return 'S';
        else if(strcmp(resName, "PRO") == 0) return 'P';
        else if(strcmp(resName, "LYS") == 0) return 'K';
        else if(strcmp(resName, "ARG") == 0) return 'R';
        else if(strcmp(resName, "ASN") == 0) return 'N';
        else if(strcmp(resName, "GLN") == 0) return 'Q';
        else if(strcmp(resName, "ASP") == 0) return 'D';
        else if(strcmp(resName, "ASH") == 0) return 'd';
        else if(strcmp(resName, "GLU") == 0) return 'E';
        else if(strcmp(resName, "GLH") == 0) return 'e';
        else if(strcmp(resName, "HID") == 0) return 'H';
        else if(strcmp(resName, "HIE") == 0) return 'h';
        else if(strcmp(resName, "HIP") == 0) return 'b';
        else if(strcmp(resName, "HIM") == 0) return 'B';
        else if(strcmp(resName, "HIS") == 0) return 'H';
        // selenocysteine
        else if(strcmp(resName, "SEC") == 0) return 'o';
        // selenomethionine
        else if(strcmp(resName, "MSE") == 0) return 'u';
        else if(strcmp(resName, "HOH") == 0) return '~';
        else if(strcmp(resName, "WAT") == 0) return '~';
        else if(strcmp(resName, "T3P") == 0) return '~';
        else {
                return '?';
                printf("Unknown residue name %s\n", resName);
                exit(1);
        }
}

#define	PI 			3.14159265358979
#define RAD_180                 3.14159265358979
#define	I80_DIV_PI             57.295779513082325

#define RAD_90                  PI / 2.0
#define RAD_60                  PI / 3.0
#define RAD_45                  PI / 4.0
#define RAD_30                  PI / 6.0
#define RAD_20          	PI / 9.0
#define RAD_15                  PI / 12.0
#define RAD_10                  PI / 18.0
#define RAD_9                   PI / 20.0
#define RAD_6                   PI / 30.0
#define RAD_5                   PI / 36.0
#define RAD_3                   PI / 60.0
#define RAD_2                   PI / 90.0
#define RAD_1                   PI / 180.0

#define RAD_360                 2.0 * PI
#define RAD_150                 5.0 * RAD_30
#define RAD_125                 25.0 * RAD_5
#define RAD_120                 2.0 * RAD_60
#define RAD_105                 7.0 * RAD_15
#define RAD_100                 10.0 * RAD_10
#define RAD_80          	4.0 * RAD_20
#define RAD_75                  5.0 * RAD_15
#define RAD_50                  5.0 * RAD_10
#define RAD_40                  2.0 * RAD_20
#define RAD_35                  7.0 * RAD_5
#define RAD_25                  5.0 * RAD_5
#define RAD_8                   4.0 * RAD_2
#define RAD_4                   2.0 * RAD_2

// for conf search
#define RAD_175                 17.5 * PI / 18.0
#define RAD_170                 17.0 * PI / 18.0
#define RAD_165                 16.5 * PI / 18.0
#define RAD_160                 16.0 * PI / 18.0
#define RAD_155                 15.5 * PI / 18.0
#define RAD_130                 13.0 * PI / 18.0
#define RAD_115                 11.5 * PI / 18.0
#define RAD_107                 107.0 * PI / 180.0
#define RAD_85                  8.5 * PI / 18.0
#define RAD_70                  7.0 * PI / 18.0
#define	RAD_65                  6.5 * PI / 18.0
#define RAD_55                  5.5 * PI / 18.0
#define	RAD_7                   7.0 * PI / 180.0

#endif	/* DEFINITIONS_H */

