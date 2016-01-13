/*
 
 Ablator++
 
 ENERGY-BASED WALL RESPONSE MODEL
 1-D TRANSIENT HEAT CONDUCTION
 1-D FINITE DIFFERENCE WAVE PROPAGATION
 TEMPERATURE-DEPENDENT MATERIAL PROPERTIES
 HEAT GENERATION FROM X-RAY DEPOSITION
 IN A DOUBLE EXPONENTIAL PULSE
 OR IN SQUARE OR GAUSSIAN PULSES
 USES 45 BIN APPROXIMATION TO SPECTRUM
 SUPPLY COLD OPACITIES IN INVERSE METERS
 EXPLICIT SCHEME
 SIMPLE RATIO ZONING
 UNITS: SI, KEV FOR BBT STUFF
 
 Translation to C++ by Michal Vasinek <michal.vasinek@vsb.cz>
 and Martin Stachon <martin.stachon@vsb.cz>
 
 */

#include <cmath>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <algorithm>

#include "fortranio.h"

using namespace std;

//  IMPLICIT DOUBLE PRECISION(A-H,O-Z) // TODO: all variables beggining with letters a-h, o-z are of type double !

#define N 100
#define NP1 101
#define NBIN 45
#define NSAVE 5
#define PI 3.14159265358979323846

#define min5(a,b,c,d,e) min((a), min((b), min((c), (d))))

double H[N]; // TOTAL INTERNAL ENERGY IN ZONE
double HG[N]; // SPECIFIC INTERNAL ENERGY FOR A ZONE
double T[NP1+1]; //TODO zero indexed in fortran! // TEMPERATURE IN A ZONE
double TNEW[NP1+1]; //TODO zero indexed in fortran! // NEW TEMPERATURE (TO BE APPLIED NEXT TIME STEP)
double ZMASS[NP1]; // MASS IN A ZONE (CONSTANT THROUGH PROBLEM)
double OPAC[NBIN]; // OPACITY IN INVERSE METERS AT A GIVEN PHOTON ENERGY
double BIN[NBIN]; // CENTRAL ENERGY (IN KEV) OF A PHOTON BIN
double BINW[NBIN]; // PHOTON BIN WIDTH
double EFRAC[N][NBIN]; // FRACTION OF INCIDENT BIN ENERGY EIN(J) ABSORBED IN A ZONE
double EZONE[N]; // FRACTION OF EIN DEPOSITED IN A ZONE
//double ELEAK[NBIN]; // RUNNING SUM OF LEAK ENERGY IN A BIN
//double EIN[NBIN]; // TOTAL ENERGY ADDED THIS STEP IN A BIN
//double ELEAKOLD[NBIN]; // SUM FROM THE LAST TIME STEP
double FBBLEAK[NBIN]; // FRACTION OF LEAK ENERGY IN A BIN
double XLOC[NP1]; // CENTER OF ORIGINAL ZONE (MICRONS)
//double ACOND[13]; // COEFFICIENTS FOR THERMAL CONDUCTIVITY VS TEMP
//double AH2T[10]; // COEFFICIENTS FOR TEMPERATURE-ETHALPY RELATIONS
//double EOSCOEF[10]; // COEFFICIENTS FOR EQUATIONS OF STATE
//double SIGMACOEF[7]; // COEFFICIENTS FOR SURFACE TENSION (AKA SIGMA)
//double YIELDCOEF[3];
double A[N+1]; //TODO zero indexed in fortran! // ACCELERATION OF A NODE
double U[NP1+1]; //TODO zero indexed in fortran! // VELOCITY OF A NODE
double X[NP1+1]; //TODO zero indexed in fortran! // LOCATION OF A NODE
double X0[NP1+1]; //TODO zero indexed in fortran! // ORIGINAL LOCATION OF A NODE
double XC[NP1]; // CURRENT CENTER LOCATION OF A ZONE
double XNEW[N+1]; //TODO zero indexed in fortran! // NEW NODE LOCATION (TO BE APPLIED NEXT TIME STEP)
double SX[NP1+1]; //TODO zero indexed in fortran! // NORMAL STRESS
double SXNEW[NP1+1]; //TODO zero indexed in fortran! // NEW NORMAL STRESS (TO BE APPLIED NEXT TIME STEP)
double P[N+1]; //TODO zero indexed in fortran! // PRESSURE
double SXD[N]; // DEVIATORIC STRESS
double SXDNEW[N]; // NEW DEVIATORIC STRESS (TO BE APPLIED NEXT TIME STEP)
double SZD[N]; // DEVIATORIC STRESS (FOR CYLINDRICAL GEOMETRIES)
double SZDNEW[N]; // NEW DEVIATORIC STRESS (TO BE APPLIED NEXT TIME STEP)
double PHI[NP1+1]; //TODO zero indexed in fortran! // DIFFERENCE BETWEEN LONGITUDINAL AND TRANSVERSE STRESS

double RHO[NP1+1]; //TODO zero indexed in fortran! // DENSITY
double RHONEW[NP1+1]; //TODO zero indexed in fortran! // NEW DENSITY (TO BE APPLIED NEXT TIME STEP)
double Q[NP1+1]; //TODO zero indexed in fortran! // ARTIFICIAL VISCOSITY
double C[N]; // LOCAL SOUND SPEED
double XQUAL[NP1]; // QUALITY (VAPOR MASS FRACTION) IN A ZONE
double VOLF[NP1]; // VAPOR VOLUME FRACTION IN A ZONE
double RATEJ[N]; // LOG OF BUBBLE NUCLEATION RATE
double TNUC[N]; // LOG OF INDUCTION TIME FOR NUCLEATION RATE DEVELOPMENT (STEADY?)
double BUBBLE[N]; // LOG OF CUMULATIVE NUMBER OF BUBBLES PER UNIT VOLUME
double QP1[N]; // QUANTITY FOR MOSS & WHITE ARTIFICIAL VISCOSITY
double QP2[N]; // QUANTITY FOR MOSS & WHITE ARTIFICIAL VISCOSITY
//double SURFCSTUFF[10]; // QUANTITIES FOR SURFACE ZONE CONDUCTION ITERATIONS
int PLASTIC[N]; // FLAG 1 FOR ZONE WHERE PLASTICITY OCCURED, OTHERWISE 0
double RATIO[N];
double TPP[N];
double TPPP[N];
double TP[N];
double HINP[N];
double HCONDP[N];
double HDEVP[N];
double HSTRESSP[N];
double HSTRESSQP[N];
double F1P[N];
double F2P[N];
double DELRHOP[N];

int GEOM, QFLAG, USEYIELDFLAG;
char TAB;
char JUNK;
string MATNAME;
double MOMSUM, MYERF;

//variable to use in DENGEOM, it is not defined in former fortran implicitly,
//but rather by [O-Z] notation
double RDEN;
//from subroutine AXIS
double DTMIN,DTMAX;
//from subroutine XSOURCE
double TSTART;
//from subroutine XPULSE
double TNS;

//tohle jsou hlavicky funkci
//EXTERNAL TSURF, PSURF, ZBRENT, CONDSURF, QTWO, ERF
//prevedeme je na ukazatele, v nekterych mistech napr EOSURFT subroutine predava
//ukazatele na funkce jako argument
//double (*TSURF)(double);
//double (*PSURF)(double);

// common blocks converted to global struct
// TODO rework this
struct PULSE {
  double TSQ;
  double ESQ;
  double EGAUSS;
  double SC;
  double FWHM;
}s_pulse;
struct MATDAT {
  double TMELT;
  double HMELTI;
  double HMELTT;
  double ACOND[16]; // COEFFICIENTS FOR THERMAL CONDUCTIVITY VS TEMP
  double AH2T[10]; // COEFFICIENTS FOR TEMPERATURE-ETHALPY RELATIONS
}s_matdat;
struct EOSDAT {
  double RHO0;
  double EOSCOEF[10]; // COEFFICIENTS FOR EQUATIONS OF STATE
  double ABOIL;
  double BBOIL;
  double PRNU;
  double YIELDCOEF[3];
}s_eosdat;
struct EINTEGRAL {
  double ELEAKOLD[NBIN]; // SUM FROM THE LAST TIME STEP
  double ELEAK[NBIN]; // RUNNING SUM OF LEAK ENERGY IN A BIN
  double EIN[NBIN]; // TOTAL ENERGY ADDED THIS STEP IN A BIN
}s_eintegral;
struct OPACBINS {
  double BIN[NBIN]; // CENTRAL ENERGY (IN KEV) OF A PHOTON BIN
  double BINW[NBIN]; // PHOTON BIN WIDTH
}s_opacbins;
struct SURFT {
  double HGSURF;
  double QUALSURFT;
}s_surft;
struct SURFP {
  double RHOSURF;
  double QUALSURFP;
  double HGVSURF;
  double HGLSURF;
  double TEMPSURF;
}s_surfp;
struct SURFC {
  double SURFCSTUFF[10]; // QUANTITIES FOR SURFACE ZONE CONDUCTION ITERATIONS
}s_surfc;
struct SURFTEN {
  double SIGMACOEF[7]; // COEFFICIENTS FOR SURFACE TENSION (AKA SIGMA)
}s_surften;

//to use in function QTWO and subroutine EOSPALL
struct QTWODAT {
  double HGTWO;
  double RHOTWO;
  double RHOV;
  double RHOL;
  double TTWO;
  double HGL;
}s_qtwodat;

//functions prototypes
void READMAT();
void ZONESET(double& rho0,double x0[], double xc[], double xloc[],double zmass[],int& geom, ifstream& file);
void XSOURCE(int& isource, int& laserflag, double fbbleak[], double& tstart,ifstream& file);
void H2T(double hjkg, double& t);
void CONDVST(double& t,double& volf, double& cond);
void EOSF(double& rho, double& f1, double& f2);
void SOUND(double& HG, double& RHO, double& P, double& F2, double& C, double& TEMP);
void XPULSE(double& tm, int& isource, double& area, double fbbleak[]);
void BETAGEOM(int& geom, double& r1, double& r2, double& betac);
void USURFIX(double& T, double& P0, double& P, double& UNEW);
void AXIS(double& x0, double& u0, double& a0, double& dtold, double& dtnew);
void PLEFTCALC(double& PSAT,double& PC,double& DXC1,double& DX1,double& RHOV,double& QUAL,double& T,double& PLEFT);
void DENGEOM(int& geom, double& r1, double& r2, double& rden);
void FSURF(double& RHO, double& RHOL, double& QUAL, double& QMIN, double& F1, double& F2);
void FVAPOR(double& RHO,double& F1, double& F2);
void F2PHASE(double& RHO,double& HG,double& QUAL, double& F1, double& F2);
double CONDSURF(double HEAST);
double ZBRENT(double (*FUNC)(double),double& X1,double& X2,double& TOL);
void BRACKET(double (*FUNC)(double),double& XNOM,double& X1,double& X2);
void EOSURF(double& RHOL,double& RHOV,double& RHOP1,double& QMIN,double& PM1,double& P,double& ALPHA,double& CSOUND);
void EOSPALL(double& HG,double& RHO,double& T,double& P,double& QUAL,double& VOLF,double& CSOUND);
void NUCLEATE(double& T,double& RHO,double& P,double& J,double& TSTAR);
void SURFVAP(int& GEOM, double& PAMB, double& TSAT, double& X, double& DT, double& DVAPMASS, double& DHVAP);
void HEATCAP(double& volf, double& t, double& hg, double& cv);
double QTWO(double QUAL);
void NEWTP(double& HG, double& P, double& RHO);
double TSURF_DEF(double d);
double PSURF_DEF(double d);
double SIGMACALC(double T);

void Binary(double x)
{
	double tot = 0;
	double px = 0;
	while (x / pow(2, px) > 1) {
		px += 1;		
		tot = px;
	}
	double rem = tot;	
	px = 1;
	while (rem <= 64) {
		cout << (int)(x / pow(2, -px));
		x -= (x/pow(2, -px))*pow(2, -px);
		px += 1;
		rem += 1;
	}
	cout << endl;
}

// DEFINE A BLACKBODY FUNCTION
double BBE(int JBB, double TBB) {
  return s_opacbins.BINW[JBB] * (pow((s_opacbins.BIN[JBB]/TBB),3)) / (exp(s_opacbins.BIN[JBB]/TBB)-1.0);
}

int main(int argc, char **argv) {
  TAB = '\t';
  ifstream file9; //initdata //TODO rename
  ifstream file8; //opacdata //TODO rename
  ofstream file10; // ldepfrac, also Pfront //TODO rename
  ofstream file11; // Pprofile
  ofstream file13; // SXProfile
  ofstream file14; // SXDprofile
  ofstream file27; // Qprofile
  ofstream file26; // Aprofile
  ofstream file15; // Uprofile
  ofstream file16; // Xprofile
  ofstream file17; // RHOprofile
  ofstream file18; // QUALprofile
  ofstream file19; // Tprofile
  ofstream file20; // HGprofile
  ofstream file21; // spall
  ofstream file22; // VOLFprofile
  ofstream file23; // summary
  ofstream file24; // TNUCprofile
  ofstream file25; // Jprofile
  ofstream file35; // BUBprofile
  ofstream file28; // PLASTICprofile
  ofstream file29; // RATIOprofile
  ofstream file30; // TPPprofile
  ofstream file31; // TPPPprofile
  ofstream file32; // EINprofile
  ofstream file33; // TPprofile
  ofstream file34; // ESUMDIFprofile
  ofstream file36; // ESUMprofile
  ofstream file37; // HINprofile
  ofstream file38; // HCONDprofile
  ofstream file39; // HDEVprofile
  ofstream file40; // HSTRESSprofile
  ofstream file41; // F1profile
  ofstream file42; // F2profile
  ofstream file43; // DELRHOprofile
  ofstream file44; // EKINprofile
  ofstream file45; // ESUMVAPprofile
  ofstream file46; // ESUMSLprofile
  ofstream file47; // ESUMDIFVAPprofile
  ofstream file48; // ESUMDIFSLprofile
  ofstream file49; // MOMprofile
  char line_buffer[1024]; // for sprintf lines, I prefer to use sprintf over C++ formatting
  
  string COMMENT;
  int LASERFLAG;
  int ISOURCE;
  double RADIUS, AREA;
  int ISPALL, ISURF;
  double PAMB;
  double TMELT;
  double DTVAP;
  double DVAPMASS;
  double DXQ;
  double XQ;
  double RHOV;
  double DQI;
  double VAPRATIO;
  double DVMNEXT;
  double HSUM;
  double EKIN;
  double HSTRESS;
  double DHVAP;
  double DHVAPN;
  double UBAR;
  double ESUM;
  int K1, K2;
  double DTI;
  double DXHEAT;
  double TEMPBIG;
  //unreferenced
  //double DTTTEMP;
  double NCYMOD;
  double DTTEMP;
  double CV;
  double RAWJ;
  double RAWTNUC;
  double RAWBUB;
  double DTPRESSEIN;
  double TPULSE;
  double DTPRESS;
  double DX;
  double RELDTLIMIT;
  double TMV;
  double HINSUM;
  double HCONDSUM;
  double HDEVSUM;
  double HSTRESSSUM;
  double DTSCREEN;
  double X0M1;
  double VAPDEPTH;
    
  // READ RADIATION SOURCE DATA FROM A FILE
  file9.open("initdata", ifstream::in);
  //file40.open("DH", ifstream::out);
  
  //  WRITE(*,*) 'BB X-ray or laser source? (1=Laser)'
  //  READ(*,*) LASERFLAG
  COMMENT = FortranIO::read_string(file9);
  LASERFLAG = FortranIO::read_int(file9);
  
  if (LASERFLAG != 1) {
    file8.open("opacdata", ifstream::in);
    MATNAME = FortranIO::read_string(file8);
    cout << MATNAME << endl;
    for(int j=0;j<NBIN;j++){
      s_opacbins.BIN[j] = FortranIO::read_double(file8);
      s_opacbins.BINW[j] = FortranIO::read_double(file8);
      OPAC[j] = FortranIO::read_double(file8);
	  OPAC[j] = OPAC[j] * 1e-6;
    }
    file8.close();
    //  SKIP READING ABS.COEFF. FROM THE FILE
    COMMENT = FortranIO::read_string(file9);
    COMMENT = FortranIO::read_string(file9);
    
  } else {
    //    WRITE(*,*) 'Enter absorbtion coefficient in reciprocal microns: '
    //    READ(*,*) OPACLASER

    COMMENT = FortranIO::read_string(file9);
    double OPACLASER = FortranIO::read_double(file9);
    
    for (int j=0; j<NBIN; j++) {
      s_opacbins.BIN[j] = 0;
      s_opacbins.BINW[j] = 0;
      OPAC[j] = OPACLASER * 1e-06;
    }
  }
  
  //  READ IN MATERIAL PROPERTIES
  READMAT();
  
  double& GAMMA = s_eosdat.EOSCOEF[8];
  double& RBAR = s_eosdat.EOSCOEF[9];
  
  //  SHEAR MODULUS COEFFICIENT, BASED ON POISSON'S RATIO
  double GNU = (1.0-2.0*s_eosdat.PRNU) / (2.0*(1.0-s_eosdat.PRNU));
  //  SOUND SPEED COEFFICIENT, BASED ON POISSONS' RATIO
  double SNU = sqrt(3.0* (1.0-s_eosdat.PRNU)/(1.0+s_eosdat.PRNU));
  
  //  CHOOSE ARTIFICIAL VISCOSITY METHOD
  // WRITE(*,*) 'Artificial viscosity method'
  // WRITE(*,*) ' 1 = VNR'
  // WRITE(*,*) ' 2 = M&W'
  // READ(*,*) QFLAG
  QFLAG = 1;
  
  //  INITIALIZE ARTIFICIAL VISCOSITY COEFFICIENTS
  //  B1, B2 ARE FOR SOLID, B3 FOR SPALLED MATERIAL (NOT USED)
  double B1 = 1.000;
  double B2 = 0.300;
  //  B3 = 0.060
  
  //  INITIALIZE TIME STEP CONTROL COEFFICIENTS
  double CT1 = 0.9;
  double CT2 = 1.1;
  //  CFL MUST BE < 0.5
  double CFL = 0.4;
  
  //  SET SPALL STRESS (FOR REMOVAL OF MELTED MATERIAL) (in Pa)
  double SPALL = -1.0e+10;
  
  //  SET UP ZONES
  //  GEOM IS FLAG FOR PLANAR (1), CYLINDRICAL (2) OR SPHERICAL (3)
  
  /*   POZNAMKA PRO DRUIDA, tady to bude treba rozsirit a predat jako
   posledni parametr soubor cislo devet, v zonesetu
   */
  ZONESET(s_eosdat.RHO0,X0,XC,XLOC,ZMASS,GEOM,file9);
  if (GEOM == 1)
    cout << "PLANAR" << endl;
  if (GEOM == 2)
    cout << "CYLINDRICAL" << endl;
  if (GEOM == 3)
    cout << "SPHERICAL" << endl;
  if (GEOM < 1 || GEOM > 3) {
    cout << "Invalid GEOM = " << GEOM;
    return 1;
  }
  
  //  CHOOSE SOURCE TERM
  do {
    //  WRITE(*,*) 'Choose pulse type'
    //  WRITE(*,*) '1 = square pulse'
    //  WRITE(*,*) '2 = Gaussian pulse'
    //  READ(*,*) ISOURCE
    COMMENT = FortranIO::read_string(file9);
    ISOURCE = FortranIO::read_int(file9);
  } while (ISOURCE > 2 || ISOURCE < 1);
  
  //  SET UP FOR ENERGY DEPOSITION BASED ON ISOURCE CHOICE
  XSOURCE(ISOURCE,LASERFLAG,FBBLEAK,TSTART,file9);
  
  //  WRITE(*,*) 'Enter radius (microns)'
  //  READ(*,*) RADIUS
  COMMENT = FortranIO::read_string(file9);
  RADIUS = FortranIO::read_double(file9);
  AREA = PI * pow(RADIUS * 1e-06, 2);
  
  //XUV
  /*
  COMMENT = FortranIO::read_string(file9);
  USEYIELDFLAG = FortranIO::read_int(file9);
  cout << "yield: " << USEYIELDFLAG << endl;
  */
  
  int LHYDRO = 1;
  
  double HG0 = 0.0;
  double& TINF=s_matdat.AH2T[0];
  
  //  DETEREMINE INITIAL TIME STEP
  //  ROOM TEMPERATURE MATERIAL PROPERTIES, TO GET THINGS STARTED
  //  FIRST GET HEAT CAPACITY AT ROOM TEMPERATURE
  double T0;
  H2T(HG0,T0);
  double H1 = HG0 + 0.1;
  double T1;
  H2T(H1,T1);
  double CVS = (H1 - HG0) / (T1 - T0);
  //  NOW GET HEAT CAPACITY OF LIQUID
  double CVL = 1.0 / s_matdat.AH2T[7];
  //  CHOOSE MINIMUM FOR USE IN TIME STEP CALCULATION
  double CV0 = min(CVS,CVL);
  
  //  J/kg.K
  double COND;
  double TEMP_VOLF=0.0;
  CONDVST(TINF,TEMP_VOLF,COND);
  //  W/m.K
  double ALPHA = COND / (s_eosdat.RHO0 * CV0);
  //  m2/sec
  double DTHEAT = CFL* pow(X0[1]-X0[0], 2) / ALPHA;
  
  //  SET PROBLEM TIME (sec) AND FLAGS FOR OCCASSIONAL DATA DUMPS (ns)
  //  WRITE(*,*) 'Enter Tstop (ns)'
  //  READ (*,*) TSTOP
  COMMENT = FortranIO::read_string(file9);
  double TSTOP = FortranIO::read_double(file9);
  file9.close();
  double TM = TSTART * 1.0e-09;
  TSTOP = TSTOP * 1.0e-09 + 2.0*DTHEAT;
  double TFRONT = TSTART - 1.0e-13;
  double TSCREEN = TSTART - 1.0e-13;
  
  double XMELTMAX = 0.0;
  double XTMELTMAX = 0.0;
  double ITMELT = 0.0;
  double TFMAX = TINF;
  double HGFMAX = HG0;
  double DHGMXMAX = 0.0;
  
  //  INITIALIZE ENERGIES AND TEMPERATURES
  //  ZONE ENERGY STORED IN H(), SPECIFIC ENERGY STORED IN HG()
  for (int i=0; i<N; i++) {
	//XUV
	/*
    HINP[i] = 0.0;
    HCONDP[i] = 0.0;
    HDEVP[i] = 0.0;
    HSTRESSP[i] = 0.0;
    HSTRESSQP[i] = 0.0;
	*/
    H[i] = HG0 * ZMASS[i];
    HG[i] = HG0;
    T[i+1] = TINF;
  }
  T[0] = TINF;
  T[NP1] = TINF;
  
  //  INITIALIZE HYDRO-MOTION VARIABLES
  double F1,F2;
  EOSF(s_eosdat.RHO0,F1,F2);
  double P0 = F1 + F2 * HG0;
  double C0, C0_TEMP;
  SOUND(HG0,s_eosdat.RHO0,P0,F2,C0,C0_TEMP);
  //  CONVERT SOUND TO ELASTIC FROM HYDRODYNAMIC MATERIAL
  C0 = C0 * SNU;
  for (int i=1; i<=N; i++) {
    A[i] = 0.00e+00;
    U[i] = 0.00e+00;
    X[i] = X0[i];
    RHO[i] = s_eosdat.RHO0;
    Q[i] = 0.00e+00;
    P[i] = P0;
    SX[i] = 0.00e+00;
    SXD[i-1] = 0.00e+00;
    PHI[i] = 0.00e+00;
    C[i-1] = C0;
    XQUAL[i-1] = 0.00e+00;
    VOLF[i-1] = 0.00e+00;
    QP1[i-1] = 0.0e+00;
    QP2[i-1] = 0.0e+00;
    BUBBLE[i-1] = 0.0e+00;
  }
  
  //  ASSUME FREE BOUNDARY AT FRONT
  U[0] = 0.00e+00;
  X[0] = X0[0];
  RHO[0] = 0.00e+00;
  RHONEW[0] = 0.00e+00;
  SX[0] = 0.00e+00;
  P[0] = 0.00e+00;
  PHI[0] = 0.00e+00;
  Q[0] = 0.00e+00;
  //  ASSUME NON-REFLECTING BOUNDARY AT BACK
  U[NP1] = 0.00e+00;
  double UNOLD = 0.00e+00;
  X[NP1] = X0[NP1];
  RHO[NP1] = 0.00e+00;
  SX[NP1] = 0.00e+00;
  PHI[NP1] = 0.00e+00;
  Q[NP1] = 0.00e+00;
  XQUAL[NP1-1] = 0.00e+00;
  VOLF[NP1-1] = 0.00e+00;
  
  //  INITIALIZE LIQUID DENSITY (FOR SURFACE VAPORIZATION)
  double RHOLSURF = s_eosdat.RHO0;
  
  //  SET INITIAL SPALL ZONE AND FRONT SURFACE ZONE
  ISPALL = 0;
  ISURF = 1;
  
  //  SET MIN QUALITY FOR SURFACE ZONE TO BE CONSIDERED W/O VAPOR
  double QMIN = 1.0e-06;
  
  //  SET INITIAL TIME STEP SIZE
  double DTHYDRO = (X0[1]-X0[0])/C0;
  double DT = min(DTHEAT,DTHYDRO);
  if (LHYDRO != 1) DT = DTHEAT;
  double DTOLD = DT;
  
  //  SET UP SOURCE ENERGY MULTIPLIERS
  //  WILL GIVE FRACTION OF ENERGY IN EACH ZONE, IN EACH BIN
  //  OF THE TOTAL INCIDENT FROM WALL AND LEAK X RAYS
  //  SCALE LEAK ENERGY WITH LAMBERTIAN DISTRIBUTION
  for (int j=0; j<NBIN; j++) {
    double FRACIN = 1.0;
    for (int i=0; i<N; i++) {
	  //XUV
      //double FRACOUT = exp(-X[i+1] * OPAC[j]);

	  double FRACOUT = exp(-X[i+1] / OPAC[j]);
      EFRAC[i][j] = FRACIN - FRACOUT;
      FRACIN = FRACOUT;
    }
    
    //  INITIALIZE ENERGIES/BIN FOR FIRST TIME STEP
    s_eintegral.ELEAKOLD[j] = 0.0;
    //    WRITE(*,*) EFRAC(1,J)
  }
   
  //XUV
  /*
  for (int j=0; j<NBIN; j++) {    
	double SUMEFRAC = 0.0;
    for (int i=0; i<N; i++) {
      SUMEFRAC = SUMEFRAC + EFRAC[i][j];
    }
    //    WRITE(*,*) J,SUMEFRAC
  }
  */
  
  //  SAVE ENERGY FRACTION IN EACH ZONE
  for (int i=0; i<N; i++) {
    EZONE[i] = 0.0e+00;
    for (int j=0; j<NBIN; j++) {
      EZONE[i] = EZONE[i] + EFRAC[i][j];
    }
  }
  
  //  SUMALEAK = 0.
  //  WRITE(*,*) '3.FBBLEAK(J):'
  //  DO J=1,NBIN
  //    SUMALEAK = SUMALEAK + FBBLEAK(J)
  //    WRITE(*,*) FBBLEAK(J)
  //  END DO
  //  WRITE(*,*) '3.SUMALEAK=',SUMALEAK
  
  //  FOR LASER SOURCE, RECORD ENERGY DEPOSITION PROFILE
  if (LASERFLAG == 1) {
    file10.open("ldepfrac", iostream::out);
    for (int i=1; i<=N; i++) {
      file10 << i << "\t" << X0[i]*1.0e+06 << "\t" << EFRAC[i][0] << endl;
    }
    file10.close();
  }
  
  file10.open("Pfront");
  file11.open("Pprofile");
  file13.open("SXProfile");
  //XUV
  //file14.open("SXDprofile");
  file14.open("Cprofile");
  file27.open("Qprofile");
  file26.open("Aprofile");
  file15.open("Uprofile");
  file16.open("Xprofile");
  file17.open("RHOprofile");
  file18.open("QUALprofile");
  file19.open("Tprofile");
  file20.open("HGprofile");
  file21.open("spall");
  file22.open("VOLFprofile");
  file24.open("TNUCprofile");
  file25.open("Jprofile");
  file35.open("BUBprofile");
  file40.open("RTinstab");
  file41.open("F1profile");
  file42.open("F2profile");
  file43.open("DELRHOprofile");
  file44.open("EKINprofile");
  file45.open("ESUMVAPprofile");
  file46.open("ESUMSLprofile");
  file47.open("ESUMDIFVAPprofile");
  file48.open("ESUMDIFSLprofile");
  file49.open("MOMprofile");

  ofstream ifs2("debug-X-isurf");
  ofstream ifs3("debug-tm-isurf");
  ofstream ifs4("debug-volf-isurf");
  ofstream ifs5("debug-dt-isurf");
  ofstream ifs6("debug-dtold-isurf");
  ofstream ifs7("debug-zmass-isurf");
  
  //int ISURFEND = 23;
  int EKINSIGN = 1;
  double EINMAX = 0.0;
  //double EINMAXVAP = 0.0;
  //double EINMAXSL = 0.0;
  DTMIN = 100.0;
  DTMAX = 0.0;
  double HEASTMAX = 0.0;
  DHGMXMAX = 0.0;
  int IFLAGAXIS;
  double DTEMPMAX;
  double DTEMP;
  int ITEMPMAX;
  double TMHIGH, TMLOW, XTMELT, XMELT;
  if (GEOM > 1 && X0[0] == 0.0) {
    IFLAGAXIS = 1;
  } else {
    IFLAGAXIS = 0;
  }
  
  //  INITIALIZE MAX BUBBLE NUCLEATION DEPTHS AND ZONE NUMBERS
  int I35 = 1;
  int I30 = 1;
  int I25 = 1;
  int I20 = 1;
  double X35 = 0.0;
  double X30 = 0.0;
  double X25 = 0.0;
  double X20 = 0.0;
  int IB27 = 1;
  int IB26 = 1;
  int IB25 = 1;
  int IB24 = 1;
  double XB27 = 0.0;
  double XB26 = 0.0;
  double XB25 = 0.0;
  double XB24 = 0.0;
  //XUV
  //  INIT PLASTICITY FLAG
  //int PLASTICFLAG = 0;
  
  //  SET UP TO PRINT INFO TO SCREEN EVERY NCYCLEPRINT CYCLES
  int NCYCLEPRINT = 10000;
  int NCYCLE = 1;
  double TPRESS = TM * 1.0e+09;

  //XUV
  //double TPRESSEIN = TPRESS;
  
  //  *******************************************************
  //  START MARCHING IN TM
  // MV - this variable is potentially used before assigned
  double RHOVSURF;

  while (true) {
    TM = TM + DT;
    if (TM > TSTOP) break;
    
    //  SET UP MAGNITUDE OF ENERGY PULSES
    //  RETURNS ENERGY IN EACH BIN FOR THIS TIME STEP, EIN(J) IN COMMON
    XPULSE(TM,ISOURCE,AREA,FBBLEAK);
    
    //  ADD THE TOTAL ENERGY DEPOSITED THIS STEP TO RUNNING TOTAL
    double EINDT=0.0;
    for (int J=1; J<NBIN; J++) {
      EINMAX = EINMAX + s_eintegral.EIN[J];
      EINDT = EINDT + s_eintegral.EIN[J];
    }
    
    //  RESET MAX SPECIFIC ENERGY CHANGE
    double DHGMAX = 0.0;
    //  RESET SPALL FLAG
    int ISPALLFLAG = 0;
    
    //  SET UP THERMAL CONDUCTIVITIES FOR THE FIRST ZONE
    double VFLAG = VOLF[0];
    if (ISURF == 1) VFLAG = 0.0;
    double CONDE;
    CONDVST(T[1],VFLAG,CONDE);
    
    double BETAC;
    BETAGEOM(GEOM,X[0],X[1],BETAC);
    double AE = CONDE * BETAC / (X[1] - X[0]);
    //  INSULATED FRONT SURFACE (CONDUCTION BETWEEN ZONES 0 AND 1)
    double HEAST = 0.0;
	//XUV
	/*
    //  SET PLASTICITY FLAG FOR EACH ZONE
    for (int i=0; i<N; i++) {
      PLASTIC[i] = 0;
    }
	*/
    
    if (LHYDRO == 1) {
      //  SET UP FOR FINDING NEW HYDRO-BASED MAXIMUM DT
      DTHYDRO = 100.0;
      if (ISURF == 1 && XQUAL[ISURF-1] > QMIN) {
        //  DETERMINE LEFT NODE MOTION OF SURF ZONE W/O F=MA
        double TRIGHT = T[1];
        double PRIGHT = P[1];
        double PLEFT = P[0];
        double UNEW;
        USURFIX(TRIGHT,PRIGHT,PLEFT,UNEW);
        A[0] = (UNEW - U[0]) / DT;
        U[0] = UNEW;
        XNEW[0] = X[0] + DT * U[0];
      } else {
        if (GEOM == 1) {
          //  DETERMINE ACCELERATION, VELOCITY, LOCATION FOR LEFT EDGE
          A[0] = 2.0 * (0.0 - (SX[1]+Q[1])) /
          (RHO[1] * (X[1] - X[0]) + 0.0);
          U[0] = U[0] + A[0] * (DT + DTOLD)/2.0;
          XNEW[0] = X[0] + DT * U[0];
        } else {
          //  CHECK IF INNER RADIUS HAS ALREADY HIT AXIS
          if (IFLAGAXIS == 0) {
            A[0] = 2.0 * (0.0 - (SX[1]+Q[1])) /
            (RHO[1]*(X[1]+X[0]) + 0.0) +
            2.0 * (GEOM-1) * (PHI[1] + PHI[0]) /
            (RHO[1]*(X[1] + X[0]) + 0.0);
            double UNEW = U[0] + A[0] * (DT + DTOLD)/2.0;
            double XNEW0 = X[0] + DT * UNEW;
            //  CHECK IF SURFACE HAS PASSED THROUGH AXIS
            if (XNEW0 <= 0.0) {
              //  KEEPING SAME ACCELERATION, FIGURE DT TO WHEN STUFF HITS AT CENTER
              double DTNEW;
              AXIS(X[0],U[0],A[0],DTOLD,DTNEW);			  
              DT = DTNEW;
              XNEW[0] = 0.0;
              U[0] = 0.0;
              IFLAGAXIS = 1;
            } else {
              U[0] = UNEW;
              XNEW[0] = XNEW0;
            }
          }
        }
      }
    }
    
    //  %%%%%%%%%%%%%%%%%%%%%% START OF LOOP THROUGH ZONES %%%%%%%%%%%%%%%%%%%%%%
    //  FRONT BC IS THERMALLY INSULATED
    //  REAR BC IS AT FIXED TEMPERATURE
    //  FRONT SURFACE IS FREE TO MOVE
    //  REAR SURFACE IS NON-REFLECTING
    
    for (int I=1; I<=N; I++) {		
      double HSTRESS;
      double HSTRESSQ;
      double HDEV;
      double DELRHO;
      double RHODOT;
      double HT;
      double QP;
      if (LHYDRO == 1) {
        //  COMPUTE ACCELERATION, VELOCITY, LOCATION		  
        if (I == (ISURF-1) && XQUAL[ISURF-1] > QMIN &&
            VOLF[I-1] > 0.99) {
          //  DETERMINE RIGHT NODE MOTION OF ISURF ZONE W/ ADJUSTED PRESSURE
          double PRIGHT = P[ISURF];
          double DXSURF = (X[ISURF]-X[I]) * (1.0-VOLF[ISURF-1]);
          double XINTERNAL = X[ISURF] - DXSURF;
          double DXC1 = XINTERNAL - XC[I-1];
          double DX1 = XINTERNAL - X[I];
          
          double PLEFT;
          PLEFTCALC(PRIGHT,P[I],DXC1,DX1,RHOVSURF,XQUAL[ISURF-1],T[ISURF],PLEFT);
          A[I] = 2.0 * ((SX[I]+Q[I]) - (PLEFT+Q[I+1])) /
          (RHO[I+1]*(X[I+1]-X[I]) + RHO[I]*(X[I]-X[I-1])) +
          2.0* (GEOM-1) * (PHI[I] - PHI[I+1]) /
          (RHO[I+1]*(X[I+1]+X[I]) + RHO[I]*(X[I-1]+X[I]));
          U[I] = U[I] + A[I] * (DT+DTOLD)/2.0;
          
        } else if (I < N) {
          A[I] = 2.0 * ((SX[I]+Q[I]) - (SX[I+1]+Q[I+1])) /
          (RHO[I+1]*(X[I+1]-X[I]) + RHO[I]*(X[I]-X[I-1])) +
          2.0* (GEOM-1) * (PHI[I] - PHI[I+1]) /
          (RHO[I+1]*(X[I+1]+X[I]) + RHO[I]*(X[I-1]+X[I]));
          U[I] = U[I] + A[I] * (DT+DTOLD)/2.0;
        } else {
          //  NON-REFLECTING BOUNDARY CONDITION (NRBC)
          //  BASED ON HALPERN, 1982
          double DXNRBC = X[N] - X[N-1];
          double TNRBC0 = 2.0*(pow(C[N-1]*DT,2)) / (C[N-1]*DT+DXNRBC);
          double TNRBC1 = (-1.0/DXNRBC) + (DXNRBC/pow(C[N-1]*DT,2));
          double TNRBC2 = (C[N-1]*DT-DXNRBC)/(2.0*pow(C[N-1]*DT, 2));
          double TNRBC3 = 1.0/DXNRBC;
          double UNEW = TNRBC0*(U[N]*TNRBC1+UNOLD*TNRBC2+U[N-1]*TNRBC3);
          UNOLD = U[N];
          U[N] = UNEW;
          A[N] = U[N] / DT;
        }
        XNEW[I] = X[I] + DT * U[I];
        //  COMPUTE NEW DENSITY FROM ZONE MASS
        DENGEOM(GEOM,XNEW[I-1],XNEW[I],RDEN);
        RHONEW[I] = ZMASS[I-1] / (XNEW[I] - XNEW[I-1] * RDEN);
		
        double RHOBAR = (RHONEW[I] + RHO[I]) / 2.0;
        DELRHO = 0.5 * (RHONEW[I]-RHO[I])/pow(RHOBAR,2);
        RHODOT = (RHONEW[I]-RHO[I])/(DT*RHOBAR);

        //  COMPUTE FACTORS FOR LINEAR (IN ENERGY) EOS
        double RHOEOS = RHONEW[I];
        if (I == ISURF) {
          //  SPECIAL TREATMENT FOR VAPORIZING SURFACE ZONE
          FSURF(RHOEOS,RHOLSURF,XQUAL[ISURF-1],QMIN,F1,F2);
        } else if (I > ISPALL) {
          //  SOLID/LIQUID ZONE
          EOSF(RHOEOS,F1,F2);
        } else if (VOLF[I-1] < 0.001) {
          EOSF(RHOEOS,F1,F2);
        } else if (VOLF[I-1] > 0.999) {
          FVAPOR(RHOEOS,F1,F2);
        } else {
          double XQ2 = XQUAL[I-1];
          F2PHASE(RHOEOS,HG[I-1],XQ2,F1,F2);
        }
        
        // OMPUTE ARTIFICIAL VISCOSITY

		double XXXX = (XNEW[I] - XNEW[I-1] + X[I] - X[I-1]) / 2.0;
        double QLIN = RHOBAR * B2 * XXXX * C[I-1] * RHODOT;
        double QQUAD = RHOBAR * pow(B1 * XXXX,2) * RHODOT * abs(RHODOT);
        //  CHOOSE VON NEUMAN & RICHTMEYER OR MOOS & WHITE FORM
        if (QFLAG == 2) {
          //  SET CRITERION FROM MOOS & WHITE
          QP = SX[I] * (U[I]-U[I-1])/XXXX;
          double QPTEST = QP1[I-1] + (QP1[I-1]-QP2[I-1]) * DT/DTOLD;
          double QF1 = 0.0;
          double QF2 = 0.0;
          if (QP > 0.0) QF1=1.0;
          if (QP > QPTEST) QF2=1.0;
          Q[I] = (QQUAD + QLIN) * (QF1 + QF2);
        } else {
          //  VON NEUMANN & RICHTMEYER
          if (RHODOT > 0.0) {
            Q[I] = QLIN + QQUAD;
          } else {
            Q[I] = 0.0;
          }
        }
        //  TURN OFF ART VISC FOR SURFACE ZONE
        if (I == (ISURF)) {
          Q[I] = 0.0;
        }
        
        double DELED;
        //  COMPUTE DEVIATORIC STRESSES AND ENERGIES IF SOLID
        //    IF (T(I) .LT. TMELT .AND. I .GT. ISPALL) THEN
        if (false) {
          // TODO dead code
          /*
           DXD = (U(I) - U(I-1)) / XXXX + RHODOT / 3.
           G = GNU * RHO(I) * C(I)**2
           //  COMPUTE CURRENT YIELD STRENGTH
           YIELDFLAG = YIELDCOEF(1)
           IF (YIELDFLAG .EQ. 1) THEN
           //  USING LINEAR RAMP OF YIELD STRENGTH FROM YIELD0 @ TINF DOWN TO 0 @ TMELT
           YIELD0 = YIELDCOEF(2)
           YBIG = YIELD0 * (TMELT-T(I))/(TMELT-TINF)
           ELSE IF (YIELDFLAG .EQ. 2) THEN
           //  USING EXPONENTIAL FIT
           CYIELD = YIELDCOEF(2)
           EYIELD = YIELDCOEF(3)
           YBIG = CYIELD * DEXP(EYIELD/T(I))
           END IF
           SIGMAXD = SXD(I) + 2. * DT * G * DXD
           
           //  SPLIT BASED ON GEOMETRY
           IF (GEOM .EQ. 1 .OR. GEOM .EQ. 3) THEN
           //  PLANAR OR SPHERICAL
           YIELDTEST = YBIG * 2./3.
           IF (ABS(SIGMAXD) .GT. YIELDTEST) THEN
           PLASTICFLAG = 1
           PLASTIC(I) = 1
           //          WRITE(*,*) 'PLASTICITY: ',I,TM
           IF (USEYIELDFLAG .EQ. 1) THEN
           SXDNEW(I) = ABS(SIGMAXD)/SIGMAXD * YIELDTEST
           ELSE
           SXDNEW(I) = SIGMAXD
           END IF
           ELSE
           SXDNEW(I) = SIGMAXD
           END IF
           PHI(I) = 1.5 * SXDNEW(I)
           DELED = 0.75*DT*DXD*(SXDNEW(I)+SXD(I)) / RHOBAR
           ELSE
           //  CYLINDRICAL
           DZD = RHODOT / 3.
           SIGMAZD = SZD(I) + 2. * DT * G* DZD
           FY = 2.*(SIGMAXD**2 + SIGMAXD*SIGMAZD + SIGMAZD**2)
           YIELDTEST = (2./3.) * YBIG**2
           IF (ABS(SIGMAXD) .GT. YIELDTEST) THEN
           PLASTICFLAG = 1
           PLASTIC(I) = 1
           IF (USEYIELDFLAG .EQ. 1) THEN
           FSQ = SQRT(YIELDTEST/FY)
           SXDNEW(I) = SIGMAXD * FSQ
           SZDNEW(I) = SIGMAZD * FSQ
           END IF
           ELSE
           FSQ = SQRT(YIELDTEST/FY)
           SXDNEW(I) = SIGMAXD * FSQ
           SZDNEW(I) = SIGMAZD * FSQ
           END IF
           DELED = (DT/(2.*RHOBAR)) * ((SXDNEW(I)+SXD(I))*(2.*DXD+DZD)+
           &         (SZDNEW(I)+SZD(I))*(2.*DZD+DZD))
           PHI(I) = 2.*SXDNEW(I) + SZDNEW(I)
           END IF
           */
        } else {
          //  MELTED OR SPALLED MATERIALS
		  DELED = 0.0;
          PHI[I] = 0.0;
          SXDNEW[I-1] = 0.0;
          SZDNEW[I-1] = 0.0;
        }
        if (I < ISPALL) SXDNEW[I-1] = 0.0;
        
        //  COMPUTE STRESS ENERGY CHANGE IN J/m2
        HSTRESS = (F1 + P[I] + 2.0*Q[I]) * DELRHO * ZMASS[I-1];

		//XUV
        //HSTRESSQ = (F1 + P[I]) * DELRHO * ZMASS[I-1];
        
        //  COMPUTE DEVIATORIC CONTRIBUTION IN J/m2
        HDEV = DELED * ZMASS[I-1];
      }
      //  END IF FOR LHYDRO-ONLY STEPS
      
      //  SET UP SOURCE ENERGY
      //  SUM TO GET NET ENERGY/ZONE
      //  NET ENERGY HAS UNITS OF J/m2
      double HIN = 0.0;
      for (int J=0; J<NBIN; J++) {
        HIN = HIN + EFRAC[I-1][J] * s_eintegral.EIN[J];
      }

	  //XUV
	  /*
      if (I <= ISURFEND) {
        EINMAXVAP = EINMAXVAP + HIN;
      } else {
        EINMAXSL = EINMAXSL + HIN;
      }
	  */
      
      //  HEAT CONDUCTION BALANCE ( (I-1) TO EAST = -(I TO WEST) )
	  double HWEST = HEAST;
      //  SET CONDUCTION COEFFICIENT FOR THIS ZONE
      double CONDP = CONDE;
      VFLAG = VOLF[I];
      if (I == (ISURF-1)) VFLAG = 0.0;
      CONDVST(T[I+1],VFLAG,CONDE);
      double CONDEI = 2.0 * CONDP * CONDE / (CONDP + CONDE);
      
      //  BASE HEAT CONDUCTION ON TEMP, LOCATIONS, PROPERTIES FROM PREVIOUS TM
      if (I == ISURF) {		
        //  ADJUST HEAT CONDUCTION DISTANCE TO REFLECT ONLY REMAINING LIQUID
        //  AND SOLVE IMPLICITLY
        double DXSURF = (X[I]-X[I-1]) * (1.0-VOLF[I-1]);
        double XINTERNAL = X[I] - DXSURF;
        BETAGEOM(GEOM,XINTERNAL,X[I],BETAC);
        double DXP1 = XC[I+1] - XC[I];
        //  FIGURE WHAT WILL BE HEAT COND BETWEEN (I+1) AND (I+2)
        double CONDE2;
        CONDVST(T[I+2],VOLF[I+1],CONDE2);
        double CONDEI1 = 2.0 * CONDE * CONDE2 / (CONDE + CONDE2);
        double BETAC1;
        BETAGEOM(GEOM,XC[I],XC[I+1],BETAC1);
        double AE1 = CONDEI1 * BETAC1 / DXP1;
        double HEAST1 = AE1*(T[I+1]-T[I+2]);
        //  SET TEMPORARY ENERGIES IN ZONE ISURF AND ISURF+1
		//XUV
        //  + HSTRESS
        double HTEMPP = H[ISURF-1] + HWEST * DT + HIN + HSTRESS;
        double HTEMPE = H[ISURF] - HEAST1*DT + HIN*EZONE[I]/EZONE[I-1];
        double HGTEMPP = HTEMPP / ZMASS[ISURF-1];
        double HGTEMPE = HTEMPE / ZMASS[ISURF];
        //  ENHANCE CONDUCTIVITY TO REDUCE SURFACE TEMPERATURE VARIATIONS
        CONDP = CONDP * (1.0+ 4.0*pow(XQUAL[ISURF-1],2));
        //  FILL COMMON BLOCKS FOR ITERATION ROUTINE
        s_surfc.SURFCSTUFF[0] = HTEMPP;
        s_surfc.SURFCSTUFF[1] = HTEMPE;
        s_surfc.SURFCSTUFF[2] = ZMASS[ISURF-1];
        s_surfc.SURFCSTUFF[3] = ZMASS[ISURF];
        s_surfc.SURFCSTUFF[4] = DT;
        s_surfc.SURFCSTUFF[5] = QMIN;
        s_surfc.SURFCSTUFF[6] = CONDP;
        s_surfc.SURFCSTUFF[7] = DXSURF;
        s_surfc.SURFCSTUFF[8] = BETAC;
        s_surfc.SURFCSTUFF[9] = T[I];
        s_surft.QUALSURFT = XQUAL[I-1];
        //  SET UP BRACKETING INITIAL GUESSES FOR HEAST BASED ON CURRENT TEMPS
        HEAST = CONDP * BETAC * (T[I]-T[I+1]) / DXSURF;
        
        double HE1,HE2;		
		BRACKET(CONDSURF, HEAST, HE1, HE2);
        double TOL = 1.0e-05;
        HEAST = ZBRENT(CONDSURF,HE1,HE2,TOL);
      } else {
        BETAGEOM(GEOM,XC[I-1],XC[I],BETAC);
        double DXC = XC[I] - XC[I-1];
        AE = CONDEI * BETAC / DXC;
        HEAST = AE*(T[I]-T[I+1]);
      }
      if (HEAST > HEASTMAX) HEASTMAX = HEAST;
      //  COMPUTE NET ENERGY CHANGE FROM CONDUCTION
      double HCOND = (HWEST - HEAST) * DT;
      
      //  UPDATE TOTAL ZONE AND SPECIFIC ENERGIES
	  //XUV
      //  + HSTRESS
      double DH = (HIN+HCOND+HDEV+HSTRESS) / (1.0 - F2*DELRHO);
      if (LHYDRO != 1) DH = HIN + HCOND;
      H[I-1] = H[I-1] + DH;
      HG[I-1] = H[I-1] / ZMASS[I-1];

	  //XUV
	  /*
      //  STORE EACH ENERGY CONTRIBUTION FOR PRINTING
      HINP[I] = HINP[I] + HIN / (1.0 - F2*DELRHO);
      HCONDP[I] = HCONDP[I] + HCOND / (1.0 - F2*DELRHO);
      HDEVP[I] = HDEVP[I] + HDEV / (1.0 - F2*DELRHO);
      HSTRESSP[I] = HSTRESSP[I] + HSTRESS / (1.0 - F2*DELRHO);
      HSTRESSQP[I] = HSTRESSQP[I] + HSTRESSQ / (1.0 - F2*DELRHO);      
      //  HSTRESSP(I) = (F1 + P(I) + 2.*Q(I)) * (0.5 / RHOBAR) * ZMASS(I) / (1.0D+00 - F2*DELRHO)
      */

      //  CHECK FOR MAXIMUM CHANGE IN ZONE ENERGY/kg
      double DHG = DH / ZMASS[I-1];
      if (abs(DHG) > DHGMAX) {
        DHGMAX = abs(DHG);
        int IDHGMAX = I;
        if (DHGMAX > DHGMXMAX) {
          DHGMXMAX = DHGMAX;
          int IDHGMXMAX = I;
        }
      }
      
      //  APPLY EOS INFORMATION TO GET PRESSURE, TEMP, SOUND SPEED
      //  DIFFERENT TREATMENTS DEPENDING ON EXPECTED STATE
      
      //  SPECIAL TREATMENT FOR SURFACE ZONE (PRESCRIBED QUALITY FROM SURFACE VAPORIZATION ROUTINE)
      if (I == ISURF) {
        s_surfp.TEMPSURF = T[I];
        s_surft.HGSURF = HG[I-1];
        s_surft.QUALSURFT = XQUAL[I-1];
        s_surfp.RHOSURF = RHONEW[I];
        s_surfp.QUALSURFP = XQUAL[I-1];
        double PM1 = P[I-1];
		
		EOSURF(RHOLSURF, RHOVSURF, RHONEW[I + 1], QMIN, PM1, P[I], VOLF[I - 1], C[I - 1]);

        TNEW[I] = s_surfp.TEMPSURF;
        if (TNEW[I] < 0.0) {
          cout << "ERROR IN SURF TEMP" << endl;
          cout << "T = " << TNEW[I] << endl;
          cout << "I = " << I << endl;
          cout << "TIME = " << TM*1.0e+09 << endl;
          cout << "HG (cal/g) = " << HT << endl;
          return 1;
        }
        
        //  USE SOLID/LIQUID EOS FOR UNSPALLED MATERIAL
      } else if (I > ISPALL) {
        //  DETERMINE WHAT THE NEW TEMPERATURE WILL BE
        HT = HG[I-1];
        H2T(HT,TNEW[I]);
        if (TNEW[I] < 0.0) {
          cout << "ERROR IN SOLID/LIQUID TEMP" << endl;
          cout << "T = " << TNEW[I] << endl;
          cout << "I = " << I << endl;
          cout << "TIME = " << TM*1.0e+09 << endl;
          cout << "HG (cal/g) = " << HT << endl;
          return 1;
        }
        
        //  COMPUTE NEW PRESSURE FROM EOS
        P[I] = F1 + F2 * HG[I-1];
        
        //  COMPUTE NEW SOUND SPEED
        double HEOS = HG[I-1];
        double REOS = RHONEW[I];
        double PEOS = P[I];
        double CEOS;
        SOUND(HEOS,REOS,PEOS,F2,CEOS,TM);
        if (TNEW[I] < s_matdat.TMELT) {
          //  APPLY CORRECTION FOR ELASTIC MATERIAL
          C[I-1] = SNU * CEOS;
        } else {
          C[I-1] = CEOS;
        }
        
        //  USE "SPALL" OR "VAPOR" EOS FOR SPALLED MATERIAL
      } else {
        double TS = T[I];
        if (HG[I-1] > 0.0 && HG[I-1] < 1.0e+15) {
          EOSPALL(HG[I-1],RHONEW[I],TS,P[I],XQUAL[I-1],VOLF[I-1],C[I-1]);
		  if(RHONEW[I] >= 1e5){
			  cout << 'T = ' << TNS << " I = " << I << " RHONEW =  " << RHONEW[I] << " RHO = " << RHO[I] << endl;
			  cout << "F1 = " << F1 << " F2 = " << F2 << " HG = " << HG[I-1] << "P = " << P[I] << endl;
		  }
          TNEW[I] = TS;
          if (TNEW[I] < 0.0) {
            cout << "ERROR IN SPALL TEMP" << endl;
            cout << "T = " << TNEW[I] << endl;
            cout << "I = " << I << endl;
            cout << "TIME = " << TM*1.0e+09 << endl;
            cout << "HG (cal/g) = " << HG[I-1] << endl;
            cout << "RHO = " << RHONEW[I] << endl;
            return 1;
          }
        } else {
          cout << "ERROR IN ZONE ENERGY" << endl;
          cout << "HG (cal/g) = " << HG[I-1]/4184 << endl;
          cout << "I = " << I << endl;
          cout << "TIME = " << TM*1.0e+09 << endl;
          cout << "T = " << T[I] << endl;
          cout << "RHO = " << RHONEW[I] << endl;
          cout << "QUAL = " << XQUAL[I-1] << endl;
          return 1;
        }
      }
      
      //  COMPUTE NEW SIGMA-X
      SXNEW[I] = P[I] - SXDNEW[I-1];
      
      //  COMPUTE MAXIMUM HYDRO TIME STEP FOR THIS ZONE
      double DELTAX = XNEW[I] - XNEW[I-1];	  
      double DTI = CT1 * DELTAX /
      (B2*C[I-1] + 2.0*pow(B1,2) * abs(RHODOT) * DELTAX +
       sqrt(pow(B2*C[I-1] + 2.0*pow(B1,2) * abs(RHODOT) * DELTAX, 2) +
            pow(C[I-1],2)));
      if (DTI < DTHYDRO) {
        DTHYDRO = DTI;
        double IHYDRO = I;
      }
      
      TNS = TM * 1.0e+09;
      
      //  UPDATE POWERS FOR ARTIFICIAL VISCOSITY SWITCH
      QP2[I-1] = QP1[I-1];
      QP1[I-1] = QP;
      
      if (LHYDRO == 1) {
        //  CHECK FOR SPALLED ZONE (LIQUID AT NEGATIVE STRESS)
        double TTEST = s_matdat.TMELT + 1.0e-05;
        if (TNEW[I] >= TTEST && I > ISPALL) {
          if (SXNEW[I] < SPALL) {
            ISPALL = I;
            cout << "ISPALL melt=" << ISPALL << endl;
            file21 << TM*1.0e+09 << TAB << ISPALL << TAB << X0[ISPALL]*1.0e+06 << endl;
            ISURF = ISPALL+1;
            ISPALLFLAG = 1;
          }
        }
      }

	  //XUV
	  /*
      F1P[I] = F1;
      F2P[I] = F2;
      DELRHOP[I] = DELRHO;
	  */
    }
    
    
    //  %%%%%%%%%%%%%%%%%%%%%%% END OF LOOP THROUGH ZONES %%%%%%%%%%%%%%%%%%%%%%%	
    //  UPDATE OLD T, X, XC,RHO, AND SIGMA VALUES
    X[0] = XNEW[0];
    T[0] = T[1];
    DTEMPMAX = 0.0;
    for (int i=1; i<=N; i++) {
      //  STORE FOR LIMITING TIME STEP
      DTEMP = abs(T[i]-TNEW[i]);
      if (DTEMP > DTEMPMAX) {
        DTEMPMAX = DTEMP;
        ITEMPMAX = i;
      }
      T[i] = TNEW[i];
      X[i] = XNEW[i];
      XC[i-1] = (X[i] + X[i-1]) / 2.0;
      RHO[i] = RHONEW[i];
      SX[i] = SXNEW[i];
      SXD[i-1] = SXDNEW[i-1];
      //  *** TO JSEM PRIDAL JA ***
      SZD[i-1] = SZDNEW[i-1];
    }
    T[N] = TINF;
    
    //  DETERMINE MELT DEPTH, IF ANY
    TMLOW = s_matdat.TMELT - 0.1;
    TMHIGH = s_matdat.TMELT + 0.1;
    if (T[1] >= TMLOW) {
      for (int I=ISPALL+1; I<=N; I++) {
        if (T[I] < TMLOW && T[I-1] >= TMLOW) {
          XMELT = XLOC[I-1];
		  ITMELT = I-1;
        }
        if (T[I] < TMHIGH && T[I-1] >= TMHIGH) {
          XTMELT = XLOC[I-1];
		  ITMELT = I-1;
        }
      }
      if (XMELT > XMELTMAX) XMELTMAX = XMELT;
      if (XTMELT > XTMELTMAX) {
        XTMELTMAX = XTMELT;
      }
    } else {
      XMELT = 0.0;
      XTMELT = 0.0;
	  ITMELT = 0.0;
    }
    
    //  SURFACE VAPORIZATION SECTION
    //  COMPUTE VAPORIZED MASS AND ENERGY IN THIS TIME STEP
    //  FIRST SET AMBIENT PRESSURE
	//  PERFORM ONLY IF SURFACE IS LIQUID
	if(TNEW[ISURF] > 1.1*s_matdat.TMELT){	
		if (ISURF == 1) {
		  PAMB = 1.33e-04;
		} else {
		  PAMB = P[ISURF-1];
		}
		//DTVAP = 100.0;
		SURFVAP(GEOM,PAMB,TNEW[ISURF],XNEW[ISURF],DT,DVAPMASS,DHVAP);
    
		//  GET NEW QUALITY (MASS FRACTION OF VAPOR) IN SURFACE ZONE
		DXQ = DVAPMASS / ZMASS[ISURF-1];
		XQ = XQUAL[ISURF-1] + DXQ;
    
		if (XQ < 1.0) {
		  XQUAL[ISURF-1] = XQ;
		} else {
		  //  SPALL OFF THIS ZONE AND ADVANCE ISURF
		  DQI = 1.0 - XQUAL[ISURF];
		  XQUAL[ISURF-1] = 1.0;
		  VOLF[ISURF-1] = 1.0;
		  DVAPMASS = ZMASS[ISURF-1] * (XQ - 1.0);
		  ISPALL = ISURF;
		  ISURF = ISURF + 1;
		  cout << "ISURFvapor =" << ISPALL << endl;
		  cout << "TIME =" << TM*1.0 << endl;
		  cout << "CYCLE =" << NCYCLE << endl;
		  sprintf(line_buffer, "%F10.4,%c,%i5,%c,%F10.4", TM*1.0e+09,TAB,ISPALL,TAB,X0[ISPALL]*1.0e+06);
		  file21 << line_buffer << endl;
      
		  //  SET UP NEXT SURFACE ZONE PROPERLY
		  DXQ = DVAPMASS / ZMASS[ISURF-1];
		  XQUAL[ISURF-1] = DXQ;
		  RHOV = P[ISURF] / (RBAR * T[ISURF]);
		  //  *** TO JSEM ZMENIL JA (VOLF = QUAL * RHO / RHOV) ***
		  VOLF[ISURF-1] = XQUAL[ISURF-1] * RHO[ISURF] / RHOV;
		  //    VOLF = QUAL * RHO / RHOV
      
		  RHOLSURF = RHO[ISURF];
		  ISPALLFLAG = 1;
		}
	}
    //  END OF SURFACE VAPORIZATION SECTION
    
    //  SET NEW TIME STEP TO PREVENT EXCESSIVE VAPORIZATION NEXT STEP
	DTVAP = 100.0;
    if (ISPALLFLAG == 1) {
      SURFVAP(GEOM,PAMB,TNEW[ISURF],XNEW[ISURF],DT,DVMNEXT,DHVAPN);
      VAPRATIO = DVMNEXT / ZMASS[ISURF-1];
      if (VAPRATIO > 0.02) {
        DTVAP = DT * 0.02 / VAPRATIO;
      }
    }
    
    //  KEEP TRACK OF MAXIMUM FRONT SURFACE ENERGY
    if (HG[0] > HGFMAX) {
      HGFMAX = HG[0];
    }
    
    //  *** TO JSEM PRIDAL JA, PRI HG(UP) MUZE T(DOWN) (VYPAROVANI) ***
    //  KEEP TRACK OF MAXIMUM FRONT SURFACE TEMPERATURE
    if (T[1] > TFMAX) {
      TFMAX = T[1];
    }
    
    //  WRITE SURFACE ZONE PRESSURE TO A FILE FOR A WHILE
    TNS = TM * 1.0e+09;
    if (TNS < 3.0) {
      sprintf(line_buffer, "%F12.7,%c,%i5,%c,%E12.5,%c,%E12.5,%c,%E12.5", TNS,TAB,ISURF,TAB,P[ISURF],TAB,T[ISURF],TAB,T[ISURF+1]);
      file10 << line_buffer << endl;
    }
    
    //  CHECK ENERGY BALANCE
    HSUM = 0.0;
    EKIN = 0.0;
    for (int I=1; I<=N; I++) {
	  HSUM = HSUM + H[I-1];
	  //XUV
      //    HSUM = HSUM + H(I) - HSTRESS(I)
      //HSUM = HSUM + H[I] - HSTRESS;
      UBAR = (U[I] + EKINSIGN*U[I-1]) / 2.0;
      EKIN = EKIN + 0.5 * ZMASS[I-1] * pow(UBAR,2);
    }
    ESUM = HSUM + EKIN;
    
    //  OUTPUT FOR DIAGNOSTIC PURPOSES
    NCYMOD = NCYCLE % NCYCLEPRINT;
	
    if (NCYMOD == 0) {
      TNS = TM * 1.0e+09;
      cout << "TIME" << TM*1.0e+09 << endl;
      K1 = ISPALL + 1;
      K2 = ISPALL + NSAVE;
      cout << "U    ";
      for (int K=K1; K<= K2; K++) {
        sprintf(line_buffer, "%E12.5", 10*U[K]);
        cout << line_buffer << " ";
      }
      cout << endl;
      cout << "P    ";
      for (int K=K1; K<= K2; K++) {
        sprintf(line_buffer, "%E12.5", 10*P[K]);
        cout << line_buffer << " ";
      }
      cout << endl;
      cout << "T    ";
      for (int K=K1; K<= K2; K++) {
        sprintf(line_buffer, "%E12.5", 10*T[K]);
        cout << line_buffer << " ";
      }
      cout << endl;
      sprintf(line_buffer, "Energy present & input (J/cm2) = %F12.3 %F12.3", ESUM/10000.0, EINMAX/10000.0);
      cout << line_buffer << endl;
    }
    NCYCLE = NCYCLE + 1;
    
    //  DETERMINE NEW HEAT TIME STEP FROM CFL CONSTRAINT
    DTHEAT = 1.0e+02;
    for (int I=1; I<=N; I++) {
      if (T[I] < 0.0) {
        cout << "cfl T =" << T[I] << endl;
        cout << "I =" << I << endl;
        cout << "ISURF=" << ISURF << endl;
        cout << "CYCLE=" << NCYCLE << endl;
        cout << "HG =" << HG << endl;
        return 1;
      }
      VFLAG = VOLF[I-1];
      if (I == ISURF) VFLAG = 0.0;
      CONDVST(T[I],VFLAG,COND);
      HEATCAP(VOLF[I-1],T[I],HG[I-1],CV);
      ALPHA = COND / (RHO[I] * CV);	  
      DXHEAT = X[I] - X[I-1];
      DTI = CFL * pow(DXHEAT,2) / ALPHA;
      if (DTI < DTHEAT) {
        DTHEAT = DTI;
      }
    }
    
    //  LIMIT TIME STEP IF MAX TEMP CHANGE IN LAST TIME STEP WAS TOO BIG
    TEMPBIG = 200;
    if (DTEMPMAX > TEMPBIG) {
      DTTEMP = DT * TEMPBIG / DTEMPMAX;
    } else {
      DTTEMP = 100.0;
    }
    
    //  SET NEXT TIME STEP FROM HYDRO AND THERMAL CONSTRAINTS	
    DTOLD = DT;
	
    DT = min5(DTHEAT,DTHYDRO,CT2*DT,DTVAP,DTTEMP);
    if (LHYDRO != 1) DT = DTHEAT;
    if (DT > DTMAX) DTMAX = DT;
    if (DT < DTMIN) {
      DTMIN = DT;
    }
    
    //  TIME STEP DIAGNOSTIC SECTION
    if (DT < 1.0e-17) {
      sprintf(line_buffer, "low DT %E12.4 %E12.4 %E12.4 %E12.4 %E12.4 %E12.4 %E12.4", DT,DTHEAT,DTHYDRO,DTVAP,DTTEMP,DTEMPMAX, ITEMPMAX);
      cout << line_buffer << endl;
      //    WRITE(*,*) ITEMPMAX,T(ITEMPMAX),ABS(T(ITEMPMAX) - TNEW(ITEMPMAX))
    }
    
    //  COMPUTE NUCLEATION RATE (LOG10(J) = RATEJ)
    for (int K=1; K<=N; K++) {
      if (K > ISURF && T[K] > s_matdat.TMELT) {
        NUCLEATE(T[K],RHO[K],P[K],RAWJ,RAWTNUC);
        
        RATEJ[K-1] = log10(RAWJ);
        TNUC[K-1] = log10(RAWTNUC);
        RAWBUB = pow(10,BUBBLE[K-1]) + RAWJ*DT;
        BUBBLE[K-1] = log10(RAWBUB);
        
        //  TEST VARIOUS MAXIMUM THRESHOLD DEPTHS
        if (RATEJ[K-1] > 35.0 && K > I35) {
          I35 = K;
          X35 = XLOC[K-1];
        }
        if (RATEJ[K-1] > 30.0 && K > I30) {
          I30 = K;
          X30 = XLOC[K-1];
        }
        if (RATEJ[K-1] > 25.0 && K > I25) {
          I25 = K;
          X25 = XLOC[K-1];
        }
        if (RATEJ[K-1] > 20.0 && K > I20) {
          I20 = K;
          X20 = XLOC[K-1];
        }
        if (BUBBLE[K-1] > 27.0 && K > IB27) {
          IB27 = K;
          XB27 = XLOC[K-1];
        }
        if (BUBBLE[K-1] > 26.0 && K > IB26) {
          IB26 = K;
          XB26 = XLOC[K-1];
        }
        if (BUBBLE[K-1] > 25.0 && K > IB25) {
          IB25 = K;
          XB25 = XLOC[K-1];
        }
        if (BUBBLE[K-1] > 24.0 && K > IB24) {
          IB24 = K;
          XB24 = XLOC[K-1];
        }
      } else {
        RATEJ[K-1] = -99.0;
      }
    }
    
	//XUV
	/*
    //  WRITE PRESSURE AND STUFF TO FILES EVERY DTPRESS NANOSECONDS
    TNS = TM * 1.0e+09;
    
    //  WRITE THE PULSE ENERGY TO FILE
    DTPRESSEIN=1e-5;
    if (s_pulse.TSQ > 0.0) {
      TPULSE = s_pulse.TSQ + 0.04 * s_pulse.TSQ;
    } else if (s_pulse.FWHM > 0.0) {
      TPULSE = s_pulse.FWHM * 3;
    }
    if (TNS < TPULSE) {
      if (TNS > TPRESSEIN) {
        sprintf(line_buffer, "%F10.5 %c %E11.5", TNS,TAB,EINDT*(DTPRESSEIN/(DTOLD*1.0e+09))/10000.0);
        file32 << line_buffer << endl;
        TPRESSEIN=TPRESSEIN+DTPRESSEIN;
      }
    }
    */

    if (TNS < 0.001) {
      DTPRESS = 1E-4;
    } else if (TNS < 0.01) {
      DTPRESS = 5E-4;
      //  ELSE IF (TNS .LT. 0.1) THEN
      //    DTPRESS = 5E-3;
    } else if (TNS < 0.5) {
	  //XUV
      //DTPRESS = 1E-3;
		DTPRESS = 5E-3;
    } else if (TNS < 1.0) {
      DTPRESS = 0.01;
    } else if (TNS < 3.0) {
      DTPRESS = 0.05;
    } else if (TNS < 20.0) {
      DTPRESS = 0.5;
    } else {
      DTPRESS = 5.0;
    }
    
    if (TNS > TPRESS) {
      sprintf(line_buffer, "%F10.5", TNS);
      file11 << line_buffer;
      for (int K=1; K<= N; K++) {
        sprintf(line_buffer, "%E12.5", P[K]);
        file11 << TAB << line_buffer;
      }
      sprintf(line_buffer, "%F10.5", TNS);
      file13 << line_buffer;
      for (int K=1; K<= N; K++) {
        sprintf(line_buffer, "%E12.5", SX[K]);
        file13 << TAB << line_buffer;
      }
      sprintf(line_buffer, "%F10.5", TNS);
      file14 << line_buffer;
      for (int K=1; K<= N; K++) {
		sprintf(line_buffer, "%E12.5", C[K-1]);
		//XUV
        //sprintf(line_buffer, "%E12.5", SXD[K]);
        file14 << TAB << line_buffer;
      }
      sprintf(line_buffer, "%F10.5", TNS);
      file27 << line_buffer;
      for (int K=1; K<= N; K++) {
        sprintf(line_buffer, "%E12.5", Q[K]);
        file27 << TAB << line_buffer;
      }
      sprintf(line_buffer, "%F10.5", TNS);
      file26 << line_buffer;
      for (int K=0; K<= N; K++) {
        sprintf(line_buffer, "%E12.5", A[K]);
        file26 << TAB << line_buffer;
      }
      sprintf(line_buffer, "%F10.5", TNS);
      file15 << line_buffer;
      for (int K=0; K<= N; K++) {
        sprintf(line_buffer, "%E12.5", U[K]);
        file15 << TAB << line_buffer;
      }
      sprintf(line_buffer, "%F10.5", TNS);
      file16 << line_buffer;
      for (int K=0; K<= N; K++) {
        sprintf(line_buffer, "%E12.5", XNEW[K]);
        file16 << TAB << line_buffer;
      }
      sprintf(line_buffer, "%F10.5", TNS);
      file17 << line_buffer;
      for (int K=1; K<= N; K++) {
        sprintf(line_buffer, "%E12.5", RHONEW[K]);
        file17 << TAB << line_buffer;
      }
      sprintf(line_buffer, "%F10.5", TNS);
      file18 << line_buffer;
      for (int K=1; K<= N; K++) {
        sprintf(line_buffer, "%E12.5", XQUAL[K-1]);
        file18 << TAB << line_buffer;
      }
      sprintf(line_buffer, "%F10.5", TNS);
      file19 << line_buffer;
      for (int K=1; K<= N; K++) {
        sprintf(line_buffer, "%E12.5", T[K]);
        file19 << TAB << line_buffer;
      }
      sprintf(line_buffer, "%F10.5", TNS);
      file20 << line_buffer;
      for (int K=1; K<= N; K++) {
        sprintf(line_buffer, "%E12.5", HG[K-1]);
        file20 << TAB << line_buffer;
      }
      sprintf(line_buffer, "%F10.5", TNS);
      file22 << line_buffer;
      for (int K=1; K<= N; K++) {
        sprintf(line_buffer, "%E12.5", VOLF[K-1]);
        file22 << TAB << line_buffer;
      }
      sprintf(line_buffer, "%F10.5", TNS);
      file24 << line_buffer;
      for (int K=1; K<= N; K++) {
        sprintf(line_buffer, "%E12.5", TNUC[K-1]);
        file24 << TAB << line_buffer;
      }
      sprintf(line_buffer, "%F10.5", TNS);
      file25 << line_buffer;
      for (int K=1; K<= N; K++) {
        sprintf(line_buffer, "%E12.5", RATEJ[K-1]);
        file25 << TAB << line_buffer;
      }
      sprintf(line_buffer, "%F10.5", TNS);
      file35 << line_buffer;
      for (int K=1; K<= N; K++) {
        sprintf(line_buffer, "%E12.5", BUBBLE[K-1]);
        file35 << TAB << line_buffer;
      }

	  //CALCULATE RT INSTABILITY PARAMS
	  double ZMASS1 = 0;
	  double A1 = A[ISPALL + 1];
	  T1 = 0;
	  double DMELT = 0;
	  //cout << ITMELT << endl;

	  for(int K=ISPALL+2;K<=ITMELT;K++){
		  ZMASS1 = ZMASS1 + ZMASS[K-1];
		  A1 = A1 + A[K];
		  double ZONEDIAM = X[K] - X[K-1];
		  DMELT = DMELT + ZONEDIAM;
		  T1 = T1 + T[K] * ZONEDIAM;
	  }

	  double REYNOLDS;
	  double WEBER;

	  if(DMELT > 0.0){
		  double RHO1 = ZMASS1 + DMELT;
		  A1 = A1 / (ITMELT - (ISPALL+1));
		  double F = A1 * RHO1;
		  T1 = T1 / DMELT;
		  //ETA UNITS: KG.M-1.S-1 (=1000 CENTIPOISE)
		  double ETA = 5.5e-04;
		  double SIGMA = SIGMACALC(T1);
		  if(F > 0.0 && SIGMA > 0){
			  REYNOLDS = DMELT* sqrt(F*RHO1*DMELT) / ETA;
			  WEBER = (F*DMELT*DMELT)/SIGMA;
		  }
		  else{
			  REYNOLDS = 0;
			  WEBER = 0;
		  }
		  file40 << TNS << "\t" << DMELT << "\t" << RHO1 << "\t" << F << "\t" << SIGMA << "\t" << T1 << "\t" << A1 << endl;
	  }

	  //XUV
	  /*
      sprintf(line_buffer, "%F10.5", TNS);
      file28 << line_buffer;
      for (int K=1; K<= N; K++) {
        sprintf(line_buffer, "%i1", PLASTIC[K]);
        file28 << TAB << line_buffer;
      }
      sprintf(line_buffer, "%F10.5", TNS);
      file37 << line_buffer;
      for (int K=1; K<= N; K++) {
        sprintf(line_buffer, "%E12.5", HINP[K]/10000.0);
        file37 << TAB << line_buffer;
      }
      sprintf(line_buffer, "%F10.5", TNS);
      file38 << line_buffer;
      for (int K=1; K<= N; K++) {
        sprintf(line_buffer, "%E12.5", HCONDP[K]/10000.0);
        file38 << TAB << line_buffer;
      }
      sprintf(line_buffer, "%F10.5", TNS);
      file39 << line_buffer;
      for (int K=1; K<= N; K++) {
        sprintf(line_buffer, "%E12.5", HDEVP[K]/10000.0);
        file39 << TAB << line_buffer;
      }
      sprintf(line_buffer, "%F10.5", TNS);
      file40 << line_buffer;
      for (int K=1; K<= N; K++) {
        sprintf(line_buffer, "%E12.5", HSTRESSP[K]/10000.0);
        file40 << TAB << line_buffer;
      }
      sprintf(line_buffer, "%F10.5", TNS);
      file41 << line_buffer;
      for (int K=1; K<= N; K++) {
        sprintf(line_buffer, "%E12.5", F1P[K]);
        file41 << TAB << line_buffer;
      }
      sprintf(line_buffer, "%F10.5", TNS);
      file42 << line_buffer;
      for (int K=1; K<= N; K++) {
        sprintf(line_buffer, "%E12.5", F2P[K]);
        file42 << TAB << line_buffer;
      }
      sprintf(line_buffer, "%F10.5", TNS);
      file43 << line_buffer;
      for (int K=1; K<= N; K++) {
        sprintf(line_buffer, "%E12.5", ZMASS[K]/6*(pow(U[K],2)+pow(U[K-1],2) + U[K]*U[K-1]));
        file43 << TAB << line_buffer;
      }
      sprintf(line_buffer, "%F10.5", TNS);
      file44 << line_buffer;
      for (int K=1; K<= N; K++) {
        sprintf(line_buffer, "%E12.5", (1/6)*ZMASS[K]*(pow(U[K],2) + U[K]*U[K-1] + pow(U[K-1],2)));
        file44 << TAB << line_buffer;
      }
      */

	  /*
      //  3RD AND 2ND T DERIVATIVE RATIO FOR DIAGNOSTICS
      DX = X[N] - X[N-1];
      RELDTLIMIT = 0.001;
      //  FIRST DERIVATIVE
      for (int I=2; I<=N-1; I++) {
        TMV = (T[I+1] + T[I-1]) / 2;
        if (abs((T[I+1]-T[I-1]) / TMV) > RELDTLIMIT) {
          TP[I] = (T[I+1] - T[I-1]) / (2 * DX);
	  
      } else {
          TP[I] = 0.0;
        }
      }
      TMV = (T[2] + T[1]) / 2;
      if (abs((T[2]-T[1]) / TMV) > RELDTLIMIT) {
        TP[1] = (T[2] - T[1]) / DX;
      } else {
        TP[1] = 0.0;
      }
      TMV = (T[N] + T[N-1]) / 2;
      if (abs((T[N]-T[N-1]) / TMV) > RELDTLIMIT) {
        TP[N-1] = (T[N] - T[N-1]) / DX;
      } else {
        TP[N-1] = 0.0;
      }
      
      //  SECOND DERIVATIVE
      for (int I=2; I<=N-1; I++) {
        TPP[I] = (TP[I+1] - TP[I-1]) / (2 * DX);
      }
      TMV = (T[2] + T[1]) / 2;
      if (abs((T[2]-T[1]) / TMV) > RELDTLIMIT) {
        TPP[1] = (T[2] - T[1]) / pow(DX,2);
      } else {
        TPP[1] = 0.0;
      }
      TMV = (T[N] + T[N-1]) / 2;
      if (abs((T[N]-T[N-1]) / TMV) > RELDTLIMIT) {
        TPP[N-1] = (-T[N] + T[N-1]) / pow(DX,2);
      } else {
        TPP[N-1] = 0.0;
      }
      
      //  THIRD DERIVATIVE
      for (int I=2; I<=N-1;I++) {
        TPPP[I] = (TPP[I+1] - TPP[I-1]) / (2 * DX);
      }
      TPPP[1] = TPP[1] / DX;
      TPPP[N-1] = -TPP[N-1] / DX;
      
      //  RATIO
      for (int I=1; I<=N; I++) {
        if (abs(TPP[I]) > 0) {
          RATIO[I] = abs((TPPP[I] * DX) / (2 * TPP[I]));
        } else {
          RATIO[I] = 0.0;
        }
      }
      
      sprintf(line_buffer, "%F10.5", TNS);
      file29 << line_buffer;
      for (int K=1; K<= N; K++) {
        sprintf(line_buffer, "%E12.5", RATIO[K]);
        file29 << TAB << line_buffer;
      }
      sprintf(line_buffer, "%F10.5", TNS);
      file33 << line_buffer;
      for (int K=1; K<= N; K++) {
        sprintf(line_buffer, "%E12.5", TP[K]);
        file33 << TAB << line_buffer;
      }
      sprintf(line_buffer, "%F10.5", TNS);
      file30 << line_buffer;
      for (int K=1; K<= N; K++) {
        sprintf(line_buffer, "%E12.5", (2 * TPP[K]));
        file30 << TAB << line_buffer;
      }
      sprintf(line_buffer, "%F10.5", TNS);
      file31 << line_buffer;
      for (int K=1; K<= N; K++) {
        sprintf(line_buffer, "%E12.5", (TPPP[K] * DX));
        file31 << TAB << line_buffer;
      }
      
      //  ENERGY BALANCE
      //  SUM FOR DIFFERENT TYPES + INPUT ENERGY
      //  FOR VAPORS
      HINSUM = 0.0;
      HCONDSUM = 0.0;
      HDEVSUM = 0.0;
      HSTRESSSUM = 0.0;
      EKIN = 0.0;
      for (int I=1; I<=ISURFEND; I++) {
        HINSUM = HINSUM + HINP[I];
        HCONDSUM = HCONDSUM + HCONDP[I];
        HDEVSUM = HDEVSUM + HDEVP[I];
        HSTRESSSUM = HSTRESSSUM + HSTRESSP[I];
        UBAR = (pow(U[N],2) + U[N]*U[N-1] + pow(U[N-1],2))/6;
        EKIN = EKIN + ZMASS[I] * UBAR;
      }
      sprintf(line_buffer, "%F10.5 %c %F10.6 %c %F10.6 %c %F10.6 %c %F10.6 %c %F10.6", TNS,TAB,HINSUM/10000.0,TAB,HCONDSUM/10000.0,TAB,HDEVSUM/10000.0,TAB,HSTRESSSUM/10000.0,TAB,EKIN/10000.0);
      file47 << line_buffer << endl;
      //  FOR SOLID/LIQUID
      HINSUM = 0.0;
      HCONDSUM = 0.0;
      HDEVSUM = 0.0;
      HSTRESSSUM = 0.0;
      EKIN = 0.0;
      for (int I=ISURFEND+1; I<=N; I++) {
        HINSUM = HINSUM + HINP[I];
        HCONDSUM = HCONDSUM + HCONDP[I];
        HDEVSUM = HDEVSUM + HDEVP[I];
        HSTRESSSUM = HSTRESSSUM + HSTRESSP[I];
        UBAR = (pow(U[N],2) + U[N]*U[N-1] + pow(U[N-1],2))/6;
        EKIN = EKIN + ZMASS[I] * UBAR;
      }
      sprintf(line_buffer, "%F10.5 %c %F10.6 %c %F10.6 %c %F10.6 %c %F10.6 %c %F10.6", TNS,TAB,HINSUM/10000.0,TAB,HCONDSUM/10000.0,TAB, HDEVSUM/10000.0,TAB,HSTRESSSUM/10000.0,TAB,EKIN/10000.0);
      file48 << line_buffer << endl;
      //  VAPOR + SOLID/LIQUID
      HINSUM = 0.0;
      HCONDSUM = 0.0;
      HDEVSUM = 0.0;
      HSTRESSSUM = 0.0;
      EKIN = 0.0;
      for (int I=1; I<=N; I++) {
        HINSUM = HINSUM + HINP[I];
        HCONDSUM = HCONDSUM + HCONDP[I];
        HDEVSUM = HDEVSUM + HDEVP[I];
        HSTRESSSUM = HSTRESSSUM + HSTRESSP[I];
        UBAR = (pow(U[N],2) + U[N]*U[N-1] + pow(U[N-1],2))/6;
        EKIN = EKIN + ZMASS[I] * UBAR;
      }
      sprintf(line_buffer, "%F10.5 %c %F10.6 %c %F10.6 %c %F10.6 %c %F10.6 %c %F10.6",TNS,TAB,HINSUM/10000.0,TAB,HCONDSUM/10000.0,TAB,HDEVSUM/10000.0,TAB,HSTRESSSUM/10000.0,TAB,EKIN/10000.0);
      file34 << line_buffer << endl;
      
      //  SUM OF ALL TYPES + INPUT ENERGY
      //  FOR VAPORS
      ESUM = 0.0;
      for (int I=1; I<=ISURFEND; I++) {
        UBAR = (U[I] + EKINSIGN*U[I-1]) / 2.0;
        EKIN = 0.5 * ZMASS[I] * pow(UBAR,2);
        ESUM = ESUM + H[I] - HSTRESSP[I];
      }
      sprintf(line_buffer, "%F10.5 %c %F10.6 %c %F10.6", TNS,TAB,EINMAXVAP/10000.0,TAB,ESUM/10000.0);
      file45 << line_buffer << endl;
      //  FOR SOLID/LIQUID
      ESUM = 0.0;
      for (int I=ISURFEND+1; I<=N; I++) {
        UBAR = (U[I] + EKINSIGN*U[I-1]) / 2.0;
        EKIN = 0.5 * ZMASS[I] * pow(UBAR,2);
        ESUM = ESUM + H[I] - HSTRESSP[I];
      }
      sprintf(line_buffer, "%F10.5 %c %F10.6 %c %F10.6", TNS,TAB,EINMAXSL/10000.0,TAB,ESUM/10000.0);
      file46 << line_buffer << endl;
      //  VAPOR + SOLID/LIQUID
      ESUM = 0.0;
      for (int I=1; I<=N; I++) {
        UBAR = (U[I] + EKINSIGN*U[I-1]) / 2.0;
        EKIN = 0.5 * ZMASS[I] * pow(UBAR,2);
        ESUM = ESUM + H[I] - HSTRESSP[I];
      }
      sprintf(line_buffer, "%F10.5 %c %F10.6 %c %F10.6", TNS,TAB,EINMAX/10000.0,TAB,ESUM/10000.0);
      file46 << line_buffer << endl;
      
      //  CHECK THE MOMENTUM CONSERVATION
      MOMSUM = 0.0;
      for (int I=1; I<=N; I++) {
        MOMSUM = MOMSUM + (U[I]+U[I-1]/2) * ZMASS[I];
      }
      sprintf(line_buffer, "%F10.5 %c %E11.5", TNS,TAB,MOMSUM);
      file49 << line_buffer << endl;

	  */
      
      TPRESS = TPRESS + DTPRESS;
    }
    
    //  WRITE SOME INFO TO SCREEN EVERY FEW NS
    DTSCREEN = 2.0;
    if (TNS > 50.0) DTSCREEN = 10.0;
    if (TNS > TSCREEN) {
      K1 = ISURF;
      K2 = ISURF + 4;
      cout << TNS << " " << DT << " ";
      for (int K=K1; K<=K2; K++) {
        cout << T[K] << " ";
      }
      cout << endl;
      TSCREEN = TSCREEN + DTSCREEN;
    }
    
  } // end of marching through zones
  
  //  *****************************************************************
  //  TIME LIMIT HAS BEEN REACHED
  
  file10.close();
  file11.close();
  file13.close();
  file14.close();
  file15.close();
  file16.close();
  file17.close();
  file26.close();
  file27.close();
  file18.close();
  file19.close();
  file20.close();
  file21.close();
  file22.close();
  file25.close();
  file35.close();

  ifs2.close();
  ifs3.close();
  ifs4.close();
  ifs5.close();
  ifs6.close();
  ifs7.close();

  //XUV
  /*
  file28.close();
  file29.close();
  file30.close();
  file31.close();
  file32.close();
  file33.close();
  file34.close();
  file36.close();
  file37.close();
  file38.close();
  file39.close();
  */
  file40.close();
  /*
  file41.close();
  file42.close();
  file43.close();
  file44.close();
  file45.close();
  file46.close();
  file47.close();
  file48.close();
  file49.close();
  */

  X0M1 = X0[ISURF-1];
  VAPDEPTH = 1.0e+06*(X0M1+XQUAL[ISURF-1]*(X0[ISURF]-X0M1));
  file23.open("summary");
  EINMAX = EINMAX / 10000.0;
  sprintf(line_buffer, "Total fluence (J/cm2) = %F12.3", EINMAX);
  cout << line_buffer << endl;
  sprintf(line_buffer, "Maximum melt depth (micron) = %F12.3", XMELTMAX);
  cout << line_buffer << endl;
  sprintf(line_buffer, "Maximum total melt depth (micron) = %F12.3", XTMELTMAX);
  cout << line_buffer << endl;
  sprintf(line_buffer, "Spalled depth (micron) = %F12.6", X0[ISPALL] * 1.0e+06);
  cout << line_buffer << endl;
  sprintf(line_buffer," Maximum front surface temp = %F12.2",TFMAX);
  cout << line_buffer << endl;
  sprintf(line_buffer,"Maximum front surface enthalpy (cal/g) = %F12.4",HGFMAX / 4184.0);
  cout << line_buffer << endl;
  sprintf(line_buffer,"Minimum and maximum time steps = %E12.3 %E12.3",10*DTMIN, 10*DTMAX);
  cout << line_buffer << endl;
  sprintf(line_buffer,"Maximum enthalpy change (cal/g) = %F12.4",DHGMXMAX / 4184.0);
  cout << line_buffer << endl;
  sprintf(line_buffer,"Maximum heat flux (W/cm2) = %E12.4",10*HEASTMAX / 10000.0);
  cout << line_buffer << endl;
  sprintf(line_buffer,"Fraction vaporized in surface zone %F10.5",XQUAL[ISURF-1]);
  cout << line_buffer << endl;
  sprintf(line_buffer,"Energy present & input (J/cm2) = %F12.3 %F12.3",ESUM/10000.0, EINMAX);
  cout << line_buffer << endl;
  sprintf(line_buffer,"Vaporized depth (micron) = %F12.3",VAPDEPTH);
  cout << line_buffer << endl;
  sprintf(line_buffer,"J-depths (35,30,25,20) = %F10.3 %F10.3 %F10.3 %F10.3",X35,X30,X25,X20);
  cout << line_buffer << endl;
  sprintf(line_buffer,"Bubble-depths (27,26,25,24) = %F10.3 %F10.3 %F10.3 %F10.3",XB27,XB26,XB25,XB24);
  cout << line_buffer << endl;

  //XUV
  /*
  if (PLASTICFLAG == 1) {
    cout << "Plastic deformation has occured!" << endl;
  }
  */
  
  sprintf(line_buffer, "Total fluence (J/cm2) = %F12.3", EINMAX);
  file23 << line_buffer << endl;
  sprintf(line_buffer, "Maximum melt depth (micron) = %F12.3", XMELTMAX);
  file23 << line_buffer << endl;
  sprintf(line_buffer, "Maximum total melt depth (micron) = %F12.3", XTMELTMAX);
  file23 << line_buffer << endl;
  sprintf(line_buffer, "Spalled depth (micron) = %F12.3", X0[ISPALL] * 1.0e+06);
  file23 << line_buffer << endl;
  sprintf(line_buffer, "Maximum front surface temp = %F12.2", TFMAX);
  file23 << line_buffer << endl;
  sprintf(line_buffer, "Maximum front surface enthalpy (cal/g) = %F12.4", HGFMAX / 4184.0);
  file23 << line_buffer << endl;
  sprintf(line_buffer, "Minimum and maximum time steps = %E12.3 %E12.3", 10*DTMIN, 10*DTMAX);
  file23 << line_buffer << endl;
  sprintf(line_buffer, "Maximum enthalpy change (cal/g) = %F12.4", DHGMXMAX / 4184.0);
  file23 << line_buffer << endl;
  sprintf(line_buffer, "Maximum heat flux (W/cm2) = %E12.4", 10*HEASTMAX / 10000.0);
  file23 << line_buffer << endl;
  sprintf(line_buffer, "Fraction vaporized in surface zone %F10.5", XQUAL[ISURF-1]);
  file23 << line_buffer << endl;
  sprintf(line_buffer, "Energy present & input (J/cm2) = %F12.3 %F12.3", ESUM/10000.0, EINMAX);
  file23 << line_buffer << endl;
  sprintf(line_buffer, "Vaporized depth (micron) = %F12.3", VAPDEPTH);
  file23 << line_buffer << endl;
  sprintf(line_buffer, "J-depths (35,30,25,20) = %F10.3 %F10.3 %F10.3 %F10.3", X35,X30,X25,X20);
  file23 << line_buffer << endl;
  sprintf(line_buffer, "Bubble-depths (27,26,25,24) = %F10.3 %F10.3 %F10.3 %F10.3", XB27,XB26,XB25,XB24);
  file23 << line_buffer << endl;

  //XUV
  /*
  if (PLASTICFLAG == 1) {
    cout << "Plastic deformation has occured!" << endl;
  }
  */
  
  file23.close();
  
  return 0;
  //  THIS IS THE END OF ABLATOR ROUTINE
}

// **********************************************************
// DEFINE THE SPECTRUM FUNCTION
double SPECTRUM(int JBB)
{		
	double SPECTRUM;
	if(BIN[JBB] <= 101.225){
		SPECTRUM = 54.83 * BIN[JBB];
	}
	else{
		SPECTRUM = BINW[JBB] * 
    			(479.64453*(1240/BIN[JBB] -0.01) +
    			25190.299*exp(-pow(1240/BIN[JBB] - 1.1547276,2) /pow(0.3190276,2)) +
    			15228.539*exp(-pow(1240/BIN[JBB] - 1.395925,2) /pow(0.13990381,2)) +
     			9344.5156*exp(-pow(1240/BIN[JBB] - 1.7601666,2) /pow(0.21062079,2)) +
     			4744.1201*exp(-pow(1240/BIN[JBB] - 2.2334344,2) /pow(0.47139305,2)) +
     			3530.1509*exp(-pow(1240/BIN[JBB] - 3.1980155,2) /pow(0.63428837,2)) +
     			5704.8501*exp(-pow(1240/BIN[JBB] - 4.7492857,2) /pow(1.4205953,2)) +
     			2281.9927*exp(-pow(1240/BIN[JBB] - 7.3327885,2) /pow(2.3396976,2)));
	}

	return SPECTRUM;
}

//  **********************************************************
//  READ IN MATERIAL PROPERTY DATA
//  UNITS INPUT AS CONVENIENT, BUT CONVERT TO SI IN THIS SUBROUTINE
void READMAT()
{
  ifstream ifs("matdata");
  MATNAME = FortranIO::read_string(ifs);
  cout << MATNAME << endl;
  FortranIO::skip_line(ifs);
  //  DENSITY IN kg/m**3
  s_eosdat.RHO0 = FortranIO::read_double(ifs);
  FortranIO::skip_line(ifs);
  //  MELT TEMPERATURE IN K
  s_matdat.TMELT = FortranIO::read_double(ifs);
  FortranIO::skip_line(ifs);
  //  CURVE FIT COEFFICIENT FOR THERMAL COND VS TEMP
  //  FIFTH-ORDER POLYNOMIALS FOR SOLID AND LIQUID
  //  AND ONE PARAMETER FOR VAPOR
  //  UNITS: cal/m.s.K and K
  for(int i=0;i<6;i++){
    s_matdat.ACOND[i] = FortranIO::read_double(ifs);
  }
  FortranIO::skip_line(ifs);
  for(int i=6;i<12;i++)
  {
    s_matdat.ACOND[i] = FortranIO::read_double(ifs);
  }
  FortranIO::skip_line(ifs);
  //xuv
  //s_matdat.ACOND[12] = FortranIO::read_double(ifs);
  for(int i=12;i<16;i++){
	   s_matdat.ACOND[i] = FortranIO::read_double(ifs);
  }
  //  CONVERT TO SI UNITS (W/m.K)
  for(int i=0;i<16;i++)
  {
    s_matdat.ACOND[i] *= 4.184;
  }
  //  CURVE FIT COEFFICIENTS FOR TEMPERATURE VS ENTHALPY
  //  FIRST 6 FOR SOLID, THEN 2 EACH FOR LIQUID AND VAPOR
  //  UNITS: K and cal/g
  FortranIO::skip_line(ifs);
  for(int i=0;i<6;i++)
  {
    s_matdat.AH2T[i] = FortranIO::read_double(ifs);
  }
  FortranIO::skip_line(ifs);
  for(int i=6;i<8;i++)
  {
    s_matdat.AH2T[i] = FortranIO::read_double(ifs);;
  }
  FortranIO::skip_line(ifs);
  for(int i=8;i<10;i++)
  {
    s_matdat.AH2T[i] = FortranIO::read_double(ifs);
  }
  //  CONVERT TO SI UNITS (J/kg)
  for(int i=1;i<6;i++)
  {
    s_matdat.AH2T[i] = s_matdat.AH2T[i] / pow(4184.0,i);
  }
  s_matdat.AH2T[7] /= 4184.0;
  s_matdat.AH2T[9] /= 4184.0;
  //  Tboil vs. P COEFFICIENTS
  //  FROM DUSHMAN: LOG(P) = A - B/T
  //  UNITS: P in Pa, T in K
  FortranIO::skip_line(ifs);
  s_eosdat.ABOIL = FortranIO::read_double(ifs);
  s_eosdat.BBOIL = FortranIO::read_double(ifs);
  //  CONVERT TO SI UNITS (Pa)
  //s_eosdat.ABOIL -= log10(7.5);
  //  EQUATION OF STATE COEFFICIENTS
  FortranIO::skip_line(ifs);
  for(int i=0;i<10;i++)
  {
    s_eosdat.EOSCOEF[i] = FortranIO::read_double(ifs);
  }
  //  POISSON'S RATIO
  FortranIO::skip_line(ifs);
  s_eosdat.PRNU = FortranIO::read_double(ifs);
  //  YIELD STRENGTH
  FortranIO::skip_line(ifs);
  for(int i=0;i<3;i++)
  {
    s_eosdat.YIELDCOEF[i] = FortranIO::read_double(ifs);
  }
  //  SURFACE TENSION PARAMETERS
  FortranIO::skip_line(ifs);
  for(int i=0;i<7;i++)
  {
    s_surften.SIGMACOEF[i] = FortranIO::read_double(ifs);
  }
  ifs.close();
  
  //  COMPUTE ETHALPIES FOR INCIPIENT ANT TOTAL MELT
  //  START BY ASSUMING A LINEAR FIT FOR LIQUID T VS. H
  //  UNITS ARE: J/kg and K
  s_matdat.HMELTT = (s_matdat.TMELT - s_matdat.AH2T[6]) / s_matdat.AH2T[7];
  //  NOW GET INCIPIENT MELT ENERGY VIA NEWTON'S METHOD
  //double hm = 0.6 * s_matdat.HMELTT;
  double hm = s_matdat.HMELTT * 0.6;
  double hmxx = s_matdat.HMELTT * 3;
  double ah2tzero = s_matdat.AH2T[0] - s_matdat.TMELT;
  for(int i=0;i<6;i++)
  {
    hm -= (ah2tzero + hm*(s_matdat.AH2T[1] + hm*(s_matdat.AH2T[2] + hm*(s_matdat.AH2T[3] + (hm*(s_matdat.AH2T[4] + hm*s_matdat.AH2T[5])))))) /
          (s_matdat.AH2T[1] + hm*(2.0*s_matdat.AH2T[2] + hm*(3.0*s_matdat.AH2T[3] + hm*(4.0*s_matdat.AH2T[4] + hm*5.0*s_matdat.AH2T[5]))));
  }
  s_matdat.HMELTI = hm;
  cout << std::setprecision(12) << "H-MELT-I (CAL/G) = " << hm/4184.0 << endl;
  cout << "H-MELT-I (CAL/G) = " << s_matdat.HMELTT/4184.0 << endl;
}

//  **********************************************************
//  SET ZONE EDGE LOCATIONS, ZONE CENTERS AND ZONE MASSES
//  GEOM IS FLAG FOR PLANAR (1), CYLINDRICAL (2) OR SPHERICAL (3)
void ZONESET(double& rho0,double x0[], double xc[], double xloc[],double zmass[],int& geom, ifstream& file)
{
  //  PLANAR GEOMETRY
  geom = 1;
  //  FIRST ZONE GRID SIZE DX METERS
  //  DX = 5E-10 
  
  FortranIO::skip_line(file);  
  double dx = FortranIO::read_double(file);
  dx *= 1e-9;

  //XUV
  //double ratio = 1.02;
  double ratio = 1.08;

  //  THEESE PARAMETERS GIVE TOTAL X = 56 MICRONS FOR 100 ZONES
  //  ZONE MASS IN kg/m2
  //  ZONE LOCATIONS IN meters
  //  INITIAL ZONE CENTER LOCATIONS IN XLOC ARRAY IN microns
  x0[0] = 0.0;
  for(int i=0;i<NP1;i++){
    x0[i+1] = x0[0] + dx * (pow(ratio,i+1) - 1.0) / (ratio-1.0);
    DENGEOM(geom,x0[i],x0[i+1],RDEN);
    zmass[i] = rho0 * (x0[i+1] - x0[i]) * RDEN;	
    xc[i] = (x0[i+1] + x0[i]) / 2.0;
    xloc[i] = xc[i] * 1.0e6;
  }
  //??? TODO -
  //130 FORMAT(A)
}

// **********************************************************
// RETURN GEOMETRY FACTOR TO GET CORRECT ZONE MASSES
void DENGEOM(int& geom, double& r1, double& r2, double& rden)
{
  if(geom == 1)
    rden = 1.0;
  else if(geom == 2)
    rden = PI * (r1 + r2);
  else
    rden = (4.0/3.0) * PI * (pow(r1,2.0) + r1*r2 + pow(r2,2.0));
}

/*C **********************************************************
 C RETURN GEOMETRY FACTOR TO GET CORRECT HEAT CONDUCTION
 C Q = BETAC * K * DT / DX
 C FOR W/m2 planar, W/m cylindrical and W spherical*/
void BETAGEOM(int& geom, double& r1, double& r2, double& betac)
{
  if(geom == 1)
    betac = 1.0;
  else if(geom == 2)
    betac = PI * (r1 + r2);
  else
    betac = PI * pow((r1 + r2),2.0);
}

/*C **********************************************************
 C RETURN AREA AT A GIVEN RADIUS DEPENDING ON GEOMETRY*/
void AREAGEOM(int& geom, double& r, double& area)
{
  if(geom == 1)
    area = 1.0;
  else if(geom == 2)
    area = 2.0 * PI * r;
  else
    area = 4.0 * PI * pow(r,2.0);
}

/*C **********************************************************
 C FIND NEW TIME STEP SIZE TO HAVE X0 POINT JUST HIT AXIS
 C FOR CYLINDRICAL AND SPHERICAL GEOMETRIES
 C ** CAUTION!! UNTESTED ROUTINE **/
void AXIS(double& x0, double& u0, double& a0, double& dtold, double& dtnew)
{
  double a = a0/2.0;
  double b = u0 + a0*dtold/2.0;
  double c = x0;
  
  double b24ac = pow(b,2.0) - 4*a*c;
  if(b24ac <= 0)
    cout << "ERROR IN ASXIS ROUTINE" << endl;
  
  double dtp = (-b + sqrt(b24ac)) / (2.0*a);
  double dtm = (-b - sqrt(b24ac)) / (2.0*a);
  
  //unreferenced
  //double dtmin;
  
  //C CHOOSE NEW TIME STEP AS SMALLEST POSITIVE ROOT
  if(dtm <= 0.0){
    if(dtp >= 0.0)
      dtnew = dtp;
    else
      cout << "ERROR IN AXIS ROUTINE" << endl;
  }
  else{
    DTMIN = dtm;
  }
}

/*C **********************************************************
 C SET UP FOR ENERGY DEPOSITION BASED ON ISOURCE CHOICE*/
void XSOURCE(int& isource, int& laserflag, double fbbleak[], double& tstart,ifstream& file9)
{
  double bbt_leak;
  //in fortran version in this place the black body function bbe was redefined
  //it makes no sense to do such thing in c++
  if(isource == 1){
    FortranIO::skip_line(file9);
    s_pulse.TSQ = FortranIO::read_double(file9);
    FortranIO::skip_line(file9);
    s_pulse.ESQ = FortranIO::read_double(file9);
    
    if(laserflag != 1){
      FortranIO::skip_line(file9);
      bbt_leak = FortranIO::read_double(file9);
    }
    else{
      // SKIP READING BLACKBODY TEMPERATURE FROM THE FILE
      FortranIO::skip_line(file9);
      FortranIO::skip_line(file9);
    }
  }
  else if(isource == 2){
    FortranIO::skip_line(file9);
    s_pulse.FWHM = FortranIO::read_double(file9);
    FortranIO::skip_line(file9);
    s_pulse.EGAUSS = FortranIO::read_double(file9);
    if(laserflag != 1){
      FortranIO::skip_line(file9);
      bbt_leak = FortranIO::read_double(file9);
    }
    else{
      // SKIP READING BLACKBODY TEMPERATURE FROM THE FILE
      FortranIO::skip_line(file9);
      FortranIO::skip_line(file9);
    }
  }
  else{
    cerr << "Invalid source term number entered" << endl;
    exit(EXIT_FAILURE);
  }
  
  if(laserflag != 1){
    // FRACTION OF BLACKBODY ENERGY IN EACH BIN
    // ASSUME SQUARE AND GAUSSIAN PULSES ARE PURE LAMBERTIAN
    double sumleak = 0.0;
    for(int j=0;j<NBIN;j++){
		fbbleak[j] = SPECTRUM(j);
      //fbbleak[j] = BBE(j,bbt_leak);
      sumleak += fbbleak[j];
    }
    // NORMALIZE TO MAKE THESE FRACTIONS OF TOTAL ENERGY
    for(int j=0;j<NBIN;j++){
      fbbleak[j] /= sumleak;
    }
  }
  else{
    // FOR LASER USE MONOCHROMATIC SPECTRUM
    for(int j=0;j<NBIN;j++){
      fbbleak[j] = 0.0;
    }
    fbbleak[0] = 1.0;
  }
  // SET PROBLEM START TIME, TSTART IN ns
  TSTART = 0.0;
}

/*
 SET MAGNITUDE OF ENERGY PULSES
 RETURNS ENERGY IN EACH BIN FOR THIS TIME STEP
 */
void XPULSE(double& tm, int& isource, double& area, double fbbleak[45])
{
  double tns = tm * 1.0e9;
  double ftmleak;
  if(isource == 1){
    // SQUARE PULSE, STARTS AT TIME = 0
    // USE 2% RAMP AT START AND END
    double tramp = s_pulse.TSQ * 0.02;
    if(tns < tramp){
      ftmleak = s_pulse.ESQ * pow(tns,2.0) / (2.0 * s_pulse.TSQ * tramp);
    }
    else if(tns <= s_pulse.TSQ){
      ftmleak = s_pulse.ESQ * (tns-tramp/2.0) / s_pulse.TSQ;
    }
    else if(tns < (s_pulse.TSQ+tramp)){
      double tau = tns - s_pulse.TSQ;
      ftmleak = (s_pulse.ESQ/s_pulse.TSQ)*(s_pulse.TSQ - tramp/2.0 + tau - pow(tau,2.0) / (2.0*tramp));
    }
    else{
      ftmleak = s_pulse.ESQ;
    }
  }
  else{
    // GAUSSIAN PULSE
    double tstar = (tns - 1.5 * s_pulse.FWHM) * s_pulse.SC / sqrt(2.0);
    cout << tns << endl << s_pulse.SC << endl << tstar << endl;
    double a1 = 0.278393;
    double a2 = 0.230389;
    double a3 = 0.000972;
    double a4 = 0.078108;
    double xin = tstar;
    double xg = abs(xin);
    cout << xg << endl;
    
    double myerf = abs(1.0 - 1.0 / pow(1.0 + xg * (a1 + xg * (a2 + xg * (a3 + xg * a4))),4.0));
    cout << "MYERF(TSTAR) = " << myerf << "(" << xin << ")" << endl;
    ftmleak = 0.5 * (0.0 + myerf) * s_pulse.EGAUSS;
    cout << "FTMLEAK = " << ftmleak << endl;
  }
  // SET UP INCIDENT ENERGY IN EACH BIN IN THIS TIME STEP
  // EIN(J) IS IN J/m2; ELEAK ARE IN JOULES
  for(int i=0;i<NBIN;i++){
    s_eintegral.ELEAK[i] = fbbleak[i] * ftmleak;
    double elnet = s_eintegral.ELEAK[i] - s_eintegral.ELEAKOLD[i];
    s_eintegral.EIN[i] = elnet / area;
    s_eintegral.ELEAKOLD[i] = s_eintegral.ELEAK[i];
  }
}

/**********************************************************
 C DEFINE THERMAL CONDUCTIVITY VS TEMPERATURE
 C USE UP TO A FIFTH-ORDER POLYNOMIAL
 C UNITS ARE: W/m.K and K IN CURVE FITS
 */
void CONDVST(double& t,double& volf, double& cond)
{
  if(t < 0.0){
    cerr << "T = "  << t << endl;
    exit(EXIT_FAILURE);
  }
  
  double conds = s_matdat.ACOND[0] + t * (s_matdat.ACOND[1] + t * (s_matdat.ACOND[2] + t * (s_matdat.ACOND[3] + t * (s_matdat.ACOND[4] + t * s_matdat.ACOND[5]))));
  double condl = s_matdat.ACOND[6] + t * (s_matdat.ACOND[7] + t * (s_matdat.ACOND[8] + t * (s_matdat.ACOND[9] + t * (s_matdat.ACOND[10] + t * s_matdat.ACOND[11]))));
  double condv = s_matdat.ACOND[12] + t * (s_matdat.ACOND[13] + t * (s_matdat.ACOND[14] + t * (s_matdat.ACOND[15])));  
  
  double T2K = 20000;

  //XUV
  // MONOATOMIC GAS CONDUCTIVITY
  //double condv = s_matdat.ACOND[12] * sqrt(t);

  if(t > T2K){
	  condl = s_matdat.ACOND[6] + T2K * (s_matdat.ACOND[7] + T2K * (s_matdat.ACOND[8] + T2K * (s_matdat.ACOND[9] + T2K * (s_matdat.ACOND[10] + T2K * s_matdat.ACOND[11]))));
	  condv = s_matdat.ACOND[12] + T2K * (s_matdat.ACOND[13] + T2K * (s_matdat.ACOND[14] + T2K * (s_matdat.ACOND[15])));  
  }
  
  if(volf <= 0.1){
    if (t < (s_matdat.TMELT - 0.1))
      cond = conds;
    else if(t > (s_matdat.TMELT + 0.1))
      cond = condl;
    else
      // TWO PHASE LIQUID/SOLID REQION
      cond = 2.0 * conds * condl / (conds + condl);
  }
  else if(volf >= 0.99){
    cond = condv;
  }
  else{
    // TWO PHASE LIQUID/VAPOR REGION
    cond = condv;
  }
}

/*
 C CONVERT ENTHALPY TO A TEMPERATURE
 C FOR BOTH SOLID AND LIQUID
 C USE UP TO A FIFTH-ORDER POLYNOMIALS
 C UNITS ARE: J/kg and K IN CURVE FITS
 */
void H2T(double hjkg, double& t)
{
  double h = hjkg;
  double ts = s_matdat.AH2T[0] + h*(s_matdat.AH2T[1] + h*(s_matdat.AH2T[2] + h*(s_matdat.AH2T[3] + h*(s_matdat.AH2T[4] + h*s_matdat.AH2T[5]))));
  double tl = s_matdat.AH2T[6] + h*s_matdat.AH2T[7];
  
  if(tl > s_matdat.TMELT)
    t = tl;
  else if(ts < s_matdat.TMELT)
    t = ts;
  else
    t = s_matdat.TMELT;
}


/*
 C DETERMINE HEAT CAPACITY (Cv) FOR SOLID, LIQUID AND VAPOR
 C VAPOR IS BASED ON IDEAL GAS, USSING GAMMA & RBAR
 C LIQUID AND SOLID USE T vs. ENTHALPY RELATIONS(Cp~Cv)
 C DECIDE STATE BASED ON VOLUME FRACTION VAPOR (VOLF)
 C USE UP TO A FIFTH-ORDER POLYNOMIALS
 C UNITS ARE: J/kg and K IN CURVE FITS
 */
void HEATCAP(double& volf, double& t, double& hg, double& cv)
{
  double& gamma = s_eosdat.EOSCOEF[8];
  double& rbar = s_eosdat.EOSCOEF[9];
  
  double cv_vapor = rbar / (gamma - 1.0);
  double cv_liquid = 1.0 / s_matdat.AH2T[7];
  double cv_solid = 1.0 / (s_matdat.AH2T[1] + hg*(2.0*s_matdat.AH2T[2] + hg*(3.0*s_matdat.AH2T[3] + hg*(4.0*s_matdat.AH2T[4] + hg*(5.0*s_matdat.AH2T[5])))));
  
  if(volf > 0.999999)
    // VAPOR
    cv = cv_vapor;
  else if(volf < 0.001){
    if(t > s_matdat.TMELT){
      // LIQIUD
      cv = cv_liquid;
    }
    else if( t < s_matdat.TMELT){
      // SOLID
      cv = cv_solid;
    }
    else{
      // TWO-PHASE SOLID/LIQUID
      if(cv_liquid < cv_solid)
        cv = cv_liquid;
      else
        cv = cv_solid;
    }
  }
  else{
    // TWO-PHASE VAPOR/LIQUID
    if(cv_liquid < cv_vapor)
      cv = cv_liquid;
    else
      cv = cv_vapor;
  }
}

/*
 C APPROXIMATION TO ERF(X) FROM ABRAMOWITZ & STEGUN
 */
double ERF(double& xin)
{
  const double  A1 = 0.278393,A2 = 0.230389,A3 = 0.000972,A4 = 0.078108;
  double erf;
  
  double x = abs(xin);
  if(xin >= 0)
    erf = abs(1.0 - 1.0 / pow(1.0 + x *(A1 + x * (A2 + x * (A3 + x * A4))),4.0));
  else
    erf = -1;
  
  return erf;
}

/*
 C COMPUTE EQUATION OF STATE COEFFICIENTS
 C F1 = F1(RHO) AND F2 = F2(RHO)
 C FOR LINEAR EOS: PRESSURE = F1 + F2 * ENERGY
 */
void EOSF(double& rho, double& f1, double& f2)
{
  double u = rho/s_eosdat.RHO0 - 1.0;
  // EOS TYPE
  int itype = int(s_eosdat.EOSCOEF[0]);
  double c0,s1,s2,s3,g0,b;
  double c[7];
  
  if(itype == 1){
    // STIENBERG'S MODEL (DYNA2D EOS #4)
    c0 = s_eosdat.EOSCOEF[1];
    s1 = s_eosdat.EOSCOEF[2];
    s2 = s_eosdat.EOSCOEF[3];
    s3 = s_eosdat.EOSCOEF[4];
    g0 = s_eosdat.EOSCOEF[5];
    b = s_eosdat.EOSCOEF[6];
    if(u <= 0.0)
      f1 = s_eosdat.RHO0 * pow(c0,2.0) * u;
    else
      f1 = s_eosdat.RHO0 * pow(c0,2.0) * u *
      (1.0 + (1.0 -g0/2.0) * u + b/2.0 * pow(u,2.0)) /
      (1.0 - (s1 - 1.0) * u - s2 * pow(u,2.0)/(u+1.0) - s3*pow(u,3)/pow(u+1.0,2.0));
    f2 = (g0 + b * u)*rho;
  }
  else if(itype == 2){
    // LINEAR POLYNOMIAL MODEL (DYNA2D EOS #1 AND #6)
    for(int i=0;i<=6;i++){
      //TODO: zkontrolovat, ze to tady pocita spravne, protoze jsem snizil index v eoscoeff o jedna
      c[i] = s_eosdat.EOSCOEF[i+1];
    }
    
    if( u <= 0){
      f1 = c[0] + c[1]*u + c[3]*pow(u,3.0);
      f2 = (c[4] + c[5]*u) * rho;
    }
    else{
      f1 = c[0] + c[1]*u + c[2]*pow(u,2.0) + c[3]*pow(u,3.0);
      f2 = (c[4] + c[5]*u + c[6]*pow(u,2.0)) * rho;
    }
  }
  else{
    cerr << "NOT A VALID EOS NUMBER" << endl;
    exit(EXIT_FAILURE);
  }
}

/*
 C COMPUTE EQUATION OF STATE COEFFICIENT DERIVATIVES
 C F1 = F1(RHO) AND F2 = F2(RHO)
 C FOR LINEAR EOS: PRESSURE = F1 + F2 * ENERGY
 */
void EOSDF(double& RHO,double& df1, double& df2)
{
  double U = RHO/s_eosdat.RHO0 - 1.0;
  double& RHO0 = s_eosdat.RHO0;
  // EOS TYPE
  int itype = int(s_eosdat.EOSCOEF[0]);
  double C0,S1,S2,S3,G0,B,xnumer,xdenom;
  //unreferenced
  //double c[7];
  
  if(itype == 1){
    // STIENBERG'S MODEL (DYNA2D EOS #4)
    C0 = s_eosdat.EOSCOEF[1];
    S1 = s_eosdat.EOSCOEF[2];
    S2 = s_eosdat.EOSCOEF[3];
    S3 = s_eosdat.EOSCOEF[4];
    G0 = s_eosdat.EOSCOEF[5];
    B = s_eosdat.EOSCOEF[6];
    if(U <= 0.0)
      df1 = pow(C0,2.0);
    else{
      xnumer = pow(C0,2.0)*RHO*(2*pow(RHO,3.0) + 4.*(1.-G0/2.)*pow(RHO,3.0)*U +
                                3.*B*pow(RHO,3.0)*pow(U,2.0) - 2.*(1.-G0/2.)*(S1-1.)*pow(RHO,3.0)*pow(U,2.0) +
                                2.*pow(RHO,2.0)*RHO0*S2*pow(U,2.0) - 2.*B*(S1-1.)*pow(RHO,3.0)*pow(U,3.0) -
                                2.*RHO*pow(RHO0,2.0)*S2*pow(U,3.0) + 4.*RHO*pow(RHO0,2.0)*S3*pow(U,3.0) -
                                B*pow(RHO,2)*RHO0*S2*pow(U,4.0) -
                                2.*(1-G0/2.)*RHO*pow(RHO0,2.0)*S2*pow(U,4.0) +
                                2.*(1-G0/2.)*RHO*pow(RHO0,2.0)*S3*pow(U,4.0) -
                                4.*pow(RHO0,3.0)*S3*pow(U,4.0) - B*RHO*pow(RHO0,2.0)*S2*pow(U,5.0) -
                                4.*(1.-G0/2.)*pow(RHO0,3.0)*S3*pow(U,5.0) - 2.*B*pow(RHO0,3.0)*S3*pow(U,6.0));
      xdenom = 2.*pow((-pow(RHO,2.0) + (S1-1.)*pow(RHO,2.0)*U +
                       RHO*RHO0*S2*pow(U,2.0) + pow(RHO0,2.0)*S3*pow(U,3.0)),2.0);
      df1 = xnumer/xdenom;
    }
    df2 = G0 - B + 2. * B * RHO / RHO0;
  }
  else if(itype == 2){
    // LINEAR POLYNOMIAL MODEL (DYNA2D EOS #1 AND #6)
    for(int k=0;k<=6;k++){
      C[k] = s_eosdat.EOSCOEF[k+1];
    }
    
    if(U <= 0.0){
      df1 = (C[1] + 3.*C[3]*pow(U,2.0)) / RHO0;
      df2 = C[4] - C[5] + 2. * C[5] * RHO / RHO0;
    }
    else{
      df1 = (C[1] + 2.*C[2]*U + 3.*C[3]*pow(U,2.0)) / RHO0;
      df2 = C[4] - C[5] + 2.*C[5]*RHO/RHO0 + C[6]*(1.0+RHO/RHO0*(-4.0e+3*RHO/RHO0));
    }
  }
  else{
    cerr << "NOT A VALID EOS NUMBER" << endl;
    exit(EXIT_FAILURE);
  }
}

/*
 C COMPUTE EQUATION OF STATE COEFFICIENTS FOR PURE VAPOR
 C F1 = F1(RHO) AND F2 = F2(RHO)
 C FOR LINEAR EOS: PRESSURE = F1 + F2 * ENERGY
 C ASSUMES IDEAL GAS EOS: P = RHO * RBAR * T
 */
void FVAPOR(double& RHO,double& F1, double& F2)
{
  double& RHO0 = s_eosdat.RHO0;
  double& AV0 = s_matdat.AH2T[8];
  double& AV1 = s_matdat.AH2T[9];
  double& GAMMA = s_eosdat.EOSCOEF[8];
  
  F1 = 0.0;
  F2 = (GAMMA - 1.0) * RHO;
}

/*
 C COMPUTE SOUND SPEED IN SOLID/LIQUID
 C USING DERIVATIVES OF EQUATION OF STATE COEFFICIENTS
 C F1 = F1(RHO) AND F2 = F2(RHO)
 C FOR LINEAR EOS: PRESSURE = F1 + F2 * ENERGY
 */
void SOUND(double& HG, double& RHO, double& P, double& F2, double& C, double& TM)
{
	double DF1,DF2;
	EOSDF(RHO,DF1,DF2);

	double DET = DF1 + HG*DF2 + P*F2/pow(RHO,2.0);
	if(DET >= 0){
		C = sqrt(DET);
	}
	else{
		cout << "'PROBLEM WITH SOUND SPEED: DET='" << DET << endl;
	}
}

/*
 C EQUATION OF STATE COEFFICIENT FOR SPALLED MATERIAL (LIQUID AND/OR VAPOR)
 C P = F1(RHO) + F2(RHO) * HG
 C FOR LIQUID PART, USE SOLID (UNSPALLED) EOS
 C FOR VAPOR PART, USE IDEAL GAS RELATION
 C BASE TWO-PHASE FITTING ON Tsat vs. Psat CURVE
 C UNITS: AL,AV COEFFICIENTS IN J/kg, PRESSURE IN Pa
 */
void F2PHASE(double& RHO,double& HG,double& QUAL, double& F1, double& F2)
{
  // CALL EOS ROUTINE TWICE WITH SLIGHTLY DIFFERENT ENERGIES
  double HG1 = HG;
  double DELTAHG = HG * 0.01;
  double HG2 = HG + DELTAHG;
  
  double T,P1,P2,VOLF,CSOUND;
  EOSPALL(HG1,RHO,T,P1,QUAL,VOLF,CSOUND);
  EOSPALL(HG2,RHO,T,P2,QUAL,VOLF,CSOUND);
  F2 = (P2-P1) / DELTAHG;
  F1 = P1 - F2 * HG1;
}

/*
 C EQUATION OF STATE FOR SPALLED MATERIAL (LIQUID AND/OR VAPOR)
 C ALSO COMPUTE SOUND SPEED
 C FOR LIQUID PART, USE SOLID (UNSPALLED) EOS
 C FOR VAPOR PART, USE IDEAL GAS RELATION
 C BASE TWO-PHASE FITTING ON Tsat vs. Psat CURVE
 C NOTE THAT QUAL (OLD VALUE) IS SENT AS AN INITIAL GUESS
 C UNITS: AL,AV COEFFICIENTS IN J/kg, PRESSURE IN Pa
 */
void EOSPALL(double& HG,double& RHO,double& T,double& P,double& QUAL,double& VOLF,double& CSOUND)
{
  double& ABOIL = s_eosdat.ABOIL;
  double& BBOIL = s_eosdat.BBOIL;
  double& HGTWO = s_qtwodat.HGTWO;
  double& RHOTWO = s_qtwodat.RHOTWO;
  double& RHOV = s_qtwodat.RHOV;
  double& RHOL = s_qtwodat.RHOL;
  double& TTWO = s_qtwodat.TTWO;
  double& HGL = s_qtwodat.HGL;
  // MINIMUM PRESSURE IN CHAMBER IS 10^-6 TORR
  double PMIN = 1.33e-4;
  double& GAMMA = s_eosdat.EOSCOEF[8];
  double& RBAR = s_eosdat.EOSCOEF[9];
  
  double& TINF = s_matdat.AH2T[0];
  double& AL0 = s_matdat.AH2T[6];
  double& AL1 = s_matdat.AH2T[7];
  double& AV0 = s_matdat.AH2T[8];
  double& AV1 = s_matdat.AH2T[9];
  // FIRST SEE IF LIQUID STATE WORKS OUT FOR LOW-QUALITY ZONE
  double TL = AL0 + HG * AL1;
  double F1,F2;
  EOSF(RHO,F1,F2);
  double PL = F1 + F2 * HG;
  double PSAT = pow(10.0,ABOIL-BBOIL/TL);
  if (PL > PSAT && QUAL < 0.05){
    // PURE LIQUID
    P = PL;
    T = TL;
	double CAS = 0.0;
    SOUND(HG,RHO,PL,F2,CSOUND, CAS);
    QUAL = 0.0;
    VOLF = 0.0;
  }
  else{
    // SEE IF PURE VAPOR WORKS FOR HIGH-QUALITY ZONE
    double TV = AV0 + HG * AV1;
    double PV = RHO * RBAR * TV;
    PSAT = pow(10.0,ABOIL-BBOIL/TL);
    if(PV < PSAT && QUAL > 0.95){
      // PURE VAPOR WORKS
      T = TV;
      P = PV;
      CSOUND = sqrt(RBAR*GAMMA*T);
      QUAL = 1.0;
      VOLF = 1.0;
    }
    else{
      // TWO-PHASE REGION
      // SEE IF TWO-PHASE BASED ON CURRENT QUALITY WILL WORK
      HGTWO = HG;
      RHOTWO = RHO;
      double QLOW = 0.0 < (QUAL-0.05) ? (QUAL-0.05) : 0.0;
      double QHIGH = 1.0 < (QUAL+0.05) ? 1.0 : (QUAL+0.05);
      double TESTLOW = QTWO(QLOW);
      double TESTHIGH = QTWO(QHIGH);
      double TESTPROD = TESTLOW * TESTHIGH;
      if(TESTPROD < 0.0){
        // EQUILIBRUM TWO-PHASE REGION
        // ITERATE ON QUALITY USING BRENT'S ROUTINE
        double TOL = 1.0e-08;
        QUAL = ZBRENT(QTWO,QLOW,QHIGH,TOL);
      }
      else{
        // TWO-PAHES BRACKETING FAILED, BASED ON LATEST QUALITY
        // NON-EQUILIBRUM REGION
        if (abs(TESTHIGH) < abs(TESTLOW))
          QUAL = QHIGH;
        else if (abs(TESTHIGH) > abs(TESTLOW))
          QUAL = QLOW;
        else
          QUAL = QUAL;
        
        // USE TTWO AND HGL FROM QTWO
        double TEST = QTWO(QUAL);
        PSAT = pow(10.0,ABOIL-BBOIL/TTWO);
        // GET RHOL FROM ITERATIONS OF NONLINEAR FUNCTIONS F1,F2
        RHOL = s_eosdat.RHO0;
        NEWTP(HGL,PSAT,RHOL);
        // DETERMINE NON-EQUILIBRUM VAPOR DENSITY
        if (QUAL < 0.0)
          RHOV = RHOL / (1.0-(1.0-(RHOL/RHO))/QUAL);
        else
          RHOV = PSAT / (RBAR * TTWO);
      }
      // COMPUTE THE REST OF THE TWO-PHASE PROPERTIES
      T = TTWO;
      // ZONE PRESSURE WILL BE SATURATION PRESSURE
      P = pow(10.0,ABOIL - BBOIL/T);
      VOLF = QUAL * RHO / RHOV;
      if (QUAL > 0.999)
        VOLF = 1.0 < VOLF ? 1.0 : VOLF;
      
      if(VOLF >= 0.0 && VOLF <= 1.0){
        // NOW COMPUTE SOUND SPEDD
        EOSF(RHOL,F1,F2);
		double CAS = 0.0;
        double CLIQ;
        SOUND(HGL,RHOL,P,F2,CLIQ,CAS);
        double CVAP = sqrt(RBAR*GAMMA*T);
        // VOLUME FRACTION OF VAPOR FOR CRUDE INITIAL GUESS AT SOUND SPEED
        CSOUND = VOLF*CVAP + (1.-VOLF)*CLIQ;
      }
      else{
        // THERE IS A PROBLEM WITH VOLF
        cerr << "PROBLEM WITH VOLF IN EOSTWO" << endl;
        cerr << "VOLF = " << VOLF << endl;
        cerr << "QUAL = " << QUAL << endl;
        cerr << "RHO = " << RHO << endl;
        cerr << "RHOV = " << RHOV << endl;
        exit(EXIT_FAILURE);
      }
    }
  }
}

/*
 C EQUATION OF STATE FOR 2-PHASE SPALLED MATERIAL (LIQUID AND VAPOR)
 C ALSO COMPUTE SOUND SPEED
 C FOR LIQUID PART, USE CURRENT LIQUID/SOLID EOS
 C FOR VAPOR PART, USE IDEAL GAS RELATION
 C BASE TWO-PHASE FITTING ON Tsat vs. Psat CURVE
 C UNITS: AL,AV COEFFICIENTS IN J/kg, PRESSURE IN Pa
 */
double QTWO(double QUAL)
{
  //UNSURE
  double& HG = s_qtwodat.HGTWO;
  double& RHO = s_qtwodat.RHOTWO;
  double& RHOV = s_qtwodat.RHOV;
  double& RHOL = s_qtwodat.RHOL;
  double& T = s_qtwodat.TTWO;
  double& HGL = s_qtwodat.HGL;
  
  double& AL0 = s_matdat.AH2T[6];
  double& AL1 = s_matdat.AH2T[7];
  double& AV0 = s_matdat.AH2T[8];
  double& AV1 = s_matdat.AH2T[9];
  double& GAMMA = s_eosdat.EOSCOEF[8];
  double& RBAR = s_eosdat.EOSCOEF[9];
  
  // DETERMINE TEMPERATURE BASED ON ENERGY BALANCE
  double Z1 = HG + QUAL*AV0/AV1 + (1.0-QUAL)*AL0/AL1;
  double Z2 = (QUAL/AV1) + ((1.0-QUAL)/AL1);
  T = Z1 / Z2;
  
  // ENTHALPIES FOLLOW DIRECTLY
  double HGV = (T-AV0) / AV1;
  HGL = (T-AL0) / AL1;
  
  // ZONE PRESSURE WILL BE SATURATION PRESSURE
  double PSAT = pow(10.0,s_eosdat.ABOIL - s_eosdat.BBOIL/T);
  // VAPOR DENSITY FROM IDEAL GAS RELATION
  RHOV = PSAT / (RBAR * T);
  // GET RHOL FROM ITERATIONS ON NONLINEAR FUNCTIONS F1,F2
  RHOL = s_eosdat.RHO0;
  NEWTP(HGL,PSAT,RHOL);
  
  // GET QUALITY IMPLIED BY MASS BALANCE
  double QUALRHO = (1.0-RHOL/RHO) / (1.0-RHOL/RHOV);
  
  // QUALITIES MUST BE EQUAL IN CONVERGED SOLUTION
  return QUAL - QUALRHO;
}

/*
 C NEWTON'S METHOD ITERATION TO FIND DENSITY,
 C GIVEN PRESSURE AND SPECIFIC ENERGY
 C F1 = F1(RHO) AND F2 = F2(RHO)
 C FOR LINEAR EOS: PRESSURE = F1 + F2 * ENERGY
 */
void NEWTP(double& HG, double& P, double& RHO)
{
  double TOL = 1.0e-06;
  int K = 0;
  double RHONEW = RHO;
  //TO CHECK
  do{
    RHO = RHONEW;
    double F1,F2;
    EOSF(RHO,F1,F2);
    double DF1,DF2;
    EOSDF(RHO,DF1,DF2);
    double F = F1 + F2*HG - P;
    double FP = DF1 + HG*DF2;
    RHONEW = RHO - F/FP;
    
    K++;
    if(K < 10){
      cerr << "TOO MANY ITERATIONS ON RHOL" << endl;
      exit(EXIT_FAILURE);
    }
  }while(abs(RHONEW-RHO) > TOL);
  RHO = RHONEW;
}

/*
 C CALCULATION OF SURFACE EVAPORATION MASS AND ENERGY LOSS
 C FOR CRUDE ESTIMATE USING MAXIMUM FLUX (PAMB=0)
 C LIMIT FROM KELLEY & ROTHENBERG, NUC INST & METH, 1985, p 755
 */
void SURFVAP(int& GEOM, double& PAMB, double& TSAT, double& X, double& DT, double& DVAPMASS, double& DHVAP)
{
  double& RBAR = s_eosdat.EOSCOEF[9];
  // VAPOR PRESSURE FUNCTION FROM DUSHMAN, PRESSURE IN Pa
  double PSAT = pow(10.0,s_eosdat.ABOIL - s_eosdat.BBOIL/TSAT);
  // FLUX IN kg/m2.sec
  double FLUX = PSAT / sqrt(2. * PI * RBAR * TSAT);
  double AREA;
  if (GEOM == 1)
    AREA = 1.0;
  else if (GEOM == 2)
    AREA = 2.0 * PI * X;
  else
    AREA = 4.0 * PI * pow(X,2.0);
  
  DVAPMASS = FLUX * DT * AREA;
  
  // ENERGY IN THIS MASS
  // ASSUME VAPOR AT CURRENT ZONE TEMPERATURE
  // USES LINEAR FIT TO T VS. H CURVE (in K vs. J/kg)
  double HG = (TSAT - s_matdat.AH2T[8]) / s_matdat.AH2T[9];
  double HL = (TSAT - s_matdat.AH2T[6]) / s_matdat.AH2T[7];
  DHVAP = (HG - HL) * DVAPMASS;
}

/*
 C COMPUTE EQUATION OF STATE COEFFICIENTS FOR SURFACE ZONE
 C F1 = F1(RHO) AND F2 = F2(RHO)
 C FOR LINEAR EOS: PRESSURE = F1 + F2 * ENERGY
 C ASSUMES IDEAL GAS EOS: P = (GAMMA-1)*RHO*HG
 C OR JUST USE SOLID EOS FOR LOW VAPOR FRACTION
 */
void FSURF(double& RHO, double& RHOL, double& QUAL, double& QMIN, double& F1, double& F2)
{
  if(QUAL <= QMIN){
    // JUST LIQUID/SOLID IN ZONE, SO IGNORE VAPOR
    EOSF(RHO,F1,F2);
  }
  else{
    // APPROXIMATION BASED ON CONSTANT RHOL FROM LAST TIME STEP
    double RHOV = QUAL / ((1./RHO) + (QUAL - 1.0)/RHOL);
    FVAPOR(RHOV,F1,F2);
  }
}

/*
 C ROUTINE TO DETERMINE NEW VELOCITY OF OUTER NODE
 C BASED ON PSEUDO-STEADY ISENTROPIC IDEAL GAS FLOW REALTIONS
 C ASSUME VAPORIZING SURFACE IS STATIONARY
 C VAPORIZING SURFACE AT P0=Psat, OUTER SURFACE AT P=P(0)
 */
void USURFIX(double& T, double& P0, double& P, double& UNEW)
{
  double& GAMMA = s_eosdat.EOSCOEF[8];
  double& RBAR = s_eosdat.EOSCOEF[9];
  double GM1 = GAMMA - 1.0;
  
  double C0 = sqrt(RBAR*GAMMA*T);
  
  UNEW = C0*sqrt((2.0/GM1)*abs(1.0-pow(P/P0,GM1/GAMMA)));
  
  // KEEP WITH SIGN CONVENTION
  if (P < P0)
    UNEW = -UNEW;
}

/*
 C ROUTINE TO DETERMINE NEW PRESSURE AT OUTER NODE OF ISURF
 C BASED ON PSEUDO-STEADY ISENTROPIC IDEAL GAS FLOW REALTIONS
 C C/C0 = 1 - B*X LINEAR ACCROS EXPANSION FAN
 C BUT ASYMPTOTE TO IDEAL GASS PRESSURE IN SPALLED ZONE
 */
void PLEFTCALC(double& PSAT,double& PC,double& DXC1,double& DX1,double& RHOV,double& QUAL,double& T,double& PLEFT)
{
  double& GAMMA = s_eosdat.EOSCOEF[8];
  double& RBAR = s_eosdat.EOSCOEF[9];
  double GM1 = GAMMA - 1.0;
  //prejmenovani kvuli konfliktu s define N
  double N0 = 4.0;
  
  // ASSUME AVERAGE DROP ACROSS SONIC KNUDSEN LAYER (KNIGHT)
  // P0 = 0.235 * PSAT
  // ASSUME NO DROP ACROSS KNUDSEN LAYER
  double P0 = PSAT;
  double B = (1.0-pow(PC/P0,GM1/(2.*GAMMA)))/DXC1;
  double PEXP = P0 * pow(1.0 - B*DX1,(2.0*GAMMA/GM1));
  PEXP = PEXP > PC ? PEXP : PC;
  
  // AVERAGE VAPOR PRESSURE
  double PVAP = RHOV * RBAR * T;
  
  // INTERPOLATE BETWEEN THESE PRESSURES WITH Nth ORDER FUNCTION
  double QSTAR = pow(1.0 - QUAL,N0);
  PLEFT = PEXP * QSTAR + (1.0 - QSTAR) * PVAP;
}

/*
 C EQUATION OF STATE FOR (IN GENERAL) 2-PHASE SURFACE ZONE
 C FOR LIQUID PART, USE CURRENT LIQUID/SOLID EOS
 C FOR VAPOR PART, USE IDEAL GAS RELATION
 C UNITS: AL,AV COEFFICIENTS IN J/kg, PRESSURE IN Pa
 C NOTE THAT T PROVIDES BOTH AN INITIAL TEMPERATURE GUESS
 C  AND RETURNS THE FINAL TEMPERATURE RESULT
 */
void EOSURF(double& RHOL,double& RHOV,double& RHOP1,double& QMIN,double& PM1,double& P,double& ALPHA,double& CSOUND)
{
  double (*TSURF)(double);
  TSURF = &TSURF_DEF;

  double (*PSURF)(double);
  PSURF = &PSURF_DEF;

  double& AL0 = s_matdat.AH2T[6];
  double& AL1 = s_matdat.AH2T[7];
  double& AV0 = s_matdat.AH2T[8];
  double& AV1 = s_matdat.AH2T[9];
  double& GAMMA = s_eosdat.EOSCOEF[8];
  double& RBAR = s_eosdat.EOSCOEF[9];
  double& QUAL = s_surft.QUALSURFT;
  double& HG = s_surft.HGSURF;
  double& RHO = s_surfp.RHOSURF;
  double& T = s_surfp.TEMPSURF;
  double& HGV = s_surfp.HGVSURF;
  double& HGL = s_surfp.HGLSURF;
  
  // ITERATE ON TEMPERATURE CHOICE TO MATCH VALUES OBTAINED FOR LIQUID AND VAPOR
  // SOLVE ISING BRENT'S METHOD (SEE NUMERICAL RECIPES)
  // HG1 AND HG2 MUST BRACKET THE ROOT TO Tvapor-Tliquid=0,
  // WHICH IS IN THE EXTERNAL FUNCTION TSURF(HGV)
  // WHERE HGV IS THE SPECIFIC ENERGY OF THE VAPOR PHASE
  // FIRST GET TWO HGV VALUES THAT BRACKET FUNCTION ZERO
  double HGNOM = (T - AV0) / AV1;

  double HG1,HG2;
  BRACKET(TSURF,HGNOM,HG1,HG2);
  double TOL = 1.0e-05;
  
  HGV = ZBRENT(TSURF,HG1,HG2,TOL);
  HGL = (HG - QUAL*HGV) / (1.0 - QUAL);
  
  double TV = AV0 + AV1 * HGV;
  double TL;
  H2T(HGL,TL);
  
  double TRATIO = 2.0*abs(TV-TL)/(TV+TL);
  if (TRATIO > 0.01){
    cerr << "PROBLEM IN ISURF TEMP MATCH" << endl;
    cerr << "TL, TV = " << TL << " , " << TV << endl;
    exit(EXIT_FAILURE);
  }
  else
    T = (TL + TV) / 2.00;
  
  // NOW DETERMINE PRESSURE IN THE SURFACE ZONE
  // FOR SMALL VAPOR CONTENT, ASSUME FULLY LIQUID/SOLID AS IN FSURF
  // OTHERWISE, FIX PRESSURE AT Psat
  double F1,F2;
  if(QUAL <= QMIN){
    RHOL = RHO;    
    EOSF(RHO,F1,F2);
    P = F1 + F2 * HG;
	double CAS = 0.0;
    SOUND(HG,RHO,P,F2,CSOUND,CAS);
	
    ALPHA = 0.0;
    RHOV = P / (RBAR * T);
  }
  else{
    // SET PRESSURE EQUAL TO THE SATURATION VAPOR PRESSURE
    // OR THE PRESSURE IN THE NEXT VAPOR ZONE, WHICHEVER IS HIGHER
    double PSAT = pow(10.0,s_eosdat.ABOIL - s_eosdat.BBOIL/T);
    P = PSAT > PM1 ? PSAT :PM1;
    // UPDATE LIQUID DENSITY BASED ON THIS PRESSURE
    // MUST ITERATE WITH BRENT'S METHOD
    // FIRST GET TWO RHOL VALUES THAT BRACKET FUNCTION ZERO
    double RNOM = RHOL;
    double RL1,RL2;
	
    BRACKET(PSURF,RNOM,RL1,RL2);
    double TOL = 1.0e-05;
    RHOL = ZBRENT(PSURF,RL1,RL2,TOL);
    // UNDER SOME CONDITIONS
    // RHOL MAY CALCULATE TO BE LESS THAN RHO
    // SET A FLOOR DENSITY ON WHICH TO BASE A RHOV CALCULATION
    // ADJUST ZONE PRESSURE ACCORDINGLY
    if (RHOL < RHO){
      if (RHOP1 > RHO)
        RHOL = RHOP1;
      else
        RHOL = RHO + 1.0;
      
      EOSF(RHOL,F1,F2);
      P = F1 + F2 * HGL;
    }
    
    // COMPUTE VOLUME FRACTION VAPOR
    RHOV = QUAL / ((1.0/RHO) + (QUAL - 1.0)/RHOL);
    ALPHA = QUAL * RHO / RHOV;
    // NOW COMPUTE SOUND SPEED
    EOSF(RHOL,F1,F2);
	double CAS=0.0;
    double CLIQ;
    SOUND(HGL,RHOL,P,F2,CLIQ,CAS);
    double CVAP = sqrt(RBAR*GAMMA*T);
    // VOLUME FRACTION OF VAPOR FOR CRUDE GUES AT SOUND SPEED
    CSOUND = ALPHA*CVAP + (1.-ALPHA)*CLIQ;
  }
}

/*
 C TEMPERATURE ITERATION FOR 2-PHASE SURFACE ZONE
 C NOTE THAT T PROVIDES BOTH AN INITIAL TEMPERATURE GUESS
 C  AND RETURNS THE FINAL TEMPERATURE RESULT
 */
void EOSURFT(double& QMIN,double& T)
{
  double (*TSURF)(double);
  TSURF = &TSURF_DEF;

  double& AL0 = s_matdat.AH2T[6];
  double& AL1 = s_matdat.AH2T[7];
  double& AV0 = s_matdat.AH2T[8];
  double& AV1 = s_matdat.AH2T[9];
  double& GAMMA = s_eosdat.EOSCOEF[8];
  double& RBAR = s_eosdat.EOSCOEF[9];
  double& QUAL = s_surft.QUALSURFT;
  double& HG = s_surft.HGSURF; 
  
  if (QUAL < QMIN){
    // ALL LIQUID OR SOLID
    
    H2T(HG,T);
  }
  else if (QUAL > 0.99999){
    // ESSENTIALLY ALL VAPOR
    T = AV0 + AV1 * HG;
  }
  else{
    // ITERATE ON TEMPERATURE CHOICE TO MATCH VALUES OBTAINED FOR LIQUID AND VAPOR
    // SOLVE ISING BRENT'S METHOD (SEE NUMERICAL RECIPES)
    // HG1 AND HG2 MUST BRACKET THE ROOT TO Tvapor-Tliquid=0,
    // WHICH IS IN THE EXTERNAL FUNCTION TSURF(HGV)
    // WHERE HGV IS THE SPECIFIC ENERGY OF THE VAPOR PHASE
    // FIRST GET TWO HGV VALUES THAT BRACKET FUNCTION ZERO
    double HGNOM = (T - AV0) / AV1;
    double HG1,HG2;
    BRACKET(TSURF,HGNOM,HG1,HG2);
    double TOL = 1.0e-05;
    double HGV = ZBRENT(TSURF,HG1,HG2,TOL);
    double HGL = (HG - QUAL*HGV) / (1.0 - QUAL);
    double TV = AV0 + AV1 * HGV;
    double TL;
    H2T(HGL,TL);
    double TRATIO = 2.0*abs(TV-TL)/(TV+TL);
    if(TRATIO > 0.01){
      cerr << "PROBLEM IN ISURF TEMP MATCH" << endl;
      cerr << "TL, TV = " << TL << ", " << TV << endl;
	  //return;
      exit(EXIT_FAILURE);	  
    }
    else{
      T = (TL + TV) / 2.00;
    }
  }
}

/*
 C TEMP MATCHING FUNCTION FOR (IN GENERAL) 2-PHASE SURFACE ZONE
 C FORM IS TSURF = Tvap - Tliq -> 0 WITH BRENT'S METHOD
 C TSURF IS A FUNCTION OF VAPOR ENERGY (HGV)
 C ROOT SOLVER WILL VARY ENERGY PARTITION VAPOR/LIQUID TO MATCH TEMPS
 */

double TSURF_DEF(double HGV)
{
  double& AV0 = s_matdat.AH2T[8];
  double& AV1 = s_matdat.AH2T[9];
  
  double& QUAL = s_surft.QUALSURFT;
  double& HG = s_surft.HGSURF;
  
  double HGL = (HG - QUAL*HGV) / (1.0 - QUAL);
  double TV = AV0 + AV1 * HGV;
  double TL;
  H2T(HGL,TL);
  
  return TV - TL;
}

/*
 C PRESSURE MATCHING FUNCTION FOR 2-PHASE SURFACE ZONE
 C FORM IS PSURF = Psat - Tliq -> 0 WITH BRENT'S METHOD
 C PSURF IS A FUNCTION OF LIQUID DENSITY (RHOL)
 C ROOT SOLVER WILL VARY LIQUID DENSITY TO MATCH PRESSURES
 */
double PSURF_DEF(double RHOL)
{
  double& HGL = s_surfp.HGLSURF;
  double& ABOIL = s_eosdat.ABOIL;
  double& BBOIL = s_eosdat.BBOIL;
  double& T = s_surfp.TEMPSURF;
  
  double PSAT = pow(10.0,ABOIL - BBOIL/T);
  
  double F1L,F2L;
  EOSF(RHOL,F1L,F2L);
  
  double PL = F1L + F2L * HGL;
  
  return PSAT - PL;
}

/*
 C BRENT'S METHOD TO FIND THE ROOT OF A FUNCTION FUNC KNOWN TO
 C LIE BETWEEN X1 AND X2. THE RROT, ZBRENT, WILL BE REFINED UNTIL
 C ITS ACCURACY IS TOL
 C PARAMETERS: MAX # OF ITERATIONS, MACHINE FLOATING-POINT PRECISION
 C REFERENCE: NUMERICAL RECIPES, SECTION 9.3
 C IS RECURSIVE
 */
double ZBRENT(double (*FUNC)(double),double& X1,double& X2,double& TOL)
{
  int ITMAX = 100;
  double EPS = 1.e-15;
  
  //unreferenced
  //int ITER;
  
  double A,B,C,D,E,FA,FB,FC,P,Q,R,S,TOL1,XM;
  
  A = X1;
  B = X2;
  FA = (*FUNC)(A);
  FB = (*FUNC)(B);
  if ((FA > 0.0 && FB > 0.0) || (FA < 0. && FB < 0.0)){
    cout << "ROOT MUST BE BRACKETED FOR ZBRENT" << endl;
    cout << "X1,X2 =" << X1 << "," << X2 << endl;
    cout << "F(X1),F(X2) =" << endl << FA << "," << FB << endl;
  }
  C = B;
  FC = FB;
  for(int i=1;i<=ITMAX;i++){
    if((FB > 0.0 && FC > 0.0) || (FB < 0.0 && FC < 0.0)){
      C = A;
      FC = FA;
      D = B - A;
      E = D;
    }
    if (abs(FC) < abs(FB)){
      A = B;
      B = C;
      C = A;
      FA = FB;
      FB = FC;
      FC = FA;
    }
    // CONVERGENCE CHECK
    TOL1 = 2. * EPS * abs(B) + 0.5 * TOL;
    XM = 0.5 * (C - B);
    if (abs(XM) <= TOL1 || FB == 0.0){
      return B;
    }
    
    if(abs(E) >= TOL1 && abs(FA) > abs(FB)){
      // ATTEMPT INVERSE QUADRATIC INTERPOLATION
      S = FB / FA;
      if(A == C){
        P = 2. * XM * S;
        Q = 1.0 - S;
      }
      else{
        Q = FA / FC;
        R = FB / FC;
        P = S * (2.*XM*Q*(Q-R) - (B-A)*(R-1.0));
        Q = (Q-1.0) * (R-1.0) * (S-1.0);
      }
      
      if(P > 0.0)
        Q = -Q;
      // CHECK WHETHER IN BOUNDS
      P = abs(P);
	  
      if(2.0*P < min(3.0*XM*Q-abs(TOL1*Q), abs(E*Q))){
        // ACCEPT INTERPOLATION
        E = D;
        D = P / Q;
      }
      else{
        // INTERPOLATION FAILED, USE BISECTION
        D = XM;
        E = D;
      }
    }
    else{
      // BOUNDS DECREASING TOO SLOWLY, USE BISECTION
      D = XM;
      E = D;
    }
    
    // MOVE LAST BEST GUESS TO A
    A = B;
    FA = FB;
    // EVALUEATE NEW TRIAL ROOT
    if(abs(D) > TOL1)
      B = B + D;
    else{
		double s;
		if(XM >= 0){
			s = 1;
		}
		else{
			s = -1;
		}
		B = B + abs(TOL1)*s;

	}
    
    FB = (*FUNC)(B);
  }
  
  cout << "ZBRENT EXCEEDED MAXIMUM ITERATIONS" << endl;
  return B;
}

/*
 C ROUTINE TO BRACKET ROOT OF FUNCTIONS, TO GIVE INPUT TO
 C BRENT'S METHOD TO FIND THE ROOT OF A FUNCTION FUNC
 C SINGLE-SIDED SEARCH, WITH XMIN AS LOWER LIMIT ALLOWED
 */

/*
 * this function is not used in the former fortran code and will not be translated yet
 */

/*
 SUBROUTINE BRKTMIN(FUNC,XMIN,XNOM,X1,X2)
 DOUBLE PRECISION XMIN,XNOM,X1,X2,FUNC
 EXTERNAL FUNC
 DOUBLE PRECISION A,B,FA,FB
 
 A = MAX(XMIN,XNOM*.85)
 B = A * 1.02
 FA = FUNC(A)
 
 DO 10, ICNT=1,20
 FB = FUNC(B)
 IF ((FA .GT. 0. .AND. FB .LT. 0.) .OR.
 &      (FA .LT. 0. .AND. FB .GT. 0.)) THEN
 X1 = A
 X2 = B
 RETURN
 END IF
 B = B * 1.02
 10  CONTINUE
 WRITE (*,*) 'TOO MANY BRACKET ITERATIONS WITH XMIN'
 B = A * 1.02
 FB = FUNC(B)
 WRITE (*,*) 'A,B=',A,B
 WRITE (*,*) 'FA,FB=',FA,FB
 STOP
 END
 */


/*
 C ROUTINE TO BRACKET ROOT OF FUNCTIONS, TO GIVE INPUT TO
 C BRENT'S METHOD TO FIND THE ROOT OF A FUNCTION FUNC
 C BASED ON USING A NOMINAL X AS STARTING POINT
 C IS RECURSIVE
 */
void BRACKET(double (*FUNC)(double),double& XNOM,double& X1,double& X2)
{
  double A,B,C,FA,FB,FC;
  int IFLAG = 0;
  
  if(XNOM == 0.0){
    A = -1.0;
    B = 1.0;
  }
  else{
    A = XNOM * 0.9999;
    B = XNOM * 1.02;
  }
  FA = (*FUNC)(A);
  FB = (*FUNC)(B);
  
  do{
    if((FA > 0.0 && FB < 0.0) || (FA < 0.0 && FB > 0.0)){
      X1 = B;
      X2 = A;
      return;
    }
    // PICK A NEW POINT WITH LINEAR INTERPOLATION + 2%
    for(int ICNT=1;ICNT<=20;ICNT++){
      C = B - 1.02 * (A-B) * FB / (FA-FB);
      FC = (*FUNC)(C);
      if((FA > 0.0 && FC < 0.0) || (FA < 0.0 && FC > 0.0)){
        if(abs(A-C) < abs(B-C)){
          X1 = C;
          X2 = A;
          return;
        }
        else{
          X1 = C;
          X2 = B;
          return;
        }
      }
      // THROW AWAY POINT FARTHEST FROM NEW POINT, AND TRY AGAIN
      if(abs(A-C) < abs(B-C)){
        B = C;
        FB = FC;
      }
      else{
        A = C;
        FA = FC;
      }
      
      if(IFLAG == 1){
        cout << A << "," << FA << endl;
        cout << B << "," << FB << endl;
      }
    }
    
    // TOO MANY BRACKET ITERATIONS
    A = XNOM * 10.000;
    B = XNOM * (-10.00);
    FA = (*FUNC)(A);
    FB = (*FUNC)(B);
    IFLAG = IFLAG+1;
  }while(IFLAG == 1);
  cout << "TOO MANY BRACKET ITERATIONS" << endl;
  cout << "A,B=" << A << "," << B << endl;
  cout << "FA,FB=" << FA << "," << FB << endl;
  exit(EXIT_FAILURE);
}


/*
 C FUNCTION FOR IMPLICIT HEAT CONDUCTION BETWEEN ZONES
 C ISURF AND (ISURF+1), TO BE SOLVED WITH BRENT'S METHOD
 C RETURNS DIFFERENCE BETWEEN HEAT CONDUCTION (INPUT)
 C AND FOURIER CONDUCTION BASED ON TEMPS AT END OF TIME STEP
 C HPN & HEN ARE ENERGIES AT STEP N, REST OF H & T AT N+1
 */
double CONDSURF(double HEAST)
{
  
  double& HPN = s_surfc.SURFCSTUFF[0];
  double& HEN = s_surfc.SURFCSTUFF[1];
  double& ZMP = s_surfc.SURFCSTUFF[2];
  double& ZME = s_surfc.SURFCSTUFF[3];
  double& DT = s_surfc.SURFCSTUFF[4];
  double& QMIN = s_surfc.SURFCSTUFF[5];
  double& CONDP = s_surfc.SURFCSTUFF[6];
  double& DX = s_surfc.SURFCSTUFF[7];
  double& BETAC = s_surfc.SURFCSTUFF[8];
  double& TP = s_surfc.SURFCSTUFF[9];
  
  double& HGP = s_surft.HGSURF;
  double& QUAL = s_surft.QUALSURFT;
  
  // COMPUTE NEW ZONE ENERGIES BASED ON CURRENT HEAST
  double HP = HPN - HEAST * DT;
  double HE = HEN + HEAST * DT;
  HGP = HP / ZMP;
  double HGE = HE / ZME;
  // UPDATE TEMPERATURES
  EOSURFT(QMIN,TP);  
  double TE;
  H2T(HGE,TE);
  // COMPUTE FOURIER CONDUCTION BASED ON NEW TEMPS
  double HEASTNEW = (TP - TE) * CONDP * BETAC / DX;
  
  return HEASTNEW - HEAST;
}

// **********************************************************
// CALCULATE SURFACE TENSION FROM TEMPERATURE
// USE CUBIC OR EXPONENTIAL FIT TO DATA
double SIGMACALC(double T)
{
	double* EOSCOEF = s_eosdat.EOSCOEF;
	double* SIGMACOEF = s_surften.SIGMACOEF;
	double SIGMACALC;

	int JFLAG = int(SIGMACOEF[0]);
	if(JFLAG == 1){
// FIT PARAMETERS FOR FIT TO MURR (CUBIC EQUATION)
		double A = SIGMACOEF[1];
		double TCRIT = SIGMACOEF[2];
		double SIGMA0 = SIGMACOEF[3];
		double ASIG = SIGMACOEF[4];
		double BSIG = SIGMACOEF[5];
		double CSIG = SIGMACOEF[6];
		double T1 = (T - s_matdat.TMELT);
		SIGMACALC = SIGMA0 + T1*(ASIG+T1*(BSIG+T1*CSIG));
	}
	else if (JFLAG == 2){
// FIT PARAMETERS FOR EXPONENTIAL EQUATION
		double A = SIGMACOEF[1];
		double TCRIT = SIGMACOEF[2];
		double AEXP = SIGMACOEF[3];
		double BEXP = SIGMACOEF[4];
		SIGMACALC = AEXP * pow(10.0,(BEXP * T));
		}
	else{
		cout << "INVALID MATERIAL (SURF TENSION)" << endl;
		exit(EXIT_FAILURE);
	}

	return SIGMACALC;
}


/*
 C BUBBLE NUCLEATION RATE ROUTINE
 C FROM CAREY, EQN 5.74
 C MATERIAL PROPERTIES HARDWIRED FOR AL,ALUMINA,SILICA
 C UNITS: J in m-3.s-1, SIGMA in N/m
 */
void NUCLEATE(double& T,double& RHO,double& P,double& J,double& TSTAR)
{
  double ISTAR;
  
  double* EOSCOEF = s_eosdat.EOSCOEF;
  double* SIGMACOEF = s_surften.SIGMACOEF;
  double A = SIGMACOEF[1];
  double TCRIT = SIGMACOEF[2];
  double SIGMA0 = SIGMACOEF[3];
  double ASIG = SIGMACOEF[4];
  double BSIG = SIGMACOEF[5];
  double CSIG = SIGMACOEF[6];
  
  double& RBAR = s_eosdat.EOSCOEF[9];
  int JFLAG = int(SIGMACOEF[0]);
  
  if (JFLAG == 0){
    // NO MATERIAL SURFACE TENSION INFORMATION -> SKIP J CALCULATION
    J = 2.0e-97;
    TSTAR = 0.1;
  }
  else{
    double SIGMA0,ASIG,BSIG,CSIG;
	
    double AVOGADRO = 6.022e+26;
    double GAMMAE = 0.5772;
    double BOLTZ = 1.38e-23;
    
    if(T < s_matdat.TMELT){
      cout << "T OUTSIDE RANGE FOR NUCLEATE" << endl;
      J = 1.0e-98;
      TSTAR = 0.1;
    }
    //else if (T >= TCRIT){	
	else if (T >= 0.0) {
      if(P > 0.0){
        J = 1.0e-96;
        TSTAR = 0.1;
      }
      else{
        J = 1.0e+43;
        TSTAR = 1.0e-19;
      }
    }
	else {
		// COMPUTE SATURATION PRESSURE
		double PSAT = pow(10.0, s_eosdat.ABOIL - s_eosdat.BBOIL / T);

		//WARNING: zajimalo by me co se stane kdyz nebude SIGMA nastavena
		double SIGMA = -1;
		//WARNING
		double ETA = -1;

		if (P > PSAT) {
			J = 1.0e-96;
			TSTAR = 0.1;
		}
		else {
			// COMPUTE SURFACE TENSION      
			SIGMA = SIGMACALC(T);

			double ETA = exp((P - PSAT) / (RHO*RBAR*T));
			double X1 = -1.8195e+24 * pow(SIGMA, 3.0) / (T*pow(ETA*PSAT - P, 2.0));
			double X2 = sqrt(pow(RHO, 2.0) * SIGMA / A);

			J = 1.44e+40 * X2 * exp(X1);
			J = J > 1.0e-97 ? J : 1.e-97;
		}

		// COMPUTE INDUCTION TIME REQUIRED TO REACH STEADY J-VALUE
		// FROM MODIFICATION OF VAPOR CONDENSATION MEO BY WILEMSKI
		// RSTAR IS THE CRITICAL BUBBLE RADIUS (meters)
		if (SIGMA == -1.0) {
			cout << "ERR: NUCLEATE - value of sigma not set" << endl;
		}

		double RSTAR = 2.9 * SIGMA / (ETA*PSAT - P);
		if (RSTAR <= 0.0)
			TSTAR = 1.0e-02;
		else {
			double U = sqrt(8.*RBAR*T / PI);
			double RHOV = PSAT / (RBAR * T);
			double SMALLM = A / AVOGADRO;
			double SMALLN = ETA * PSAT / (BOLTZ * T);
			ISTAR = (4. / 3.) * PI * pow(RSTAR, 3.0) * RHOV / SMALLM;
			double WSTAR = (16.0*PI / 2.)*pow(SIGMA, 3.0) / (BOLTZ*T*pow((ETA*PSAT - P), 2.0));
			TSTAR = RSTAR*ISTAR*RHOV / (U*SMALLM*SMALLN*WSTAR) * (log(WSTAR / 3.0) + GAMMAE);
			if (TSTAR <= 0.0) {
				TSTAR = 1.0e-18;
			}
		}
	}
  }
}
  
  //TODO
  //101 FORMAT(7(1PE12.3))
  
  /*
   * Ta sekce pokracuje za returnem, nechapu smysl
   *
   RETURN
   ENDR = 2. * SIGMA / (ETA*PSAT - P)
   IF (RSTAR .LE. 0.) THEN
   TSTAR = 1.0D-02
   ELSE
   U = SQRT(8.*RBAR*T/PI)
   RHOV = PSAT / (RBAR * T)
   SMALLM = A / AVOGADRO
   SMALLN = ETA * PSAT / (BOLTZ * T)
   ISTAR = (4./3.) * PI * RSTAR**3 * RHOV / SMALLM
   WSTAR = (16.*PI/2.)*SIGMA**3/(BOLTZ*T*(ETA*PSAT-P)**2)
   TSTAR = RSTAR*ISTAR*RHOV / (U*SMALLM*SMALLN*WSTAR) * 
   &          (DLOG(WSTAR/3.) + GAMMAE)
   IF (TSTAR .LE. 0.) TSTAR = 1.0D-18
   END IF
   END 
   */
