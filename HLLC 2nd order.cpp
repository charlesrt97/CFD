// Riemann solver implemented: HLLC 2nd order, along with the Godunov scheme
// Slope limiter used: minmod

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <math.h>
#include <time.h>
using namespace std;

/******************************************************************************/

// constant parameters used for the simulation
const int    NX = 400;         // mesh size
const double XL = 0.0;         // left phyisical coordinate
const double XR = 1.0;         // right phyisical coordinate
const double TFIN = 0.2;       // time of integration
const double CFL = 0.9;        // courant time 
const double dtprint = 0.005;   // time interval to write to disk

const double gam=1.4;

const int ieq=3;

// for the schock-tube: coordinate of separation between initial states
const double X0 = 0.5;

// for advection: wave's velocity
const double A = 1.0;

const double DX = (XR-XL)/NX;      // space between nodes

// global variables
double U[ieq][NX+4];       // "current" conservative variables
double UP[ieq][NX+4];      // "advanced" conservative variables
double F[ieq][NX+4];       // physical fluxes
double P[ieq][NX+4];
double Fhllc1[ieq][NX+4];
const double eta = 0;
//double prim[ieq][NX+4];

double pL[ieq];
double pR[ieq];

double fR[ieq];
double fL[ieq];

double UUL[ieq];
double UUR[ieq];

//double csr[ieq][NX+4];
//double csl[ieq][NX+4];

double dt;            // time step
double t;          // current time 
int it;               // current iteration
clock_t start;        // initial time
double tprint;        // time for the following output
int itprint;          // output number


/******************************************************************************/

// sets initial conditions
void initflow(double U[ieq][NX+4]) {

// initialize U in all the domain
// including ghost cells
double x, rho, u, pres;
const double rhol=1.0;
const double ul=0.0;
const double pres_l=1.0;

const double rhor=0.125;
const double ur=0.0;
const double pres_r=0.1;

  for (int i=0; i <= NX+3; i++){
    x=XL+i*DX;
    if (x<=X0){
      rho = rhol;
      u = ul;
      pres = pres_l;
    } else {
      rho = rhor;
      u = ur;
      pres = pres_r;
    }
    U[0][i] = rho;
    U[1][i] = rho*u;
    U[2][i] = 0.5*rho*(u*u)+pres/(gam-1);
  }


  // Initialize other variables
  t = 0;
  it = 0;
  itprint = 0;
  tprint = 0;

}

/******************************************************************************/

// writes to disk the state of the simulation
void output(double P[ieq][NX+4]) {

  // generates output file's name
  char fname[80];
  sprintf(fname, "Lax_%02i.txt", itprint);

  // opens the file
  fstream fout(fname, ios::out);

  // writes U values to disk
  double x;
  for (int i=2; i <= NX+1; i++) {
    x = XL + i*DX;
    fout << x << " " << P[0][i] << " " << P[1][i] << " " << P[2][i] << " " << endl;
  }

  // closes the file
  fout.close();

  printf("Se escribió salida %s\n", fname);

  itprint = itprint + 1;
  tprint = itprint * dtprint;

}

/******************************************************************************/

// applies boundary conditions to ghost cells
void boundary(double U[ieq][NX+4]) {

  for(int iieq=0; iieq<=ieq-1;iieq++){
    U[iieq][0]=U[iieq][2]; // 1st ghost cell
    U[iieq][1]=U[iieq][2]; // 2nd ghost cell

    U[iieq][NX+2]=U[iieq][NX+1]; 
    U[iieq][NX+3]=U[iieq][NX+1];

  }
}

/******************************************************************************/

// computes primitives, including ghost cells
void primitivas1(double U[ieq][NX+4], double P[ieq][NX+4]) {

  for (int i=0; i<=NX+3;i++) {
    //F[i]=A*U[i];
    P[0][i]=U[0][i];
    P[1][i]=U[1][i]/U[0][i];
    P[2][i]=(gam-1)*(U[2][i]-pow(U[1][i],2)/(2*U[0][i]));

  }

}

// computes physical fluxes, including ghost cells
void fluxes(double P[ieq][NX+4], double F[ieq][NX+4]) {

  for (int i=0; i<=NX+3;i++) {

    F[0][i]=P[0][i]*P[1][i];
    F[1][i]=P[0][i]*pow(P[1][i],2)+P[2][i];
    F[2][i]=P[1][i]*(0.5*P[0][i]*pow(P[1][i],2)+P[2][i]*gam/(gam-1));
  }

}

/******************************************************************************/

// computes new time step resulting from the CFL condition
double timestep(double P[ieq][NX+4]) {

  double dt;

  // for advection eq., max_u is simply A
  // double max_speed = abs(A);

  // for other cases, we shall compute the maximum value abs(vel)
  double max_speed = 0.0;
  double cs;
  double k;
  for (int i = 0; i <= NX+3; i++) {
    cs = sqrt(gam*P[2][i]/P[0][i]);
    k = abs(P[1][i])+cs;
    if (k > max_speed) max_speed = k;
  }

  dt = CFL * DX / max_speed;

  return dt;

}

/******************************************************************************/

void godunov1(double U[ieq][NX+4], double Fhllc1[ieq][NX+4], double UP[ieq][NX+4]) {

// from 1 to NX (physical fluxes)

  for (int i=2; i<=NX+1; i++){
    for (int iieq=0; iieq<=ieq-1; iieq++){

      UP[iieq][i]=U[iieq][i]-dt/(2*DX)*(Fhllc1[iieq][i]-Fhllc1[iieq][i-1]);

    }
  }
}

/******************************************************************************/

void godunov2(double U[ieq][NX+4], double Fhllc1[ieq][NX+4], double UP[ieq][NX+4]) {

// from 1 to NX (physical fluxes)

  for (int i=2; i<=NX+1; i++){
    for (int iieq=0; iieq<=ieq-1; iieq++){

      UP[iieq][i]=U[iieq][i]-dt/(DX)*(Fhllc1[iieq][i]-Fhllc1[iieq][i-1]);

    }
  }
}

/******************************************************************************/

// one time step
void step(double U[ieq][NX+4], double UP[ieq][NX+4]) {

  for (int i = 0; i <= NX+3; i++) {
    for (int iieq=0; iieq<=ieq-1; iieq++){

      U[iieq][i]=UP[iieq][i];

    }
  }

  t = t + dt;
  it = it + 1;

}

// time step including artificial viscosity (for 2nd order schemes, i.e. diffusive schemes)
void stepviscoso(double U[ieq][NX+4], double UP[ieq][NX+4]) {

  for (int i = 2; i <= NX+1; i++) {
    for (int iieq=0; iieq<=2; iieq++){
      //U[iieq][i]=UP[iieq][i];
      if (((UP[iieq][i+1]-UP[iieq][i])*(UP[iieq][i]-UP[iieq][i-1]))<0) {
        U[iieq][i]=UP[iieq][i]+eta*(UP[iieq][i+1]+UP[iieq][i-1]-2.0*UP[iieq][i]);
      }
      else {
        U[iieq][i]=UP[iieq][i];
      }

    }
  }

  t = t + dt;
  it = it + 1;

}

/******************************************************************************/

// applies HLLC scheme
void HLLC(double P[ieq][NX+4], double U[ieq][NX+4], double F[ieq][NX+4], double Fhllc1[ieq][NX+4]) {

// from 0 to NX

double sl, sr, s_star, U_star_l, U_star_r;
double csr, csl;


  for (int i=1; i<=NX+2; i++){
    for (int iieq=0; iieq<=ieq-1; iieq++){

      csl = sqrt(gam*P[2][i]/P[0][i]);
      csr = sqrt(gam*P[2][i+1]/P[0][i+1]);

      sl = min(P[1][i]-csl,P[1][i+1]-csr);
      sr = max(P[1][i]+csl,P[1][i+1]+csr);

      s_star = (P[2][i+1]-P[2][i]+P[0][i]*P[1][i]*(sl-P[1][i])-P[0][i+1]*P[1][i+1]*(sr-P[1][i+1]))/(P[0][i]*(sl-P[1][i])-P[0][i+1]*(sr-P[1][i+1]));

      if (iieq == 0){
        U_star_l=P[0][i]*(sl-P[1][i])/(sl-s_star)*1;
        U_star_r=P[0][i+1]*(sr-P[1][i+1])/(sr-s_star)*1;
      }
      else if (iieq == 1) {
        U_star_l=P[0][i]*(sl-P[1][i])/(sl-s_star)*s_star;
        U_star_r=P[0][i+1]*(sr-P[1][i+1])/(sr-s_star)*s_star;
      }
      else if (iieq == 2) {
        U_star_l=P[0][i]*(sl-P[1][i])/(sl-s_star)*((0.5*P[0][i]*P[1][i]*P[1][i]+P[2][i]/(gam-1))/P[0][i]+(s_star-P[1][i])*(s_star+P[2][i]/(P[0][i]*(sl-P[1][i]))));
        U_star_r=P[0][i+1]*(sr-P[1][i+1])/(sr-s_star)*((0.5*P[0][i+1]*P[1][i+1]*P[1][i+1]+P[2][i+1]/(gam-1))/P[0][i+1]+(s_star-P[1][i+1])*(s_star+P[2][i+1]/(P[0][i+1]*(sr-P[1][i+1]))));
      }

      // numerical intercell fluxes
      if (sl >= 0) {
        Fhllc1[iieq][i] = F[iieq][i];
      }
      else if ((sl <= 0) && (0 <= s_star)) {
        Fhllc1[iieq][i] = F[iieq][i]+sl*(U_star_l-U[iieq][i]);
      }
      else if ((s_star <= 0) && (0 <= sr)) {
        Fhllc1[iieq][i] = F[iieq][i+1]+sr*(U_star_r-U[iieq][i+1]);
      }
      else if (sr <= 0) {
        Fhllc1[iieq][i] = F[iieq][i+1];
      }
    }
  }
}

/******************************************************************************/

// applies 2nd order HLLC scheme
void MUSCL(double P[ieq][NX+4], double Fhllc1[ieq][NX+4]){


double s1, s2, deltaL1, deltaR1, deltaL2, deltaR2, delta1, delta2;
double primL, primR;

double sl, sr, s_star, U_star_l, U_star_r;
double csr, csl;

double dt;
double max_speed = 0.0;
double cs;
double k;

  for (int i=1; i<= NX+2; i++){
    // SLOPE LIMITER
    for (int iieq=0; iieq <= ieq-1; iieq++){

      deltaL1 = P[iieq][i] - P[iieq][i-1];
      deltaR1 = P[iieq][i+1] - P[iieq][i];

      deltaL2 = P[iieq][i+1] - P[iieq][i];
      deltaR2 = P[iieq][i+2] - P[iieq][i+1];

      s1 = copysign(1.0,deltaL1);
      s2 = copysign(1.0,deltaL2);

      delta1 = s1*max(0.0,min(abs(deltaL1),s1*deltaR1));
      delta2 = s2*max(0.0,min(abs(deltaL2),s2*deltaR2));

      pL[iieq] = P[iieq][i] + 0.5 * delta1;      // nueva primitiva W_L: i
      pR[iieq] = P[iieq][i+1] - 0.5 * delta2;    // nueva primitiva W_R: i+1

    }

    // updates conservatives U from primitives P 

    UUL[0] = pL[0];
    UUL[1] = pL[0]*pL[1];
    UUL[2] = 0.5*pL[0]*pL[1]*pL[1]+pL[2]/(gam-1);

    UUR[0] = pR[0];
    UUR[1] = pR[0]*pR[1];
    UUR[2] = 0.5*pR[0]*pR[1]*pR[1]+pR[2]/(gam-1);

    // new physical fluxes (using updated primitives)
    fL[0] = pL[0]*pL[1];
    fL[1] = pL[0]*pow(pL[1],2)+pL[2];
    fL[2] = pL[1]*(0.5*pL[0]*pow(pL[1],2)+pL[2]*gam/(gam-1));

    fR[0] = pR[0]*pR[1];
    fR[1] = pR[0]*pow(pR[1],2)+pR[2];
    fR[2] = pR[1]*(0.5*pR[0]*pow(pR[1],2)+pR[2]*gam/(gam-1));

    // HLLC

    csl = sqrt(gam*pL[2]/pL[0]);
    csr = sqrt(gam*pR[2]/pR[0]);

    sl = min(pL[1]-csl,pR[1]-csr);
    sr = max(pL[1]+csl,pR[1]+csr);

    s_star = (pR[2]-pL[2]+pL[0]*pL[1]*(sl-pL[1])-pR[0]*pR[1]*(sr-pR[1]))/(pL[0]*(sl-pL[1])-pR[0]*(sr-pR[1]));

    for (int iieq=0; iieq<=ieq-1; iieq++){

      if (iieq == 0){
        U_star_l=pL[0]*(sl-pL[1])/(sl-s_star)*1;
        U_star_r=pR[0]*(sr-pR[1])/(sr-s_star)*1;
      }
      else if (iieq == 1) {
        U_star_l=pL[0]*(sl-pL[1])/(sl-s_star)*s_star;
        U_star_r=pR[0]*(sr-pR[1])/(sr-s_star)*s_star;
      }
      else if (iieq == 2) {
        U_star_l=pL[0]*(sl-pL[1])/(sl-s_star)*((0.5*pL[0]*pL[1]*pL[1]+pL[2]/(gam-1))/pL[0]+(s_star-pL[1])*(s_star+pL[2]/(pL[0]*(sl-pL[1]))));
        U_star_r=pR[0]*(sr-pR[1])/(sr-s_star)*((0.5*pR[0]*pR[1]*pR[1]+pR[2]/(gam-1))/pR[0]+(s_star-pR[1])*(s_star+pR[2]/(pR[0]*(sr-pR[1]))));
      }

      // numerical intercell fluxes
      if (sl >= 0) {
        Fhllc1[iieq][i] = fL[iieq];
      }
      else if ((sl <= 0) && (0 <= s_star)) {
        Fhllc1[iieq][i] = fL[iieq]+sl*(U_star_l-UUL[iieq]);
      }
      else if ((s_star <= 0) && (0 <= sr)) {
        Fhllc1[iieq][i] = fR[iieq]+sr*(U_star_r-UUR[iieq]);
      }
      else if (sr <= 0) {
        Fhllc1[iieq][i] = fR[iieq];
      }
    }
  }
}

/******************************************************************************/

int main() {

  // initial conditions and initializes variables
  initflow(U);

  primitivas1(U,P);

  // writes initial conditions to disk
  output(P);

  // simulation's initial time
  start = clock();
  while (t <= TFIN) {

    primitivas1(U,P);

    // updates time step
    dt = timestep(P);

    printf("%f\n", dt);

    // updates physical fluxes
    fluxes(P, F);

    HLLC(P,U,F,Fhllc1);

    godunov1(U,Fhllc1,UP);

    // applies boundary conditions to UP
    boundary(UP);

    // second sub-time step

    primitivas1(UP,P); // updates UP using P

    MUSCL(P, Fhllc1); // computes numerical fluxes

    godunov2(U, Fhllc1, UP);

    boundary(UP);

    stepviscoso(U,UP);

    // writes to disk
    if (t >= tprint) {
      primitivas1(U,P); // computes primitives
      output(P); // writes to disk primitives
    }

  }

// end
cout << "\nSe calcularon " << it << " iteraciones en "
     << (double)(clock() - start)/CLOCKS_PER_SEC << "s.\n\n";
}
