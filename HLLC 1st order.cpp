// Riemann solver implemented: HLLC along with the Godunov scheme

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <math.h>
#include <time.h>
using namespace std;

/******************************************************************************/

// constant parameters used for the simulation
const int    NX = 400;         // Tamaño de la malla
const double XL = 0.0;         // Coordenada física del extremo izquierdo
const double XR = 1.0;         // Coordenada física del extremo derecho
const double TFIN = 0.2;       // Tiempo final de integración
const double CFL = 0.9;        // Parametro de Courant
const double dtprint = 0.005;   // Intervalo para escribir a disco

const double gam=1.4;

const int ieq=3;

// for the schock-tube: coordinate of separation between initial states
const double X0 = 0.5;

// for advection: wave's velocity
const double A = 1.0;

// Constantes derivadas de las anteriores
const double DX = (XR-XL)/NX;      // space between nodes

// Variables globales
double U[ieq][NX+2];       // "current" conservative variables
double UP[ieq][NX+2];      // "advanced" conservative variables
double F[ieq][NX+2];       // physical fluxes
double P[ieq][NX+2];
double Fhllc[ieq][NX+2];

double csr[ieq][NX+2];
double csl[ieq][NX+2];

double dt;            // time step
double t;          // currrent time
int it;               // current iteration
clock_t start;        // initial time
double tprint;        // time for the following output
int itprint;          // output number


/******************************************************************************/

// sets initial conditions
void initflow(double U[ieq][NX+2]) {

// initialize U in all the domain
// including ghost cells
double x, rho, u, pres;
const double rhol=1.0;
const double ul=0.0;
const double pres_l=1.0;

const double rhor=0.125;
const double ur=0.0;
const double pres_r=0.1;

  for (int i=0; i <= NX+1; i++){
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
void output(double P[ieq][NX+2]) {

  // generates output file's name
  char fname[80];
  sprintf(fname, "Lax_%02i.txt", itprint);

  // opens the file
  fstream fout(fname, ios::out);

  // writes U values to disk
  double x;
  for (int i=1; i <= NX; i++) {
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
void boundary(double U[ieq][NX+2]) {

  for(int iieq=0; iieq<=ieq-1;iieq++){
    U[iieq][0]=U[iieq][1];
    U[iieq][NX+1]=U[iieq][NX];
  }
}

/******************************************************************************/

// computes primitives, including ghost cells
void primitivas(double U[ieq][NX+2], double P[ieq][NX+2]) {

  for (int i=0; i<=NX+1;i++) {
    //F[i]=A*U[i];
    P[0][i]=U[0][i];
    P[1][i]=U[1][i]/U[0][i];
    P[2][i]=(gam-1)*(U[2][i]-pow(U[1][i],2)/(2*U[0][i]));

  }

}

// computes physical fluxes, including ghost cells
void fluxes(double P[ieq][NX+2], double F[ieq][NX+2]) {

  for (int i=0; i<=NX+1;i++) {

    F[0][i]=P[0][i]*P[1][i];
    F[1][i]=P[0][i]*pow(P[1][i],2)+P[2][i];
    F[2][i]=P[1][i]*(0.5*P[0][i]*pow(P[1][i],2)+P[2][i]*gam/(gam-1));
  }

}

/******************************************************************************/

// computes new time step resulting from the CFL condition
double timestep(double P[ieq][NX+2]) {

  double dt;

  // for advection eq., max_u is simply A
  // double max_speed = abs(A);

  // for other cases, we shall compute the maximum value abs(vel)
  double max_speed = 0.0;
  double cs;
  double k;
  for (int i = 1; i <= NX; i++) {
    cs = sqrt(gam*P[2][i]/P[0][i]);
    k = abs(P[1][i])+cs;
    if (k > max_speed) max_speed = k;
  }

  dt = CFL * DX / max_speed;

  return dt;

}

/******************************************************************************/

// applies the HLLC method to obtain the intercell numerical fluxes
void HLLC(double P[ieq][NX+2], double U[ieq][NX+2], double F[ieq][NX+2], double Fhllc[ieq][NX+2]) {

// from 0 to NX

double sl, sr, s_star, U_star_l, U_star_r;
double csr, csl;

/*
  for (int i=0; i<=NX; i++){
    for (int iieq=0; iieq<=ieq-1; iieq++){

      csl[iieq][i] = sqrt(gam*P[2][i]/P[0][i]);
      csr[iieq][i] = sqrt(gam*P[2][i+1]/P[0][i+1]);

      sl[iieq][i] = min(P[1][i]-csl[iieq][i],P[1][i+1]-csr[iieq][i]);
      sr[iieq][i] = max(P[1][i]+csl[iieq][i],P[1][i+1]+csr[iieq][i]);

      s_star[iieq][i] = (P[2][i+1]-P[2][i]+P[0][i]*P[1][i]*(sl-P[1][i])-P[0][i+1]*P[1][i+1]*(sr-P[1][i+1]))/(P[0][i]*(sl-P[1][i])-P[0][i+1]*(sr-P[1][i+1]));

    }
  }

  for (int i=0; i<=NX; i++){

    U_star[0][i]=P[0][i]*(sl[0][i]-P[1][i])/(sl[0][i]-s_star[0][i])*1;
    U_star[1][i]=
    U_star[2][i]=

  }

*/

  for (int i=0; i<=NX; i++){
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


      // intercell numerical fluxes
      if (sl >= 0) {
        Fhllc[iieq][i] = F[iieq][i];
      }
      else if ((sl <= 0) && (0 <= s_star)) {
        Fhllc[iieq][i] = F[iieq][i]+sl*(U_star_l-U[iieq][i]);
      }
      else if ((s_star <= 0) && (0 <= sr)) {
        Fhllc[iieq][i] = F[iieq][i+1]+sr*(U_star_r-U[iieq][i+1]);
      }
      else if (sr <= 0) {
        Fhllc[iieq][i] = F[iieq][i+1];
      }
    }
  }
}

/******************************************************************************/

void godunov(double U[ieq][NX+2], double Fhllc[ieq][NX+2], double UP[ieq][NX+2]) {

  for (int i=1; i<=NX; i++){
    for (int iieq=0; iieq<=ieq-1; iieq++){

      UP[iieq][i]=U[iieq][i]-dt/DX*(Fhllc[iieq][i]-Fhllc[iieq][i-1]);

    }
  }
}

/******************************************************************************/

// this represents one complete time step
void step(double U[ieq][NX+2], double UP[ieq][NX+2]) {

  for (int i = 0; i <= NX+1; i++) {
    for (int iieq=0; iieq<=ieq-1; iieq++){

      U[iieq][i]=UP[iieq][i];

    }
  }

  t = t + dt;
  it = it + 1;

}

/******************************************************************************/

int main() {

  // initial conditions and initializes variables
  initflow(U);

  primitivas(U,P);

  // writes initial conditions to disk
  output(P);

  // simulation's initial time
  start = clock();
  while (t <= TFIN) {

     // updates primitives
    primitivas(U,P);

    // updates time step
    dt = timestep(P);

    // updates phyisical fluxes
    fluxes(P, F);

    // applies HLLC method
    HLLC(P,U,F,Fhllc);
    
    // applies godunov's scheme
    godunov(U,Fhllc,UP);

    // applies boundary conditions to the UP-variables
    boundary(UP);

    // this represents one complete time step
    step(U, UP);

    // writes to disk
    if (t >= tprint) {
      primitivas(U,P);
      output(P);
    }

  }

// end
cout << "\nSe calcularon " << it << " iteraciones en "
     << (double)(clock() - start)/CLOCKS_PER_SEC << "s.\n\n";
}
