#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <complex>
#include <array>
#include <cmath>
#include <cassert>
#include <unistd.h>
#include <omp.h>

#include <boost/numeric/odeint.hpp>
#include "Hydrogen.h"


//#define NDEBUG


using namespace std;
using namespace boost::numeric::odeint;



/* Unit
   Time [ns]
   Frequency [GHz]
   Angular frequency [Grad/s]
   Length [cm]
   Power [GW]
   Total Energy [J] = [W s] = [GW ns]
   Individual Energy [GRad/s]
 */



constexpr complex< double > I(0.0, 1.0);     // Imarinary unit

// Physics Constants
constexpr double H = 6.62606957E-25;         // [J ns]
constexpr double HBAR = 1.054571726E-25;     // [J ns / rad]
constexpr double EPSILON0 = 8.854187817E-14; // [C^2 J-1 cm-1]
constexpr double C = 29.9792458;             // [cm/ns]
constexpr double NDENSITY = 2.6867805E19;    // [cm-3] number density ( 0 degC, 1 atm)


// Unit Conversion Function
double Kayser2GRad( double k /*kayser*/ )
{
  // k * 2 pi c
  return k * 188.36515685l; /*[GRad/s]*/
}

double GRad2Kayser( double omega /*GRad/s*/ )
{
  return omega / 188.36515685l; /*[kayser]*/
}

double NDensity2Amagat( double ndensity /*cm-3*/ )
{
  return ndensity / 2.6867805E19; /*[amagat]*/
}

double Intensity2Strength( double intensity /*GW/cm2*/ )
{
  return sqrt(intensity) * 8.68021098e5; /*[V/cm]*/
}



// Grobal parameters for the simulation
double AbsoluteError = 1E-5; // Absolute error used in Bulirsch-Stoer
double RelativeError = 1E-5; // Relative error used in Bulirsch-Stoer
int Xi = 100; // Half number of space divisions (-Xi -- Xi)
constexpr int EdgeXi = 5; // One side number of space divisions in the Edge( -Xi-6, ..., -Xi-1) ( Xi+1, ..., Xi+6)
int OutputXiPitch = 1;
double SimTimePitch = 0.05; // [ns]
double OutputTimePitch = 0.1; // [ns]
double NextOutputTime = 0;
double SimTime = 20;    // [ns]
constexpr int NVar = 7; // Number of variances in one division ( rho_gg, rho_ee, rho_ge, E^R+(driving), E^L-(driving), E^R-(trigger), E^L+(signal))
ofstream ofs;

enum SpatialAlgorithm { CentralDifference, WENO, };
enum TimeAlgorithm { Dopri5, BulirschStoer, };
enum EquationType { RHS0, RHS1, };




class MaxwellBloch{

private:
  // Target
  const double m_L;           // target length [cm]
  const double m_DeltaZ;      // spatial division pitch [cm]
  const double m_pressure;    // target pressure [Pa]
  const double m_temperature; // target temperature [K]
  const double m_density;     // number density [cm-3]
 
  // Polarizability [C2 cm2 J-2 ns-1 Rad-1] 
  // (Caution: The definition in this program is different fromt the masuda's pdf.  It is multiplied by epsilon0/hbar to simplify the calculation in the iteration.)
  double* m_a;     // polarizability for ground state
  double* m_b;     // polarizability for excited state
  double* m_dp;    // polarizability for the two photon transition

  // Angular Frequency [GRad/s]
  double* m_omega; // angular frequencies for each sideband [Grad/s]
  double m_delta;  // detuning [Grad/s]

  // Relaxation [Grad/s]
  double m_gamma1g;
  double m_gamma1e;
  double m_gamma2;

  // Two Photon Rabi frequency [GRad/s]
  double Rabi_gg;
  double Rabi_ee;
  complex< double > Rabi_ge;

  // Driving laser field
  double m_I_ERp; // peak intensity [GW/cm2]
  double m_I_ELm; // peak intensity [GW/cm2]
  double m_T_ERp; // pulse width (sigma) [ns]
  double m_T_ELm; // pulse width (sigma) [ns]
  double m_D_ERp; // peak time [ns]
  double m_D_ELm; // peak time [ns]

  // Trigger laser field (ERm)
  double m_I_Trig; // peak intensity [GW/cm2]
  double m_T_Trig; // pulse width (sigma) [ns]
  double m_D_Trig; // peak time [ns]

  // Flag
  bool m_resonance;  // Default False : force Omega_gg-Omega_ee+delta=0 if it is true
  bool m_absorption; // Default True  : Absorption term in the MB eq. is ignored if it is false
  bool m_relaxation; // Default True  : Relaxation term in the MB eq. is ignored if it is false
  bool m_squarewave; // Default False : Square electric wave if it is true

  // Algorithm
  EquationType m_equation;
  SpatialAlgorithm m_sAlgorithm;
  TimeAlgorithm m_tAlgorithm;

  void SetOmega( void );
  void CalculatePolarizabilities( void );
  void CalculateRabiFrequency( const complex< double > ERp, const complex< double > ELm, const complex< double > ERm, const complex< double > ELp );
  complex< double > GetERpTerminal( const double time, const int xi=-Xi-1 ); // ER+ at xi==-Xi-1
  complex< double > GetELmTerminal( const double time, const int xi=Xi+1 ); // EL- at xi==Xi+1
  complex< double > GetTrigTerminal( const double time, const int xi=-Xi-1 ); // Trig at xi==-Xi-1

public:
  typedef vector< complex< double > > state_type;
  state_type state;
  const unsigned int NStates;
  const unsigned int NSideBand;

  MaxwellBloch( double length, double pressure=6e4, double temperature=78. ); // [cm], [Pa], [K]
  void operator() ( const state_type &r, state_type &drdt, const double t );
  unsigned int ConvertIndex( int ipar, int xi ) const { return (xi+Xi+EdgeXi)*NVar + ipar; }
  static void Write( const state_type &r, const double t );
  void DumpSetting( void ) const;

  void SetDelta( const double delta );

  void SetERpIntensity( const double intensity ){ m_I_ERp = intensity; };
  void SetELmIntensity( const double intensity ){ m_I_ELm = intensity; };
  void SetERpWidth( const double width ){ m_T_ERp = width; };
  void SetELmWidth( const double width ){ m_T_ELm = width; };
  void SetERpDelay( const double delay ){ m_D_ERp = delay; };
  void SetELmDelay( const double delay ){ m_D_ELm = delay; };
  void SetTrigIntensity( const double intensity ){ m_I_Trig = intensity; };
  void SetTrigWidth( const double width ){ m_T_Trig = width; };
  void SetTrigDelay( const double delay ){ m_D_Trig = delay; };
  void SetSquareWave( const bool flag ){ m_squarewave = flag; };

  void SetResonance( const bool flag ){ m_resonance = flag; };
  void SetAbsorption( const bool flag ){ m_absorption = flag; };
  void SetRelaxation( const bool flag ){ m_relaxation = flag; };
  void SetEquationType( const EquationType eq ){ m_equation = eq; };
  EquationType GetEquationType( void ) const{ return m_equation; };
  void SetSpatialAlgorithm( const SpatialAlgorithm sa ){ m_sAlgorithm = sa; };
  SpatialAlgorithm GetSpatialAlgorithm( void ) const{ return m_sAlgorithm; };
  void SetTimeAlgorithm( const TimeAlgorithm ta ){ m_tAlgorithm = ta; };
  TimeAlgorithm GetTimeAlgorithm( void ) const{ return m_tAlgorithm; };
  

  void SetGammas( const double gamma1g, const double gamma1e, const double gamma2 ){ m_gamma1g=gamma1g; m_gamma1e=gamma1e; m_gamma2=gamma2; } // [Grad/s]

  void DumpRefractiveIndex( const double temp ) const; // [K]

  void rhs0( const state_type &r, state_type &drdt, const double t ); // time : first-derivative; space : first-derivative
  void rhs1( const state_type &r, state_type &drdt, const double t ); // time : first-derivative; space : second-derivative

  inline complex< double > CentralDifference1( const complex< double > &r0, const complex< double > &r1 ){ return (r1 - r0) / (2*m_DeltaZ); };
  inline complex< double > CentralDifference2( const complex< double > &r0, const complex< double > &r1, const complex< double > &r2 ){ return (r2 + r0 - 2.*r1) / pow(m_DeltaZ,2); };
  inline complex< double > WENO1( const complex< double > &u0, const complex< double > &u1, const complex< double > &u2, const complex< double > &u3, const complex< double > &u4, const complex< double > &u5 );
  inline complex< double > WENO2( const complex< double > &u0, const complex< double > &u1, const complex< double > &u2, const complex< double > &u3, const complex< double > &u4, const complex< double > &u5, const complex< double > &u6 );


};


MaxwellBloch::MaxwellBloch( double length, double pressure, double temperature ) : 
  m_L(length), m_DeltaZ(length/2/Xi), m_pressure(pressure), m_temperature(temperature), m_density(NDENSITY/101325.*m_pressure*273.15/m_temperature), NStates(NVar*(2*(Xi+EdgeXi)+1)), NSideBand(1)
{

  state.resize( NStates );

  // Initial conditions
  for( int xi=-Xi-EdgeXi; xi<Xi+EdgeXi+1; xi++ ){
    state[ConvertIndex(0,xi)] = complex< double >( 0.0, 0.0 ); // rho_gg (real)
    state[ConvertIndex(1,xi)] = complex< double >( 0.0, 0.0 ); // rho_ee (real)
    state[ConvertIndex(2,xi)] = complex< double >( 0.0, 0.0 ); // rho_ge (complex)
    state[ConvertIndex(3,xi)] = complex< double >( 0.0, 0.0 ); // E^R+ (complex)
    state[ConvertIndex(4,xi)] = complex< double >( 0.0, 0.0 ); // E^L- (complex)
    state[ConvertIndex(5,xi)] = complex< double >( 0.0, 0.0 ); // E^R- (complex)
    state[ConvertIndex(6,xi)] = complex< double >( 0.0, 0.0 ); // E^L+ (complex)
  }
  for( int xi=-Xi; xi<Xi+1; xi++ ){
    // Hydrogen gas locates only -Xi to +Xi.
    state[ConvertIndex(0,xi)] = complex< double >( 1.0, 0.0 ); // rho_gg (real)
  }
  
  // parameter setting
  m_omega = new double[NSideBand]();
  m_a     = new double[NSideBand]();
  m_b     = new double[NSideBand]();
  m_dp    = new double[NSideBand]();

  m_gamma1g = 2*M_PI * 0.003; // 3 MHz
  m_gamma1e = m_gamma1g;
  m_gamma2  = 2*M_PI * (76/NDensity2Amagat(m_density) + 45.4*NDensity2Amagat(m_density) ) * 1E-3 / 2.; // Hydrogen only (W. Bischel PRA 33 5 3113 (1986))

  m_resonance  = false;
  m_absorption = true;
  m_relaxation = true;
  m_squarewave = false;

  m_equation = RHS0;   // Time : first derivative; Space : first derivative
  m_sAlgorithm = WENO; // WENO
  m_tAlgorithm = BulirschStoer; // Bulirsch-store

  SetDelta( 0 );



}


void MaxwellBloch::operator() ( const state_type &r, state_type &drdt, const double t )
{
  switch( m_equation ){
  case RHS0:
    rhs0( r, drdt, t); // time : first-derivative; space : first-derivative (original MB equation for the adiabatic raman)
    break;
  case RHS1:
    rhs1( r, drdt, t); // time : first-derivative; space : second-derivative (modified MB equation for the soliton paper (PTEP) )
    break;
  }
}


void MaxwellBloch::rhs0( const state_type &r, state_type &drdt, const double t )
// time : first-derivative; space : first-derivative (original MB equation for the adiabatic raman)
{
#pragma omp parallel for
  for( int xi=-Xi-EdgeXi; xi<Xi+EdgeXi+1; xi++ ){
    CalculateRabiFrequency( r[ConvertIndex(3,xi)], r[ConvertIndex(4,xi)], r[ConvertIndex(5,xi)], r[ConvertIndex(6,xi)] );
    
    unsigned int ibase = ConvertIndex( 0, xi );
    
    // population
    drdt[ibase+0] = -2.0 * imag( Rabi_ge * conj(r[ibase+2]) ); // rho_gg w/o decay
    drdt[ibase+1] = -drdt[ibase+0];                            // rho_ee w/o decay
    if( m_relaxation ){
      drdt[ibase+0] += m_gamma1g * r[ibase+1];
      drdt[ibase+1] -= m_gamma1e * r[ibase+1];
    }
    
    // coherence
    drdt[ibase+2] = I * Rabi_ge * ( r[ibase+1]-r[ibase] ); // force resonance w/o decay
    if( !m_resonance ){
      drdt[ibase+2] += I * r[ibase+2] * (Rabi_gg-Rabi_ee+m_delta); // stark shift
    }
    if( m_relaxation ){
      drdt[ibase+2] -= m_gamma2 * r[ibase+2];
    }


    // Electric field
    if( m_sAlgorithm == CentralDifference ){ // Central Difference
      if( xi == -Xi-EdgeXi ){
	// Left terminal
	drdt[ibase+3] = - C * CentralDifference1( GetERpTerminal(t), r[ibase+3+NVar] );
	drdt[ibase+4] =   C * CentralDifference1( r[ibase+4], r[ibase+4+NVar] );
	drdt[ibase+5] = - C * CentralDifference1( GetTrigTerminal(t), r[ibase+5+NVar] );
	drdt[ibase+6] =   C * CentralDifference1( r[ibase+6], r[ibase+6+NVar] );
      }else if( xi == Xi+EdgeXi ){
	// Right terminal
	drdt[ibase+3] = - C * CentralDifference1( r[ibase+3-NVar], r[ibase+3] );
	drdt[ibase+4] =   C * CentralDifference1( r[ibase+4-NVar], GetELmTerminal(t) );
	drdt[ibase+5] = - C * CentralDifference1( r[ibase+5-NVar], r[ibase+5] );
	drdt[ibase+6] =   C * CentralDifference1( r[ibase+6-NVar], 0 );
      }else{
	drdt[ibase+3] = - C * CentralDifference1( r[ibase+3-NVar], r[ibase+3+NVar] );
	drdt[ibase+4] =   C * CentralDifference1( r[ibase+4-NVar], r[ibase+4+NVar] );
	drdt[ibase+5] = - C * CentralDifference1( r[ibase+5-NVar], r[ibase+5+NVar] );
	drdt[ibase+6] =   C * CentralDifference1( r[ibase+6-NVar], r[ibase+6+NVar] );
      }
    }else{ // WENO (linear weights)
      if( xi == -Xi-EdgeXi ){ // Left terminal
	drdt[ibase+3] = - C * WENO1( GetERpTerminal(t,-Xi-EdgeXi-3), GetERpTerminal(t,-Xi-EdgeXi-2), GetERpTerminal(t), r[ibase+3], r[ibase+3+NVar], r[ibase+3+2*NVar] );
	drdt[ibase+4] = - C * WENO1( r[ibase+4+3*NVar], r[ibase+4+2*NVar], r[ibase+4+NVar], r[ibase+4], r[ibase+4], r[ibase+4] );
	drdt[ibase+5] = - C * WENO1( GetTrigTerminal(t,-Xi-EdgeXi-3), GetTrigTerminal(t,-Xi-EdgeXi-2), GetTrigTerminal(t), r[ibase+5], r[ibase+5+NVar], r[ibase+5+2*NVar] );
	drdt[ibase+6] = - C * WENO1( r[ibase+6+3*NVar], r[ibase+6+2*NVar], r[ibase+6+NVar], r[ibase+6], r[ibase+6], r[ibase+6] );
      }else if( xi == -Xi-EdgeXi+1 ){
	drdt[ibase+3] = - C * WENO1( GetERpTerminal(t,-Xi-EdgeXi-2), GetERpTerminal(t), r[ibase+3-NVar], r[ibase+3], r[ibase+3+NVar], r[ibase+3+2*NVar] );
	drdt[ibase+4] = - C * WENO1( r[ibase+4+3*NVar], r[ibase+4+2*NVar], r[ibase+4+NVar], r[ibase+4], r[ibase+4-NVar], r[ibase+4-NVar] );
	drdt[ibase+5] = - C * WENO1( GetTrigTerminal(t,-Xi-EdgeXi-2), GetTrigTerminal(t), r[ibase+5-NVar], r[ibase+5], r[ibase+5+NVar], r[ibase+5+2*NVar] );
	drdt[ibase+6] = - C * WENO1( r[ibase+6+3*NVar], r[ibase+6+2*NVar], r[ibase+6+NVar], r[ibase+6], r[ibase+6-NVar], r[ibase+6-NVar] );
      }else if( xi == -Xi-EdgeXi+2 ){
	drdt[ibase+3] = - C * WENO1( GetERpTerminal(t), r[ibase+3-2*NVar], r[ibase+3-NVar], r[ibase+3], r[ibase+3+NVar], r[ibase+3+2*NVar] );
	drdt[ibase+4] = - C * WENO1( r[ibase+4+3*NVar], r[ibase+4+2*NVar], r[ibase+4+NVar], r[ibase+4], r[ibase+4-NVar], r[ibase+4-2*NVar] );
	drdt[ibase+5] = - C * WENO1( GetTrigTerminal(t), r[ibase+5-2*NVar], r[ibase+5-NVar], r[ibase+5], r[ibase+5+NVar], r[ibase+5+2*NVar] );
	drdt[ibase+6] = - C * WENO1( r[ibase+6+3*NVar], r[ibase+6+2*NVar], r[ibase+6+NVar], r[ibase+6], r[ibase+6-NVar], r[ibase+6-2*NVar] );
      }else if( xi == Xi+EdgeXi-2 ){
	drdt[ibase+3] = - C * WENO1( r[ibase+3-3*NVar], r[ibase+3-2*NVar], r[ibase+3-NVar], r[ibase+3], r[ibase+3+NVar], r[ibase+3+2*NVar] );
	drdt[ibase+4] = - C * WENO1( GetELmTerminal(t), r[ibase+4+2*NVar], r[ibase+4+NVar], r[ibase+4], r[ibase+4-NVar], r[ibase+4-2*NVar] );
	drdt[ibase+5] = - C * WENO1( r[ibase+5-3*NVar], r[ibase+5-2*NVar], r[ibase+5-NVar], r[ibase+5], r[ibase+5+NVar], r[ibase+5+2*NVar] );
	drdt[ibase+6] = - C * WENO1( 0, r[ibase+6+2*NVar], r[ibase+6+NVar], r[ibase+6], r[ibase+6-NVar], r[ibase+6-2*NVar] );
      }else if( xi == Xi+EdgeXi-1 ){
	drdt[ibase+3] = - C * WENO1( r[ibase+3-3*NVar], r[ibase+3-2*NVar], r[ibase+3-NVar], r[ibase+3], r[ibase+3+NVar], r[ibase+3+NVar] );
	drdt[ibase+4] = - C * WENO1( GetELmTerminal(t,Xi+EdgeXi+2), GetELmTerminal(t), r[ibase+4+NVar], r[ibase+4], r[ibase+4-NVar], r[ibase+4-2*NVar] );
	drdt[ibase+5] = - C * WENO1( r[ibase+5-3*NVar], r[ibase+5-2*NVar], r[ibase+5-NVar], r[ibase+5], r[ibase+5+NVar], r[ibase+5+NVar] );
	drdt[ibase+6] = - C * WENO1( 0, 0, r[ibase+6+NVar], r[ibase+6], r[ibase+6-NVar], r[ibase+6-2*NVar] );
      }else if( xi == Xi+EdgeXi ){ // Right terminal
	drdt[ibase+3] = - C * WENO1( r[ibase+3-3*NVar], r[ibase+3-2*NVar], r[ibase+3-NVar], r[ibase+3], r[ibase+3], r[ibase+3] );
	drdt[ibase+4] = - C * WENO1( GetELmTerminal(t,Xi+EdgeXi+3), GetELmTerminal(t,Xi+EdgeXi+2), GetELmTerminal(t), r[ibase+4], r[ibase+4-NVar], r[ibase+4-2*NVar] );
	drdt[ibase+5] = - C * WENO1( r[ibase+5-3*NVar], r[ibase+5-2*NVar], r[ibase+5-NVar], r[ibase+5], r[ibase+5], r[ibase+5] );
	drdt[ibase+6] = - C * WENO1( 0, 0, 0, r[ibase+6], r[ibase+6-NVar], r[ibase+6-2*NVar] );
      }else{
	drdt[ibase+3] = - C * WENO1( r[ibase+3-3*NVar], r[ibase+3-2*NVar], r[ibase+3-NVar], r[ibase+3], r[ibase+3+NVar], r[ibase+3+2*NVar] );
	drdt[ibase+4] = - C * WENO1( r[ibase+4+3*NVar], r[ibase+4+2*NVar], r[ibase+4+NVar], r[ibase+4], r[ibase+4-NVar], r[ibase+4-2*NVar] );
	drdt[ibase+5] = - C * WENO1( r[ibase+5-3*NVar], r[ibase+5-2*NVar], r[ibase+5-NVar], r[ibase+5], r[ibase+5+NVar], r[ibase+5+2*NVar] );
	drdt[ibase+6] = - C * WENO1( r[ibase+6+3*NVar], r[ibase+6+2*NVar], r[ibase+6+NVar], r[ibase+6], r[ibase+6-NVar], r[ibase+6-2*NVar] );
      }
    }
    if( m_absorption ){
      drdt[ibase+3] += I*m_omega[0]*m_density * (HBAR/EPSILON0) * ( (r[ibase]*m_a[0]+r[ibase+1]*m_b[0])*r[ibase+3] + 2.*m_dp[0]*conj(r[ibase+2]*r[ibase+4]) );
      drdt[ibase+4] += I*m_omega[0]*m_density * (HBAR/EPSILON0) * ( (r[ibase]*m_a[0]+r[ibase+1]*m_b[0])*r[ibase+4] + 2.*m_dp[0]*conj(r[ibase+2]*r[ibase+3]) );
      drdt[ibase+5] += I*m_omega[0]*m_density * (HBAR/EPSILON0) * ( (r[ibase]*m_a[0]+r[ibase+1]*m_b[0])*r[ibase+5] + 2.*m_dp[0]*conj(r[ibase+2]*r[ibase+6]) );
      drdt[ibase+6] += I*m_omega[0]*m_density * (HBAR/EPSILON0) * ( (r[ibase]*m_a[0]+r[ibase+1]*m_b[0])*r[ibase+6] + 2.*m_dp[0]*conj(r[ibase+2]*r[ibase+5]) );
    }

  }
}



void MaxwellBloch::rhs1( const state_type &r, state_type &drdt, const double t )
// time : first-derivative; space : second-derivative (modified MB equation for the soliton paper (PTEP) )
{
#pragma omp parallel for
  for( int xi=-Xi-EdgeXi; xi<Xi+EdgeXi+1; xi++ ){
    CalculateRabiFrequency( r[ConvertIndex(3,xi)], r[ConvertIndex(4,xi)], r[ConvertIndex(5,xi)], r[ConvertIndex(6,xi)] );
    
    unsigned int ibase = ConvertIndex( 0, xi );
    
    // population
    drdt[ibase+0] = -2.0 * imag( Rabi_ge * conj(r[ibase+2]) ); // rho_gg w/o decay
    drdt[ibase+1] = -drdt[ibase+0];                            // rho_ee w/o decay
    if( m_relaxation ){
      drdt[ibase+0] += m_gamma1g * r[ibase+1];
      drdt[ibase+1] -= m_gamma1e * r[ibase+1];
    }
    
    // coherence
    drdt[ibase+2] = I * Rabi_ge * ( r[ibase+1]-r[ibase] ); // force resonance w/o decay
    if( !m_resonance ){
      drdt[ibase+2] += I * r[ibase+2] * (Rabi_gg-Rabi_ee+m_delta); // stark shift
    }
    if( m_relaxation ){
      drdt[ibase+2] -= m_gamma2 * r[ibase+2];
    }


    // Electric field
    if( m_sAlgorithm == CentralDifference ){ // Central Difference
      if( xi == -Xi-EdgeXi ){
	// Left terminal
	drdt[ibase+3] = - C * CentralDifference1( GetERpTerminal(t), r[ibase+3+NVar] )  + I*C*C/(2.*m_omega[0]) * CentralDifference2( GetERpTerminal(t), r[ibase+3], r[ibase+3+NVar]);
	drdt[ibase+4] =   C * CentralDifference1( r[ibase+4], r[ibase+4+NVar] )         + I*C*C/(2.*m_omega[0]) * CentralDifference2( r[ibase+4], r[ibase+4], r[ibase+4+NVar] );
	drdt[ibase+5] = - C * CentralDifference1( GetTrigTerminal(t), r[ibase+5+NVar] ) + I*C*C/(2.*m_omega[0]) * CentralDifference2( GetTrigTerminal(t), r[ibase+5], r[ibase+5+NVar] );
	drdt[ibase+6] =   C * CentralDifference1( r[ibase+6], r[ibase+6+NVar] )         + I*C*C/(2.*m_omega[0]) * CentralDifference2( r[ibase+6], r[ibase+6], r[ibase+6+NVar] );
      }else if( xi == Xi+EdgeXi ){
	// Right terminal
	drdt[ibase+3] = - C * CentralDifference1( r[ibase+3-NVar], r[ibase+3] )        + I*C*C/(2.*m_omega[0]) * CentralDifference2( r[ibase+3-NVar], r[ibase+3], r[ibase+3] );
	drdt[ibase+4] =   C * CentralDifference1( r[ibase+4-NVar], GetELmTerminal(t) ) + I*C*C/(2.*m_omega[0]) * CentralDifference2( r[ibase+4-NVar], r[ibase+4], GetELmTerminal(t) );
	drdt[ibase+5] = - C * CentralDifference1( r[ibase+5-NVar], r[ibase+5] )        + I*C*C/(2.*m_omega[0]) * CentralDifference2( r[ibase+5-NVar], r[ibase+5], r[ibase+5] );
	drdt[ibase+6] =   C * CentralDifference1( r[ibase+6-NVar], 0 )                 + I*C*C/(2.*m_omega[0]) * CentralDifference2( r[ibase+6-NVar], r[ibase+6], 0 );
      }else{
	drdt[ibase+3] = - C * CentralDifference1( r[ibase+3-NVar], r[ibase+3+NVar] ) + I*C*C/(2.*m_omega[0]) * CentralDifference2( r[ibase+3-NVar], r[ibase+3], r[ibase+3+NVar] );
	drdt[ibase+4] =   C * CentralDifference1( r[ibase+4-NVar], r[ibase+4+NVar] ) + I*C*C/(2.*m_omega[0]) * CentralDifference2( r[ibase+4-NVar], r[ibase+4], r[ibase+4+NVar] );
	drdt[ibase+5] = - C * CentralDifference1( r[ibase+5-NVar], r[ibase+5+NVar] ) + I*C*C/(2.*m_omega[0]) * CentralDifference2( r[ibase+5-NVar], r[ibase+5], r[ibase+5+NVar] );
	drdt[ibase+6] =   C * CentralDifference1( r[ibase+6-NVar], r[ibase+6+NVar] ) + I*C*C/(2.*m_omega[0]) * CentralDifference2( r[ibase+6-NVar], r[ibase+6], r[ibase+6+NVar] );
      }
    }else{ // WENO (linear weight)
      if( xi == -Xi-EdgeXi ){ // Left terminal
	drdt[ibase+3] = - C * WENO1( GetERpTerminal(t,-Xi-EdgeXi-3), GetERpTerminal(t,-Xi-EdgeXi-2), GetERpTerminal(t), r[ibase+3], r[ibase+3+NVar], r[ibase+3+2*NVar] )
	  + I*C*C/(2.*m_omega[0]) * WENO2( GetERpTerminal(t,-Xi-EdgeXi-3), GetERpTerminal(t,-Xi-EdgeXi-2), GetERpTerminal(t), r[ibase+3], r[ibase+3+NVar], r[ibase+3+2*NVar], r[ibase+3+3*NVar] );
	drdt[ibase+4] = - C * WENO1( r[ibase+4+3*NVar], r[ibase+4+2*NVar], r[ibase+4+NVar], r[ibase+4], r[ibase+4], r[ibase+4] )
	  + I*C*C/(2.*m_omega[0]) * WENO2( r[ibase+4+3*NVar], r[ibase+4+2*NVar], r[ibase+4+NVar], r[ibase+4], r[ibase+4], r[ibase+4], r[ibase+4] );
	drdt[ibase+5] = - C * WENO1( GetTrigTerminal(t,-Xi-EdgeXi-3), GetTrigTerminal(t,-Xi-EdgeXi-2), GetTrigTerminal(t), r[ibase+5], r[ibase+5+NVar], r[ibase+5+2*NVar] )
	  + I*C*C/(2.*m_omega[0]) * WENO2( GetTrigTerminal(t,-Xi-EdgeXi-3), GetTrigTerminal(t,-Xi-EdgeXi-2), GetTrigTerminal(t), r[ibase+5], r[ibase+5+NVar], r[ibase+5+2*NVar], r[ibase+5+3*NVar] );
	drdt[ibase+6] = - C * WENO1( r[ibase+6+3*NVar], r[ibase+6+2*NVar], r[ibase+6+NVar], r[ibase+6], r[ibase+6], r[ibase+6] )
	  + I*C*C/(2.*m_omega[0]) * WENO2( r[ibase+6+3*NVar], r[ibase+6+2*NVar], r[ibase+6+NVar], r[ibase+6], r[ibase+6], r[ibase+6], r[ibase+6] );
      }else if( xi == -Xi-EdgeXi+1 ){
	drdt[ibase+3] = - C * WENO1( GetERpTerminal(t,-Xi-EdgeXi-2), GetERpTerminal(t), r[ibase+3-NVar], r[ibase+3], r[ibase+3+NVar], r[ibase+3+2*NVar] )
	  + I*C*C/(2.*m_omega[0]) * WENO2( GetERpTerminal(t,-Xi-EdgeXi-2), GetERpTerminal(t), r[ibase+3-NVar], r[ibase+3], r[ibase+3+NVar], r[ibase+3+2*NVar], r[ibase+3+3*NVar] );
	drdt[ibase+4] = - C * WENO1( r[ibase+4+3*NVar], r[ibase+4+2*NVar], r[ibase+4+NVar], r[ibase+4], r[ibase+4-NVar], r[ibase+4-NVar] )
	  + I*C*C/(2.*m_omega[0]) * WENO2( r[ibase+4+3*NVar], r[ibase+4+2*NVar], r[ibase+4+NVar], r[ibase+4], r[ibase+4-NVar], r[ibase+4-NVar], r[ibase+4-NVar] );
	drdt[ibase+5] = - C * WENO1( GetTrigTerminal(t,-Xi-EdgeXi-2), GetTrigTerminal(t), r[ibase+5-NVar], r[ibase+5], r[ibase+5+NVar], r[ibase+5+2*NVar] )
	  + I*C*C/(2.*m_omega[0]) * WENO2( GetTrigTerminal(t,-Xi-EdgeXi-2), GetTrigTerminal(t), r[ibase+5-NVar], r[ibase+5], r[ibase+5+NVar], r[ibase+5+2*NVar], r[ibase+5+3*NVar] );
	drdt[ibase+6] = - C * WENO1( r[ibase+6+3*NVar], r[ibase+6+2*NVar], r[ibase+6+NVar], r[ibase+6], r[ibase+6-NVar], r[ibase+6-NVar] )
	  + I*C*C/(2.*m_omega[0]) * WENO2( r[ibase+6+3*NVar], r[ibase+6+2*NVar], r[ibase+6+NVar], r[ibase+6], r[ibase+6-NVar], r[ibase+6-NVar], r[ibase+6-NVar] );
      }else if( xi == -Xi-EdgeXi+2 ){
	drdt[ibase+3] = - C * WENO1( GetERpTerminal(t), r[ibase+3-2*NVar], r[ibase+3-NVar], r[ibase+3], r[ibase+3+NVar], r[ibase+3+2*NVar] )
	  + I*C*C/(2.*m_omega[0]) * WENO2( GetERpTerminal(t), r[ibase+3-2*NVar], r[ibase+3-NVar], r[ibase+3], r[ibase+3+NVar], r[ibase+3+2*NVar], r[ibase+3+3*NVar] );
	drdt[ibase+4] = - C * WENO1( r[ibase+4+3*NVar], r[ibase+4+2*NVar], r[ibase+4+NVar], r[ibase+4], r[ibase+4-NVar], r[ibase+4-2*NVar] )
	  + I*C*C/(2.*m_omega[0]) * WENO2( r[ibase+4+3*NVar], r[ibase+4+2*NVar], r[ibase+4+NVar], r[ibase+4], r[ibase+4-NVar], r[ibase+4-2*NVar], r[ibase+4-2*NVar] );
	drdt[ibase+5] = - C * WENO1( GetTrigTerminal(t), r[ibase+5-2*NVar], r[ibase+5-NVar], r[ibase+5], r[ibase+5+NVar], r[ibase+5+2*NVar] )
	  + I*C*C/(2.*m_omega[0]) * WENO2( GetTrigTerminal(t), r[ibase+5-2*NVar], r[ibase+5-NVar], r[ibase+5], r[ibase+5+NVar], r[ibase+5+2*NVar], r[ibase+5+3*NVar] );
	drdt[ibase+6] = - C * WENO1( r[ibase+6+3*NVar], r[ibase+6+2*NVar], r[ibase+6+NVar], r[ibase+6], r[ibase+6-NVar], r[ibase+6-2*NVar] )
	  + I*C*C/(2.*m_omega[0]) * WENO2( r[ibase+6+3*NVar], r[ibase+6+2*NVar], r[ibase+6+NVar], r[ibase+6], r[ibase+6-NVar], r[ibase+6-2*NVar], r[ibase+6-2*NVar] );
      }else if( xi == Xi+EdgeXi-2 ){
	drdt[ibase+3] = - C * WENO1( r[ibase+3-3*NVar], r[ibase+3-2*NVar], r[ibase+3-NVar], r[ibase+3], r[ibase+3+NVar], r[ibase+3+2*NVar] )
	  + I*C*C/(2.*m_omega[0]) * WENO2( r[ibase+3-3*NVar], r[ibase+3-2*NVar], r[ibase+3-NVar], r[ibase+3], r[ibase+3+NVar], r[ibase+3+2*NVar], r[ibase+3+2*NVar] );
	drdt[ibase+4] = - C * WENO1( GetELmTerminal(t), r[ibase+4+2*NVar], r[ibase+4+NVar], r[ibase+4], r[ibase+4-NVar], r[ibase+4-2*NVar] )
	  + I*C*C/(2.*m_omega[0]) * WENO2( GetELmTerminal(t), r[ibase+4+2*NVar], r[ibase+4+NVar], r[ibase+4], r[ibase+4-NVar], r[ibase+4-2*NVar], r[ibase+4-3*NVar] );
	drdt[ibase+5] = - C * WENO1( r[ibase+5-3*NVar], r[ibase+5-2*NVar], r[ibase+5-NVar], r[ibase+5], r[ibase+5+NVar], r[ibase+5+2*NVar] )
	  + I*C*C/(2.*m_omega[0]) * WENO2( r[ibase+5-3*NVar], r[ibase+5-2*NVar], r[ibase+5-NVar], r[ibase+5], r[ibase+5+NVar], r[ibase+5+2*NVar], r[ibase+5+2*NVar] );
	drdt[ibase+6] = - C * WENO1( 0, r[ibase+6+2*NVar], r[ibase+6+NVar], r[ibase+6], r[ibase+6-NVar], r[ibase+6-2*NVar] )
	  + I*C*C/(2.*m_omega[0]) * WENO2( 0, r[ibase+6+2*NVar], r[ibase+6+NVar], r[ibase+6], r[ibase+6-NVar], r[ibase+6-2*NVar], r[ibase+6-3*NVar] );
      }else if( xi == Xi+EdgeXi-1 ){
	drdt[ibase+3] = - C * WENO1( r[ibase+3-3*NVar], r[ibase+3-2*NVar], r[ibase+3-NVar], r[ibase+3], r[ibase+3+NVar], r[ibase+3+NVar] )
	  + I*C*C/(2.*m_omega[0]) * WENO2( r[ibase+3-3*NVar], r[ibase+3-2*NVar], r[ibase+3-NVar], r[ibase+3], r[ibase+3+NVar], r[ibase+3+NVar], r[ibase+3+NVar] );
	drdt[ibase+4] = - C * WENO1( GetELmTerminal(t,Xi+EdgeXi+2), GetELmTerminal(t), r[ibase+4+NVar], r[ibase+4], r[ibase+4-NVar], r[ibase+4-2*NVar] )
	  + I*C*C/(2.*m_omega[0]) * WENO2( GetELmTerminal(t,Xi+EdgeXi+2), GetELmTerminal(t), r[ibase+4+NVar], r[ibase+4], r[ibase+4-NVar], r[ibase+4-2*NVar], r[ibase+4-3*NVar] );
	drdt[ibase+5] = - C * WENO1( r[ibase+5-3*NVar], r[ibase+5-2*NVar], r[ibase+5-NVar], r[ibase+5], r[ibase+5+NVar], r[ibase+5+NVar] )
	  + I*C*C/(2.*m_omega[0]) * WENO2( r[ibase+5-3*NVar], r[ibase+5-2*NVar], r[ibase+5-NVar], r[ibase+5], r[ibase+5+NVar], r[ibase+5+NVar], r[ibase+5+NVar] );
	drdt[ibase+6] = - C * WENO1( 0, 0, r[ibase+6+NVar], r[ibase+6], r[ibase+6-NVar], r[ibase+6-2*NVar] )
	  + I*C*C/(2.*m_omega[0]) * WENO2( 0, 0, r[ibase+6+NVar], r[ibase+6], r[ibase+6-NVar], r[ibase+6-2*NVar], r[ibase+6-3*NVar] );
      }else if( xi == Xi+EdgeXi ){ // Right terminal
	drdt[ibase+3] = - C * WENO1( r[ibase+3-3*NVar], r[ibase+3-2*NVar], r[ibase+3-NVar], r[ibase+3], r[ibase+3], r[ibase+3] )
	  + I*C*C/(2.*m_omega[0]) * WENO2( r[ibase+3-3*NVar], r[ibase+3-2*NVar], r[ibase+3-NVar], r[ibase+3], r[ibase+3], r[ibase+3], r[ibase+3] );
	drdt[ibase+4] = - C * WENO1( GetELmTerminal(t,Xi+EdgeXi+3), GetELmTerminal(t,Xi+EdgeXi+2), GetELmTerminal(t), r[ibase+4], r[ibase+4-NVar], r[ibase+4-2*NVar] )
	  + I*C*C/(2.*m_omega[0]) * WENO2( GetELmTerminal(t,Xi+EdgeXi+3), GetELmTerminal(t,Xi+EdgeXi+2), GetELmTerminal(t), r[ibase+4], r[ibase+4-NVar], r[ibase+4-2*NVar], r[ibase+4-3*NVar] );
	drdt[ibase+5] = - C * WENO1( r[ibase+5-3*NVar], r[ibase+5-2*NVar], r[ibase+5-NVar], r[ibase+5], r[ibase+5], r[ibase+5] )
	  + I*C*C/(2.*m_omega[0]) * WENO2( r[ibase+5-3*NVar], r[ibase+5-2*NVar], r[ibase+5-NVar], r[ibase+5], r[ibase+5], r[ibase+5], r[ibase+5] );
	drdt[ibase+6] = - C * WENO1( 0, 0, 0, r[ibase+6], r[ibase+6-NVar], r[ibase+6-2*NVar] )
	  + I*C*C/(2.*m_omega[0]) * WENO2(  0, 0, 0, r[ibase+6], r[ibase+6-NVar], r[ibase+6-2*NVar], r[ibase+6-3*NVar] );
      }else{
	drdt[ibase+3] = - C * WENO1( r[ibase+3-3*NVar], r[ibase+3-2*NVar], r[ibase+3-NVar], r[ibase+3], r[ibase+3+NVar], r[ibase+3+2*NVar] )
	  + I*C*C/(2.*m_omega[0]) * WENO2( r[ibase+3-3*NVar], r[ibase+3-2*NVar], r[ibase+3-NVar], r[ibase+3], r[ibase+3+NVar], r[ibase+3+2*NVar], r[ibase+3+3*NVar] );
	drdt[ibase+4] = - C * WENO1( r[ibase+4+3*NVar], r[ibase+4+2*NVar], r[ibase+4+NVar], r[ibase+4], r[ibase+4-NVar], r[ibase+4-2*NVar] )
	  + I*C*C/(2.*m_omega[0]) * WENO2(  r[ibase+4+3*NVar], r[ibase+4+2*NVar], r[ibase+4+NVar], r[ibase+4], r[ibase+4-NVar], r[ibase+4-2*NVar], r[ibase+4-3*NVar] );
	drdt[ibase+5] = - C * WENO1( r[ibase+5-3*NVar], r[ibase+5-2*NVar], r[ibase+5-NVar], r[ibase+5], r[ibase+5+NVar], r[ibase+5+2*NVar] )
	  + I*C*C/(2.*m_omega[0]) * WENO2( r[ibase+5-3*NVar], r[ibase+5-2*NVar], r[ibase+5-NVar], r[ibase+5], r[ibase+5+NVar], r[ibase+5+2*NVar], r[ibase+5+3*NVar] );
	drdt[ibase+6] = - C * WENO1( r[ibase+6+3*NVar], r[ibase+6+2*NVar], r[ibase+6+NVar], r[ibase+6], r[ibase+6-NVar], r[ibase+6-2*NVar] )
	  + I*C*C/(2.*m_omega[0]) * WENO2(  r[ibase+6+3*NVar], r[ibase+6+2*NVar], r[ibase+6+NVar], r[ibase+6], r[ibase+6-NVar], r[ibase+6-2*NVar], r[ibase+6-3*NVar] );
      }
    }
    if( m_absorption ){
      drdt[ibase+3] += I*m_omega[0]*m_density * (HBAR/EPSILON0) * ( (r[ibase]*m_a[0]+r[ibase+1]*m_b[0])*r[ibase+3] + 2.*m_dp[0]*conj(r[ibase+2]*r[ibase+4]) );
      drdt[ibase+4] += I*m_omega[0]*m_density * (HBAR/EPSILON0) * ( (r[ibase]*m_a[0]+r[ibase+1]*m_b[0])*r[ibase+4] + 2.*m_dp[0]*conj(r[ibase+2]*r[ibase+3]) );
      drdt[ibase+5] += I*m_omega[0]*m_density * (HBAR/EPSILON0) * ( (r[ibase]*m_a[0]+r[ibase+1]*m_b[0])*r[ibase+5] + 2.*m_dp[0]*conj(r[ibase+2]*r[ibase+6]) );
      drdt[ibase+6] += I*m_omega[0]*m_density * (HBAR/EPSILON0) * ( (r[ibase]*m_a[0]+r[ibase+1]*m_b[0])*r[ibase+6] + 2.*m_dp[0]*conj(r[ibase+2]*r[ibase+5]) );
    }


  }
}






complex< double > MaxwellBloch::WENO1( const complex< double > &u0, const complex< double > &u1, const complex< double > &u2, const complex< double > &u3, const complex< double > &u4, const complex< double > &u5 )
{
  /* WENO (linear weight) for the first derivative equation
    C. W. Shu, SIAM Review, 51:82--126, 2009 */

  /* u0: u_{i-3}; u1: u_{i-2}; u2: u_{i-1}; u3: u_{i}; u4: u_{i+1}; u5: u_{i+2} */
  //complex< double > upward = 3./128.*u1 - 5./32.*u2 + 45./64.*u3 + 15./32.*u4 - 5./128.*u5;
  //complex< double > downward = 3./128.*u0 - 5./32.*u1 + 45./64.*u2 + 15./32.*u3 - 5./128.*u4;
  //return (upward - downward) / m_DeltaZ;

  return ( -3.*u0 + 23.*u1 -110.*u2 + 30.*u3 + 65.*u4 - 5.*u5 ) * 7.8125e-3 / m_DeltaZ;

}


complex< double > MaxwellBloch::WENO2( const complex< double > &u0, const complex< double > &u1, const complex< double > &u2, const complex< double > &u3, const complex< double > &u4, const complex< double > &u5, const complex< double > &u6 )
{
  /* WENO (linear wegith) for the second derivative equation
     Y. Liu, C. W. Shu, and M. Zhang, SIAM J. Sci. Comput. 33 2 939--965 (2011) */

  /* u0: u_{i-3}; u1: u_{i-2}; u2: u_{i-1}; u3: u_{i}; u4: u_{i+1}; u5: u_{i+2}; u6: u_{i+3}; */
  return ( 2.*u0 - 27.*u1 + 270.*u2 - 490.*u3 + 270.*u4 - 27.*u5 + 2.*u6 ) / ( 180. * pow(m_DeltaZ,2) );
}





void MaxwellBloch::Write( const state_type &r, const double t )
{
  if( fabs(NextOutputTime-t) > SimTimePitch/2. ) return;
  NextOutputTime += OutputTimePitch;

  ofs << t;
  

  for( int i=0; i<NVar; i++ ){ // left edge point
    ofs << '\t' << r[i].real() << '\t' << r[i].imag();
  }
  for( int ipos = EdgeXi; ipos<2*Xi+EdgeXi+1; ipos+=OutputXiPitch ){ // target inside
    for( int i = ipos*NVar; i<ipos*NVar+NVar; i++ ){
      ofs << '\t' << r[i].real() << '\t' << r[i].imag();
    }
  }
  for( int i=2*(Xi+EdgeXi)*NVar; i<2*(Xi+EdgeXi)*NVar+NVar; i++ ){ // right edge point
    ofs << '\t' << r[i].real() << '\t' << r[i].imag();
  }
  ofs << '\n' << flush;
  cerr << t << " ns\r" << flush;
}


void MaxwellBloch::DumpSetting( void ) const
{
  ofs.setf(ios::scientific);
  ofs.precision(10);

  ofs << "# TargetLength\t\t" << m_L << " cm" << endl;
  ofs << "# DeltaZ\t\t" << m_DeltaZ << " cm" << endl;
  ofs << "# Xi\t\t" << Xi << endl;
  ofs << "# OutputXiPitch\t\t" << OutputXiPitch << endl;
  ofs << "# EdgeXi\t\t" << EdgeXi << endl;
  ofs << "# SimulationTime\t\t" << SimTime << " ns" << endl;
  ofs << "# SimTimePitch\t\t" << SimTimePitch << " ns" << endl;
  ofs << "# OutputTimePitch\t\t" << OutputTimePitch << " ns" << endl;

  ofs << "# Pressure\t\t" << m_pressure << " Pa" << endl;
  ofs << "# Temperature\t\t" << m_temperature << " K" << endl;
  ofs << "# NumberDensity\t\t" << m_density << " cm-3" << endl;
  
  for( unsigned int q=0; q<NSideBand; q++ ){
    ofs << "# omega" << q << "\t\t" << m_omega[q] << " GRad/s" << endl;
  }
  ofs << "# Detuning\t\t" << m_delta << " GRad/s" << endl;

  ofs << "# a0\t\t" << m_a[0] << endl;
  ofs << "# b0\t\t" << m_b[0] << endl;
  ofs << "# dp0\t\t" << m_dp[0] << endl;

  ofs << "# Gamma1g\t\t" << m_gamma1g << " Grad/s" << endl;
  ofs << "# Gamma1e\t\t" << m_gamma1e << " Grad/s" << endl;
  ofs << "# Gamma2\t\t" << m_gamma2 << " Grad/s" << endl;

  ofs << "# ERpIntensity\t\t" << m_I_ERp << " GW/cm2" << endl;
  ofs << "# ELmIntensity\t\t" << m_I_ELm << " GW/cm2" << endl;
  ofs << "# ERpWidth (sigma)\t\t" << m_T_ERp << " ns" << endl;
  ofs << "# ELmWidth (sigma)\t\t" << m_T_ELm << " ns" << endl;
  ofs << "# ERpDelay\t\t" << m_D_ERp << " ns" << endl;
  ofs << "# ELmDelay\t\t" << m_D_ELm << " ns" << endl;

  ofs << "# TrigIntensity\t\t" << m_I_Trig << " GW/cm2" << endl;
  ofs << "# TrigWidth (sigma)\t\t" << m_T_Trig << " ns" << endl;
  ofs << "# TrigDelay\t\t" << m_D_Trig << " ns" << endl;
  
  ofs << "# ForceResonance\t\t" << m_resonance << endl;
  ofs << "# Absorption\t\t" << m_absorption << endl;
  ofs << "# Relaxation\t\t" << m_relaxation << endl;
  ofs << "# Square\t\t" << m_squarewave << endl;

  ofs << "# AbsoluteError\t\t" << AbsoluteError << endl;
  ofs << "# RelativeError\t\t" << RelativeError << endl;

  ofs << "# Equation\t\t" << m_equation << endl;
  ofs << "# SpatialAlgorithm\t" << m_sAlgorithm << endl;
  ofs << "# TimeAlgorithm\t" << m_tAlgorithm << endl;

#ifdef _OPENMP
  ofs << "# OpenMP\t\t" << omp_get_max_threads() << endl;
#else
  ofs << "# OpenMP\t\t0" << endl;
#endif

  ofs << "# Parameter Setting" << endl;
}

void MaxwellBloch::SetDelta( const double delta )
{ 
  m_delta = delta;
  SetOmega();
}


void MaxwellBloch::SetOmega( void )
{
  m_omega[0] = ( Kayser2GRad( Hydrogen::m_omega_eg_k ) - m_delta ) / 2.;
  CalculatePolarizabilities();
  
}


void MaxwellBloch::CalculatePolarizabilities( void )
{

  for( unsigned int q=0; q<NSideBand; q++ ){
    double ws = GRad2Kayser( m_omega[q] );
    double sa=0, sb=0, sd=0;

    for( int j=0; j<Hydrogen::nOmegaB; j++ ){
      sa += Hydrogen::m_uaB[j] * ( 1./(Hydrogen::m_omegaB[j]-ws) + 1./(Hydrogen::m_omegaB[j]+ws) );
      sb += Hydrogen::m_ubB[j] * ( 1./(Hydrogen::m_omegaB[j]-Hydrogen::m_omega_eg_k-ws) + 1./(Hydrogen::m_omegaB[j]-Hydrogen::m_omega_eg_k+ws) );
      sd += sqrt( Hydrogen::m_uaB[j] * Hydrogen::m_ubB[j] ) / ( Hydrogen::m_omegaB[j]-Hydrogen::m_omega_eg_k+ws) * Hydrogen::m_saB[j] * Hydrogen::m_sbB[j];
    }
    for( int j=0; j<Hydrogen::nOmegaC; j++ ){
      sa += Hydrogen::m_uaC[j] * ( 1./(Hydrogen::m_omegaC[j]-ws) + 1./(Hydrogen::m_omegaC[j]+ws) );
      sb += Hydrogen::m_ubC[j] * ( 1./(Hydrogen::m_omegaC[j]-Hydrogen::m_omega_eg_k-ws) + 1./(Hydrogen::m_omegaC[j]-Hydrogen::m_omega_eg_k+ws) );
      sd += sqrt( Hydrogen::m_uaC[j] * Hydrogen::m_ubC[j] ) / ( Hydrogen::m_omegaC[j]-Hydrogen::m_omega_eg_k+ws) * Hydrogen::m_saC[j] * Hydrogen::m_sbC[j];
    }

    /*
    The unit of |d|^2 in my equation and that of u in Fortrun program are different. ( C^2/cm^2 or C^2/m^2 )
    The unit of the polarizability (a_0, b_0) is different. 
    */


    double ConversionFactor = 1e4; // [C2 m2 -> C2 cm2]
    ConversionFactor *= 1. / (2*M_PI*C) / (2*HBAR*HBAR); // Kayser -> GRad/s (2 pi c), 2hbar^2 is
    //ConversionFactor *= 6. / M_PI;  // Empirical?


    m_a[q] = sa * ConversionFactor * sqrt(3); // [C2 cm2 J-2 ns-1 Rad-1]
    m_b[q] = sb * ConversionFactor * sqrt(3);
    m_dp[q] = sd * ConversionFactor;
  }
   
}




void MaxwellBloch::CalculateRabiFrequency( const complex< double > ERp, const complex< double > ELm, const complex< double > ERm, const complex< double > ELp )
{

  Rabi_gg = 0.5 * m_a[0] * ( norm(ERp) + norm(ELm) + norm(ERm) + norm(ELp) ); // [GRad/s]
  Rabi_ee = Rabi_gg * m_b[0] / m_a[0];
  Rabi_ge = m_dp[0] * ( conj( ERp*ELm ) + conj( ERm*ELp ) );

}


complex< double > MaxwellBloch::GetERpTerminal( const double time, const int xi )
{
  double retime = time - (xi+Xi+EdgeXi+1)*m_DeltaZ/C;
  double Efield = Intensity2Strength(m_I_ERp); // peak electric field [V/cm]

  if( !m_squarewave ){
    Efield *= exp( - pow((retime-m_D_ERp)/m_T_ERp, 2) / 4. );
  }else if( retime<m_D_ERp || m_D_ERp+m_T_ERp<retime ){
    return complex< double >( 0, 0);
  }

  return complex< double >( Efield, 0);
}

complex< double > MaxwellBloch::GetELmTerminal( const double time, const int xi )
{
  double retime = time + (xi-Xi-EdgeXi-1)*m_DeltaZ/C; 
  double Efield = Intensity2Strength(m_I_ELm); // peak electric field [V/cm]

  if( !m_squarewave ){
    Efield *= exp( - pow((retime-m_D_ELm)/m_T_ELm, 2) /4. );
  }else if( retime<m_D_ELm || m_D_ELm+m_T_ELm<retime ){
    return complex< double >( 0, 0 );
  }
  return complex< double >( Efield, 0);
}

complex< double > MaxwellBloch::GetTrigTerminal( const double time, const int xi )
{
  double retime = time - (xi+Xi+EdgeXi+1)*m_DeltaZ/C;
  double Efield = Intensity2Strength(m_I_Trig); // peak electric field [V/cm]

  if( !m_squarewave ){
    Efield *= exp( - pow((retime-m_D_Trig)/m_T_Trig, 2) /4. );
  }else if( retime<m_D_Trig || m_D_Trig+m_T_Trig<retime ){
    return complex< double >( 0, 0 );
  }
  return complex< double >( Efield, 0);
}


void MaxwellBloch::DumpRefractiveIndex( const double temp /* [K] */ ) const
{
  double density = 273.15 / temp * 2.6867805E25 * 1E-6; // [cm-3]
  cout << "# Density: " << density << " cm-3" << endl;
  cout << "# Temperature: " << temp <<  " K" << endl;

  for( int nm=120; nm<1200; nm++ ){
    double ws = 1E7 / nm; // wavenumber

    double sa=0, sb=0, sd=0;

    for( int j=0; j<Hydrogen::nOmegaB; j++ ){
      sa += Hydrogen::m_uaB[j] * ( 1./(Hydrogen::m_omegaB[j]-ws) + 1./(Hydrogen::m_omegaB[j]+ws) );
      sb += Hydrogen::m_ubB[j] * ( 1./(Hydrogen::m_omegaB[j]-Hydrogen::m_omega_eg_k-ws) + 1./(Hydrogen::m_omegaB[j]-Hydrogen::m_omega_eg_k+ws) );
      sd += sqrt( Hydrogen::m_uaB[j] * Hydrogen::m_ubB[j] ) / ( Hydrogen::m_omegaB[j]-Hydrogen::m_omega_eg_k+ws) * Hydrogen::m_saB[j] * Hydrogen::m_sbB[j];
    }
    for( int j=0; j<Hydrogen::nOmegaC; j++ ){
      sa += Hydrogen::m_uaC[j] * ( 1./(Hydrogen::m_omegaC[j]-ws) + 1./(Hydrogen::m_omegaC[j]+ws) );
      sb += Hydrogen::m_ubC[j] * ( 1./(Hydrogen::m_omegaC[j]-Hydrogen::m_omega_eg_k-ws) + 1./(Hydrogen::m_omegaC[j]-Hydrogen::m_omega_eg_k+ws) );
      sd += sqrt( Hydrogen::m_uaC[j] * Hydrogen::m_ubC[j] ) / ( Hydrogen::m_omegaC[j]-Hydrogen::m_omega_eg_k+ws) * Hydrogen::m_saC[j] * Hydrogen::m_sbC[j];
    }
    double ConversionFactor = 1e4; // [C2 m2 -> C2 cm2]
    ConversionFactor *= 1. / (2*M_PI*C) / (2*HBAR*HBAR); // Kayser -> GRad/s (2 pi c), 2hbar^2 is
    sa *= ConversionFactor * sqrt(3); // [C2 cm2 J-2 ns-1 Rad-1]

    double nres = density * HBAR * sa / EPSILON0;
    cout << nm << "\t" << nres << endl;
  }    

}



int main(int argc, char **argv)
{

  // Arguments check
  /* Arg0 : ./MaxwellBloch
     Arg1 : input setup file name
     Arg2 : output file name
  */
  if( argc != 3 ){
    cerr << "Argument error!" << endl;
    cerr << "Usage: ./MaxwellBloch infile outfile" << endl;
    return -1;
  }
  if( access( argv[1], R_OK )!=0 ){
    cerr << argv[1] << " does not exist." << endl;
    return -1;
  }


  // Setting
  double TargetLength = 0;
  double TargetPressure = 0;
  double Detuning = 0;
  double Gamma_1g = 0;
  double Gamma_1e = 0;
  double Gamma_2 = 0;
  double ERpIntensity = 0;
  double ELmIntensity = 0;
  double ERpWidth = 0;
  double ELmWidth = 0;
  double ERpDelay = 0;
  double ELmDelay = 0;
  double TrigIntensity = 0;
  double TrigWidth = 0;
  double TrigDelay = 0;
  bool Resonance = false;
  bool Absorption = false;
  bool Relaxation = false;
  bool Square = false;
  int EType = 0;
  int SAlg = 0;
  int TAlg = 0;
  int nThreads = 0;

  ifstream ifs( argv[1] );
  string line;
  while( getline( ifs, line ) ){
    istringstream iss( line );
    string header;
    iss >> header;
    if( header == "TargetLength" ){ iss >> TargetLength; }
    else if( header == "Pressure" ){ iss >> TargetPressure; }
    else if( header == "Xi" ){ iss >> Xi; }
    else if( header == "OutputXiPitch" ){ iss >> OutputXiPitch; }
    else if( header == "SimulationTime" ){ iss >> SimTime; }
    else if( header == "SimTimePitch" ){ iss >> SimTimePitch; }
    else if( header == "OutputTimePitch" ){ iss >> OutputTimePitch; }
    else if( header == "Detuning/2pi" ){ iss >> Detuning; Detuning *= 2*M_PI; }
    else if( header == "Gamma1g/2pi" ){ iss >> Gamma_1g; Gamma_1g *= 2*M_PI; }
    else if( header == "Gamma1e/2pi" ){ iss >> Gamma_1e; Gamma_1e *= 2*M_PI; }
    else if( header == "Gamma2/2pi" ){ iss >> Gamma_2; Gamma_2 *= 2*M_PI; }
    else if( header == "ERpIntensity" ){ iss >> ERpIntensity; }
    else if( header == "ELmIntensity" ){ iss >> ELmIntensity; }
    else if( header == "ERpWidth" ){ iss >> ERpWidth; }
    else if( header == "ELmWidth" ){ iss >> ELmWidth; }
    else if( header == "ERpDelay" ){ iss >> ERpDelay; }
    else if( header == "ELmDelay" ){ iss >> ELmDelay; }
    else if( header == "TrigIntensity" ){ iss >> TrigIntensity; }
    else if( header == "TrigWidth" ){ iss >> TrigWidth; }
    else if( header == "TrigDelay" ){ iss >> TrigDelay; }
    else if( header == "ForceResonance" ){ iss >> Resonance; }
    else if( header == "Absorption" ){ iss >> Absorption; }
    else if( header == "Relaxation" ){ iss >> Relaxation; }
    else if( header == "Square" ){ iss >> Square; }
    else if( header == "Equation" ){ iss >> EType; }
    else if( header == "SpatialAlgorithm" ){ iss >> SAlg; }
    else if( header == "TimeAlgorithm" ){ iss >> TAlg; }
    else if( header == "OpenMP" ){ iss >> nThreads; omp_set_num_threads(nThreads); }
  }

  // Open MP
#ifdef _OPENMP
  cerr << "OpenMP is used. " << omp_get_max_threads() << " threads is working." << endl;
#else
  cerr << "OpenMP is not used." << endl;
#endif


  



  MaxwellBloch mb( TargetLength, TargetPressure ); // cm

  mb.SetERpIntensity( ERpIntensity ); // GW/cm2
  mb.SetELmIntensity( ELmIntensity ); // GW/cm2
  mb.SetERpWidth( ERpWidth ); // sigma of intensity [ns]
  mb.SetELmWidth( ELmWidth ); // sigma of intensity [ns]
  mb.SetERpDelay( ERpDelay ); // peak delay [ns]
  mb.SetELmDelay( ELmDelay ); // peak delay [ns]
  //mb.SetSquareWave( Suare );

  mb.SetTrigIntensity( TrigIntensity );
  mb.SetTrigWidth( TrigWidth );
  mb.SetTrigDelay( TrigDelay );

  mb.SetResonance( Resonance );
  mb.SetDelta( Detuning );
  mb.SetAbsorption( Absorption );
  mb.SetRelaxation( Relaxation );

  mb.SetEquationType( EquationType(EType) );
  mb.SetSpatialAlgorithm( SpatialAlgorithm(SAlg) );
  mb.SetTimeAlgorithm( TimeAlgorithm(TAlg) );

  mb.SetGammas( Gamma_1g, Gamma_1e, Gamma_2 );


  ofs.open( argv[2] );
  mb.DumpSetting();


  //mb.DumpRefractiveIndex( 273.15+20.0);



  ofs << std::hexfloat;

  
  // dt is the initial step size. The actual step size is changed according to error control of the stepper. For the last step, the step size will be reduced to ensure we end exactly at t1. The observer is called after each time step. (No dense output)
  //integrate_adaptive( make_dense_output( AbsoluteError, RelativeError, runge_kutta_dopri5< MaxwellBloch::state_type >() ), mb, mb.state, 0.0, SimTime, SimTimePitch, MaxwellBloch::Write );

  switch( mb.GetTimeAlgorithm() ){
  case Dopri5:
    // dt is the initial step size. The actual step size will be adjusted during integration due to error control. Observer output is interpolated value at t=t+n dt.
    integrate_const( make_dense_output( AbsoluteError, RelativeError, runge_kutta_dopri5< MaxwellBloch::state_type >() ), mb, mb.state, 0.0, SimTime, SimTimePitch, MaxwellBloch::Write );
    break;
  case BulirschStoer:
    // bulirsch_stoer usage
    bulirsch_stoer_dense_out< MaxwellBloch::state_type > stepper( AbsoluteError, RelativeError, 1.0 , 1.0 );
    integrate_const( stepper, mb, mb.state, 0.0, SimTime, SimTimePitch, MaxwellBloch::Write );
    break;
  default:
    cerr << "Time algorithm " << mb.GetTimeAlgorithm() << " is not defined." << endl;
    return -1;
  }

  //integrate_const( runge_kutta4< MaxwellBloch::state_type >(), mb, mb.state, 0.0, SimTime, SimTimePitch, MaxwellBloch::Write );

  ofs.close();
  return 0;

}

