#ifndef UNITS_H
#define UNITS_H

#include "math/math_base.h"

typedef double si_type;

// Unit conversion
constexpr si_type m = 1.0;
constexpr si_type cm = 1e-2 * m;
constexpr si_type mm = 1e-3 * m;
constexpr si_type um = 1e-6 * m;
constexpr si_type nm = 1e-9 * m;
constexpr si_type pm = 1e-12 * m;
constexpr si_type fm = 1e-15 * m;
constexpr si_type km = 1e3 * m;
constexpr si_type inch = 2.54 * cm;
constexpr si_type ft = 12 * inch;
constexpr si_type yd = 3 * ft;
constexpr si_type mi = 5280 * ft;
constexpr si_type nmi = 1852 * m;
constexpr si_type ly = 9460730472580800 * m;
constexpr si_type au = 149597870700 * m;
constexpr si_type pc = 648000 / PI * au;
constexpr si_type kpc = 1e3 * pc;
constexpr si_type mpc = 1e6 * pc;
constexpr si_type gpc = 1e9 * pc;

constexpr si_type rad = 1.0;
constexpr si_type deg = PI / 180 * rad;
constexpr si_type grad = 0.9 * deg;
constexpr si_type arcmin = deg / 60;
constexpr si_type arcsec = arcmin / 60;

constexpr si_type s = 1.0;
constexpr si_type ms = 1e-3 * s;
constexpr si_type us = 1e-6 * s;
constexpr si_type ns = 1e-9 * s;
constexpr si_type ps = 1e-12 * s;
constexpr si_type fs = 1e-15 * s;
constexpr si_type min = 60 * s;
constexpr si_type hr = 60 * min;
constexpr si_type day = 24 * hr;
constexpr si_type yr = 365.25 * day;

constexpr si_type N = 1.0;

constexpr si_type kg = 1.0;

constexpr si_type J = 1.0;
constexpr si_type eV = 1.602176634e-19 * J;
constexpr si_type KeV = 1e3 * eV;
constexpr si_type MeV = 1e6 * eV;
constexpr si_type GeV = 1e9 * eV;
constexpr si_type erg = 1e-7 * J;
constexpr si_type cal = 4.184 * J;

constexpr si_type W = 1.0;

constexpr si_type Pa = 1.0;

constexpr si_type C = 1.0;
constexpr si_type e = 1.602176634e-19 * C;

constexpr si_type V = 1.0;

constexpr si_type A = 1.0;

constexpr si_type K = 1.0;

constexpr si_type becquerel = 1.0 / s;
constexpr si_type Bq = becquerel;
constexpr si_type curie = 3.7e10 * becquerel;
constexpr si_type Ci = curie;
constexpr si_type mCi = 1e-3 * Ci;
constexpr si_type uCi = 1e-6 * Ci;
constexpr si_type gray = 1.0 * J / kg;
constexpr si_type Gy = gray;
constexpr si_type sievert = 1.0 * J / kg;
constexpr si_type Sv = sievert;

constexpr si_type mol = 1.0;

// constants
constexpr si_type c = 299792458 * m / s;
constexpr si_type G = 6.67430e-11 * m * m * m / (kg * s * s);
constexpr si_type h = 6.62607015e-34 * J * s;
constexpr si_type hbar = h / (2 * PI);
constexpr si_type mu0 = 4e-7 * PI * N / A / A;
constexpr si_type k = 1.380649e-23 * J / K;
constexpr si_type e0 = 1 / (mu0 * c * c);
constexpr si_type R = 8.31446261815324 * J / K / mol;
constexpr si_type NA = 6.02214076e23 / mol;
constexpr si_type me = 9.1093837015e-31 * kg; // electron mass
constexpr si_type mp = 1.67262192369e-27 * kg; // proton mass
constexpr si_type mn = 1.67492749804e-27 * kg; // neutron mass

#endif