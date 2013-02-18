// Copyright © 2008-2013 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#include "Serializer.h"
#include "Pi.h"
#include <map>
#include "utils.h"
#include "Game.h"
#include "Orbit.h"

double Orbit::MeanAnomalyAtTime(double time) const {
	const double e = eccentricity;
	if(e >= 0 && e < 1) { // elliptic orbit
		return 2.0*M_PI*time / Period() + orbitalPhaseAtStart; // mean anomaly
	} else {
		return -2 * time * velocityAreaPerSecond / semiMajorAxis / semiMajorAxis / sqrt(e*e-1) + orbitalPhaseAtStart; // mean anomaly
	}
}

vector3d Orbit::OrbitalPosAtTime(double t) const
{
	const double e = eccentricity;
	double r = 0, cos_v = 0, sin_v = 0;
	const double M = MeanAnomalyAtTime(t); // mean anomaly

	if(e >= 0 && e < 1) { // elliptic orbit
		// eccentric anomaly
		// NR method to solve for E: M = E-sin(E)
		double E = M;
		for (int iter=5; iter > 0; --iter) {
			E = E - (E-e*(sin(E))-M) / (1.0 - e*cos(E));
		}
		// heliocentric distance
		r = semiMajorAxis * (1.0 - e*cos(E));
		// true anomaly (angle of orbit position)
		cos_v = (cos(E) - e) / (1.0 - e*cos(E));
		sin_v = (sqrt(1.0-e*e)*sin(E))/ (1.0 - e*cos(E));

	} else if(e > 1) { // hyperbolic orbit
		// eccentric anomaly
		// NR method to solve for E: M = E-sinh(E)
		// sinh E and cosh E are solved directly, because of inherent numerical instability of tanh(k arctanh x)
		double ch, sh = 2;
		for (int iter=5; iter > 0; --iter) {
			sh = sh - (M + e*sh - asinh(sh))/(e - 1/sqrt(1 + pow(sh,2)));
		}

		ch = sqrt(1 + sh*sh);

		// heliocentric distance
		r = semiMajorAxis * (e*ch - 1.0);

		// true anomaly (angle of orbit position)
		cos_v = (ch - e) / (1.0 - e*ch);
		sin_v = (sqrt(e*e-1.0)*sh)/ (e*ch - 1.0);
	}

	vector3d pos = vector3d(-cos_v*r, sin_v*r, 0);
	pos = rotMatrix * pos;
	return pos;
}

// used for stepping through the orbit in small fractions
// mean anomaly <-> true anomaly conversion doesn't have
// to be taken into account
vector3d Orbit::EvenSpacedPosTrajectory(double angle) const
{
	const double e = eccentricity;
	vector3d pos = vector3d(0.0f,0.0f,0.0f);

	if(e < 1) {
		double v = 2*M_PI*angle +TrueAnomaly(MeanAnomalyAtTime(Pi::game->GetTime()));
		const double r = semiMajorAxis * (1 - e*e) / (1 + e*cos(v));
		pos = vector3d(-cos(v)*r, sin(v)*r, 0);
	} else {
		double v = 2*M_PI*angle +TrueAnomaly(MeanAnomalyAtTime(Pi::game->GetTime()));
		double r = semiMajorAxis * (e*e - 1) / (1 + e*cos(v));

		// planet is in infinity
		if(v <= - acos(-1/e)) {
			v = - acos(-1/e) + 0.0001;
			r =  1.5e13; // 100 AU
		}
		if(v >= acos(-1/e)) {
			v = acos(-1/e) - 0.0001;
			r =  1.5e13; // 100 AU
		}

		pos = vector3d(-cos(v)*r, sin(v)*r, 0);
	}
	pos = rotMatrix * pos;
	return pos;
}

double Orbit::Period() const {
	if(eccentricity < 1 && eccentricity >= 0) {
		return M_PI * semiMajorAxis * semiMajorAxis * sqrt(1 - eccentricity * eccentricity)/ velocityAreaPerSecond;
	} else { // hyperbola.. period makes no sense, should not be used
		assert(0);
		return 0;
	}
}

double Orbit::TrueAnomaly(double MeanAnomaly) const {
	const double e = eccentricity, M = MeanAnomaly;
	double cos_v, sin_v, v;
	if(e >= 0 && e < 1) { // elliptic orbit
		// eccentric anomaly
		// NR method to solve for E: M = E-sin(E)
		double E = M;
		for (int iter=5; iter > 0; --iter) {
			E = E - (E-e*(sin(E))-M) / (1.0 - e*cos(E));
		}

		// true anomaly (angle of orbit position)
		cos_v = (cos(E) - e) / (1.0 - e*cos(E));
		sin_v = (sqrt(1.0-e*e)*sin(E))/ (1.0 - e*cos(E));

	} else if(e > 1) { // hyperbolic orbit
		// eccentric anomaly
		// NR method to solve for E: M = E-sinh(E)
		// sinh E and cosh E are solved directly, because of inherent numerical instability of tanh(k arctanh x)
		double ch, sh = 2;
		for (int iter=5; iter > 0; --iter) {
			sh = sh - (M + e*sh - asinh(sh))/(e - 1/sqrt(1 + pow(sh,2)));
		}

		ch = sqrt(1 + sh*sh);

		// true anomaly (angle of orbit position)
		cos_v = (ch - e) / (1.0 - e*ch);
		sin_v = (sqrt(e*e-1.0)*sh)/ (e*ch - 1.0);
	}

	v = atan2(sin_v, cos_v);
	return v;
}

double Orbit::MeanAnomalyFromTrueAnomaly(double trueAnomaly) const {
	double M_t0;
	const double e = eccentricity;
	if(e > 0 && e < 1) {
		M_t0 = 2*atan(tan(trueAnomaly/2)*sqrt((1-e)/(1+e)));
		M_t0 = M_t0 - e*sin(M_t0);
		M_t0 -= Pi::game->GetTime() *2.0*M_PI/ Period();
	} else {
		// For hyperbolic trajectories, mean anomaly has opposite sign to true anomaly, therefore trajectories which go forward
		// in time decrease their true anomaly. Yes, it is confusing.
		M_t0 = 2*atanh(tan(trueAnomaly/2)*sqrt((e-1)/(1+e)));
		M_t0 = M_t0 - e*sinh(M_t0);
		M_t0 += Pi::game->GetTime()  * 2 * velocityAreaPerSecond / semiMajorAxis / semiMajorAxis / sqrt(e*e-1);
	}

	return M_t0;
}

vector3d Orbit::Apogeum() const {
	if(eccentricity < 1) {
		return semiMajorAxis * (1 + eccentricity) * (rotMatrix * vector3d(1,0,0));
	} else {
		return vector3d(0,0,0);
	}
}

vector3d Orbit::Perigeum() const {
	if(eccentricity < 1) {
		return semiMajorAxis * (1 - eccentricity) * (rotMatrix * vector3d(-1,0,0));
	} else {
		return semiMajorAxis * (eccentricity - 1) * (rotMatrix * vector3d(-1,0,0));
	}
}

double Orbit::calc_velocity_area_per_sec(double semiMajorAxis, double centralMass, double eccentricity) {
	if(eccentricity < 1)
		return M_PI * semiMajorAxis * semiMajorAxis * sqrt(1 - eccentricity * eccentricity)/ calc_orbital_period(semiMajorAxis, centralMass);
	else
		return M_PI * semiMajorAxis * semiMajorAxis * sqrt(eccentricity * eccentricity - 1)/ calc_orbital_period(semiMajorAxis, centralMass);
}

double Orbit::calc_velocity_area_per_sec_gravpoint(double semiMajorAxis, double totalMass, double bodyMass, double eccentricity) {
	if(eccentricity < 1) {
		return M_PI * semiMajorAxis * semiMajorAxis * sqrt(1 - eccentricity * eccentricity)/ calc_orbital_period_gravpoint(semiMajorAxis, totalMass, bodyMass);
	} else {
		return M_PI * semiMajorAxis * semiMajorAxis * sqrt(eccentricity * eccentricity - 1)/ calc_orbital_period_gravpoint(semiMajorAxis, totalMass, bodyMass);
	}
}

double Orbit::calc_orbital_period(double semiMajorAxis, double centralMass)
{
	return 2.0*M_PI*sqrt((semiMajorAxis*semiMajorAxis*semiMajorAxis)/(G*centralMass));
}

double Orbit::calc_orbital_period_gravpoint(double semiMajorAxis, double totalMass, double bodyMass)
{
	// variable names according to the formula in:
	// http://en.wikipedia.org/wiki/Barycentric_coordinates_(astronomy)#Two-body_problem
	//
	// We have a 2-body orbital system, represented as a gravpoint (at the barycentre),
	// plus two bodies, each orbiting that gravpoint.
	// We need to compute the orbital period, given the semi-major axis of one body's orbit
	// around the gravpoint, the total mass of the system, and the mass of the body.
	//
	// According to Kepler, the orbital period P is defined by:
	//
	// P = 2*pi * sqrt( a**3 / G*(M1 + M2) )
	//
	// where a is the semi-major axis of the orbit, M1 is the mass of the primary and M2 is
	// the mass of the secondary. But we don't have that semi-major axis value, we have the
	// the semi-major axis for the orbit of the secondary around the gravpoint, instead.
	//
	// So, this first computes the semi-major axis of the secondary's orbit around the primary,
	// and then uses the above formula to compute the orbital period.
	const double r1 = semiMajorAxis;
	const double m2 = (totalMass - bodyMass);
	const double a = r1 * totalMass / m2;
	const double a3 = a*a*a;
	return 2.0 * M_PI * sqrt(a3 / (G * totalMass));
}
