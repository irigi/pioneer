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

Orbit * Orbit::calc_orbit(Orbit * ret, vector3d pos, vector3d vel, double mass)  {
	double r_now = pos.Length() + 1e-12;
	double v_now = vel.Length() + 1e-12;

	// angular momentum
	vector3d ang = -(vel.Cross(pos));
	double LL =ang.Length() + 1e-12;

	// total energy
	double EE = vel.LengthSqr()/2 - mass*6.672e-11/r_now;

	if(mass <= 1e-3 || r_now <= 1e-3 || v_now <= 1e-3 || fabs(EE) <= 1e-12 || 1.0-ang.z*ang.z/LL/LL < 0) {
		ret->eccentricity = 0;
		ret->semiMajorAxis = 0;
		ret->velocityAreaPerSecond = 0;
		ret->orbitalPhaseAtStart = 0;
		ret->rotMatrix =  matrix3x3d::RotateY(0);
		return ret;
	}

	// http://en.wikipedia.org/wiki/Orbital_eccentricity
	ret->eccentricity = 1 + 2*EE*LL*LL/pow(mass*6.672e-11, 2);
	if(ret->eccentricity < 1e-12) ret->eccentricity = 1e-12;
	if(ret->eccentricity == 1.0) ret->eccentricity = 1-1e-6;
	ret->eccentricity = sqrt(ret->eccentricity);

	// lines represent these quantities:
	// 		(e M G)^2
	// 		M G (e - 1) / 2 EE, always positive (EE and (e-1) change sign
	// 		M G / 2 EE,
	// which is a (http://en.wikipedia.org/wiki/Semi-major_axis); a of hyperbola is taken as positive here
	ret->semiMajorAxis = 2*EE*LL*LL + pow(mass*6.672e-11, 2);
	if(ret->semiMajorAxis < 0) ret->semiMajorAxis  = 0;
	ret->semiMajorAxis = (sqrt(ret->semiMajorAxis ) - mass*6.672e-11)/2/EE;
	ret->semiMajorAxis = ret->semiMajorAxis/fabs(1.0-ret->eccentricity);

	// The formulas for rotation matrix were derived based on following assumptions:
	//	1. Trajectory follows Kepler's law and vector {-r cos(v), -r sin(v), 0}, r(t) and v(t) are parameters.
	//	2. Correct transformation must transform {0,0,LL} to ang and {-r_now cos(orbitalPhaseAtStart), -r_now sin(orbitalPhaseAtStart), 0} to pos.
	//  3. orbitalPhaseAtStart (=offset) is calculated from r = a ((e^2 - 1)/(1 + e cos(v) ))
	double angle1 = acos(Clamp(ang.z/LL ,-1 + 1e-6,1 - 1e-6)) * (ang.x > 0 ? -1 : 1),
			angle2 = asin(Clamp(ang.y / LL / sqrt(1.0 - ang.z*ang.z/LL/LL) ,-1 + 1e-6,1 - 1e-6) ) * (ang.x > 0 ? -1 : 1);

	// There are two possible solutions of the equation and the only way how to find the correct one
	// I know about is to try both and check if the position is transformed correctly. We minimize the difference
	// of the transformed  position and expected result.
	double value = 1e99, offset = 0, cc = 0;
	for(int i = -1; i <= 1; i += 2) {
		double off = 0, ccc = 0;
		matrix3x3d mat;

		if(ret->eccentricity < 1) {
			off = ret->semiMajorAxis*(1 - ret->eccentricity*ret->eccentricity) - r_now;
		} else {
			off = ret->semiMajorAxis*(-1 + ret->eccentricity*ret->eccentricity) - r_now;
		}

		// correct sign of offset is given by sign pos.Dot(vel) (heading towards apohelion or perihelion?]
		off = Clamp(off/(r_now * ret->eccentricity), -1 + 1e-6,1 - 1e-6);
		off = -pos.Dot(vel)/fabs(pos.Dot(vel))*acos(off);

		ccc = acos(-pos.z/r_now/sin(angle1)) * i;
		mat = matrix3x3d::RotateZ(angle2) * matrix3x3d::RotateY(angle1) * matrix3x3d::RotateZ(ccc - off);

		if(((mat*vector3d(-r_now*cos(off),r_now*sin(off),0)) - pos).Length() < value) {
			value = ((mat*vector3d(-r_now*cos(off),r_now*sin(off),0)) - pos).Length();
			cc = ccc;
			offset = off;
		}
	}

	// matrix3x3d::RotateX(M_PI) and minus sign before offset changes solution above, derived for orbits {-r cos(v), -r sin(v), 0}
	// to {-r cos(v), -r sin(v), 0}
	ret->rotMatrix = matrix3x3d::RotateZ(angle2) * matrix3x3d::RotateY(angle1) * matrix3x3d::RotateZ(cc - offset) * matrix3x3d::RotateX(M_PI);
	ret->velocityAreaPerSecond = Orbit::calc_velocity_area_per_sec(ret->semiMajorAxis, mass,ret->eccentricity);

	ret->orbitalPhaseAtStart = ret->MeanAnomalyFromTrueAnomaly(-offset);

	return ret;
}


