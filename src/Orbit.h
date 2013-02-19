/*
 * Orbit.h
 *
 *  Created on: Feb 18, 2013
 *      Author: olsij4am
 */

#ifndef ORBIT_H_
#define ORBIT_H_

struct Orbit {
	Orbit(): orbitalPhaseAtStart(0.0) {};
	vector3d OrbitalPosAtTime(double t) const;
	// 0.0 <= t <= 1.0. Not for finding orbital pos
	vector3d EvenSpacedPosTrajectory(double angle) const;
	/* duplicated from SystemBody... should remove probably */
	static double calc_orbital_period(double semiMajorAxis, double centralMass);
	static double calc_orbital_period_gravpoint(double semiMajorAxis, double totalMass, double bodyMass);
	static double calc_velocity_area_per_sec(double semiMajorAxis, double centralMass, double eccentricity);
	static double calc_velocity_area_per_sec_gravpoint(double semiMajorAxis, double totalMass, double bodyMass, double eccentricity);

	static Orbit *calc_orbit(Orbit * ret, vector3d current_pos, vector3d current_vel, double mass);

	double Period() const;
	double TrueAnomaly(double MeanAnomaly) const;
	double MeanAnomalyFromTrueAnomaly(double trueAnomaly) const;
	double MeanAnomalyAtTime(double time) const;
	vector3d Apogeum() const;
	vector3d Perigeum() const;
	double eccentricity;
	double semiMajorAxis;
	double orbitalPhaseAtStart; // 0 to 2 pi radians
	/* dup " " --------------------------------------- */
	double velocityAreaPerSecond; // seconds
	matrix3x3d rotMatrix;
};


#endif /* ORBIT_H_ */
