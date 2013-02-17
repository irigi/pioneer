// Copyright © 2008-2013 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#include "libs.h"
#include "DynamicBody.h"
#include "Space.h"
#include "Frame.h"
#include "Serializer.h"
#include "Game.h"
#include "Planet.h"
#include "Pi.h"

DynamicBody::DynamicBody(): ModelBody()
{
	m_flags = Body::FLAG_CAN_MOVE_FRAME;
	m_oldPos = GetPosition();
	m_oldAngDisplacement = vector3d(0.0);
	m_force = vector3d(0.0);
	m_torque = vector3d(0.0);
	m_vel = vector3d(0.0);
	m_angVel = vector3d(0.0);
	m_mass = 1;
	m_angInertia = 1;
	m_massRadius = 1;
	m_isMoving = true;
	m_atmosForce = vector3d(0.0);
	m_gravityForce = vector3d(0.0);
	m_externalForce = vector3d(0.0);		// do external forces calc instead?
	m_lastForce = vector3d(0.0);
	m_lastTorque = vector3d(0.0);
}

void DynamicBody::SetForce(const vector3d &f)
{
	m_force = f;
}

void DynamicBody::AddForce(const vector3d &f)
{
	m_force += f;
}

void DynamicBody::AddTorque(const vector3d &t)
{
	m_torque += t;
}

void DynamicBody::AddRelForce(const vector3d &f)
{
	m_force += GetOrient() * f;
}

void DynamicBody::AddRelTorque(const vector3d &t)
{
	m_torque += GetOrient() * t;
}

void DynamicBody::Save(Serializer::Writer &wr, Space *space)
{
	ModelBody::Save(wr, space);
	wr.Vector3d(m_force);
	wr.Vector3d(m_torque);
	wr.Vector3d(m_vel);
	wr.Vector3d(m_angVel);
	wr.Double(m_mass);
	wr.Double(m_massRadius);
	wr.Double(m_angInertia);
	wr.Bool(m_isMoving);
}

void DynamicBody::Load(Serializer::Reader &rd, Space *space)
{
	ModelBody::Load(rd, space);
	m_force = rd.Vector3d();
	m_torque = rd.Vector3d();
	m_vel = rd.Vector3d();
	m_angVel = rd.Vector3d();
	m_mass = rd.Double();
	m_massRadius = rd.Double();
	m_angInertia = rd.Double();
	m_isMoving = rd.Bool();
}

void DynamicBody::PostLoadFixup(Space *space)
{
	Body::PostLoadFixup(space);
	m_oldPos = GetPosition();
//	CalcExternalForce();		// too dangerous
}

void DynamicBody::SetTorque(const vector3d &t)
{
	m_torque = t;
}

void DynamicBody::SetMass(double mass)
{
	m_mass = mass;
	// This is solid sphere mass distribution, my friend
	m_angInertia = (2/5.0)*m_mass*m_massRadius*m_massRadius;
}

void DynamicBody::SetFrame(Frame *f)
{
	ModelBody::SetFrame(f);
	// external forces will be wrong after frame transition
	m_externalForce = m_gravityForce = m_atmosForce = vector3d(0.0);
}

void DynamicBody::CalcExternalForce()
{
	// gravity
	if (!GetFrame()) return;			// no external force if not in a frame
	Body *body = GetFrame()->GetBody();
	if (body && !body->IsType(Object::SPACESTATION)) {	// they ought to have mass though...
		vector3d b1b2 = GetPosition();
		double m1m2 = GetMass() * body->GetMass();
		double invrsqr = 1.0 / b1b2.LengthSqr();
		double force = G*m1m2 * invrsqr;
		m_externalForce = -b1b2 * sqrt(invrsqr) * force;
	}
	else m_externalForce = vector3d(0.0);
	m_gravityForce = m_externalForce;

	// atmospheric drag
	if (GetFrame()->IsRotFrame() && body->IsType(Object::PLANET))
	{
		Planet *planet = static_cast<Planet*>(body);
		double dist = GetPosition().Length();
		double speed = m_vel.Length();
		double pressure, density;
		planet->GetAtmosphericState(dist, &pressure, &density);
		const double radius = GetClipRadius();		// bogus, preserving behaviour
		const double AREA = radius;
		// ^^^ yes that is as stupid as it looks
		const double DRAG_COEFF = 0.1; // 'smooth sphere'
		vector3d dragDir = -m_vel.NormalizedSafe();
		vector3d fDrag = 0.5*density*speed*speed*AREA*DRAG_COEFF*dragDir;

		// make this a bit less daft at high time accel
		// only allow atmosForce to increase by .1g per frame
		vector3d f1g = m_atmosForce + dragDir * GetMass();
		if (fDrag.LengthSqr() > f1g.LengthSqr()) m_atmosForce = f1g;
		else m_atmosForce = fDrag;

		m_externalForce += m_atmosForce;
	}
	else m_atmosForce = vector3d(0.0);

	// centrifugal and coriolis forces for rotating frames
	if (GetFrame()->IsRotFrame()) {
		vector3d angRot(0, GetFrame()->GetAngSpeed(), 0);
		m_externalForce -= m_mass * angRot.Cross(angRot.Cross(GetPosition()));	// centrifugal
		m_externalForce -= 2 * m_mass * angRot.Cross(GetVelocity());			// coriolis
	}
}

void DynamicBody::TimeStepUpdate(const float timeStep)
{
	m_oldPos = GetPosition();
	if (m_isMoving) {
		m_force += m_externalForce;

		m_vel += double(timeStep) * m_force * (1.0 / m_mass);
		m_angVel += double(timeStep) * m_torque * (1.0 / m_angInertia);

		double len = m_angVel.Length();
		if (len > 1e-16) {
			vector3d axis = m_angVel * (1.0 / len);
			matrix3x3d r = matrix3x3d::Rotate(len * timeStep, axis);
			SetOrient(r * GetOrient());
		}
		m_oldAngDisplacement = m_angVel * timeStep;

		SetPosition(GetPosition() + m_vel * double(timeStep));

//if (this->IsType(Object::PLAYER))
//printf("pos = %.1f,%.1f,%.1f, vel = %.1f,%.1f,%.1f, force = %.1f,%.1f,%.1f, external = %.1f,%.1f,%.1f\n",
//	pos.x, pos.y, pos.z, m_vel.x, m_vel.y, m_vel.z, m_force.x, m_force.y, m_force.z,
//	m_externalForce.x, m_externalForce.y, m_externalForce.z);

		m_lastForce = m_force;
		m_lastTorque = m_torque;
		m_force = vector3d(0.0);
		m_torque = vector3d(0.0);
		CalcExternalForce();			// regenerate for new pos/vel
	} else {
		m_oldAngDisplacement = vector3d(0.0);
	}
}

void DynamicBody::UpdateInterpTransform(double alpha)
{
	m_interpPos = alpha*GetPosition() + (1.0-alpha)*m_oldPos;

	double len = m_oldAngDisplacement.Length() * (1.0-alpha);
	if (len > 1e-16) {
		vector3d axis = m_oldAngDisplacement.Normalized();
		matrix3x3d rot = matrix3x3d::Rotate(-len, axis);		// rotate backwards
		m_interpOrient = rot * GetOrient();
	}
	else m_interpOrient = GetOrient();
}

void DynamicBody::SetMassDistributionFromModel()
{
	CollMesh *m = GetCollMesh();
	// XXX totally arbitrarily pick to distribute mass over a half
	// bounding sphere area
	m_massRadius = m->GetRadius()*0.5f;
	SetMass(m_mass);
}

vector3d DynamicBody::GetAngularMomentum() const
{
	return m_angInertia * m_angVel;
}

DynamicBody::~DynamicBody()
{
}

vector3d DynamicBody::GetVelocity() const
{
	return m_vel;
}

void DynamicBody::SetVelocity(const vector3d &v)
{
	m_vel = v;
}

vector3d DynamicBody::GetAngVelocity() const
{
	return m_angVel;
}

void DynamicBody::SetAngVelocity(const vector3d &v)
{
	m_angVel = v;
}

#define KINETIC_ENERGY_MULT	0.00001f
bool DynamicBody::OnCollision(Object *o, Uint32 flags, double relVel)
{
	// don't bother doing collision damage from a missile that will now explode, or may have already
	// also avoids an occasional race condition where destruction event of this could be queued twice
	// returning true to ensure that the missile can react to the collision
	if (o->IsType(Object::MISSILE)) return true;

	double kineticEnergy = 0;
	if (o->IsType(Object::DYNAMICBODY)) {
		kineticEnergy = KINETIC_ENERGY_MULT * static_cast<DynamicBody*>(o)->GetMass() * relVel * relVel;
	} else {
		kineticEnergy = KINETIC_ENERGY_MULT * m_mass * relVel * relVel;
	}
	// damage (kineticEnergy is being passed as a damage value) is measured in kilograms
	// ignore damage less than a gram
	if (kineticEnergy > 1e-3) OnDamage(o, float(kineticEnergy));
	return true;
}

// return parameters for orbit of any body, gives both elliptic and hyperbolic trajectories
Orbit *DynamicBody::ReturnOrbit() {
	Orbit *ret = &(this->orbit);

	Frame *fram = this->GetFrame();
	if(fram->IsRotFrame()) fram = fram->GetNonRotFrame();
	double mass = fram->GetSystemBody()->GetMass();

	// current velocity and position with respect to non-rotating frame
	vector3d vel = this->GetVelocityRelTo(fram);
	vector3d pos = this->GetPositionRelTo(fram);
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

		// correct sign of offset is given by sign pos.Dot(vel) (heading towards apogeum or porigeum?]
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

	ret->orbitalPhaseAtStart = -offset; // in time t = 0

	double M_t0; // mean anomaly at t = 0
	const double e = ret->eccentricity;
	if(e > 0 && e < 1) {
		M_t0 = 2*atan(tan(ret->orbitalPhaseAtStart/2)*sqrt((1-e)/(1+e)));
		M_t0 = M_t0 - e*sin(M_t0);
		M_t0 -= Pi::game->GetTime() * ret->velocityAreaPerSecond / ret->semiMajorAxis / ret->semiMajorAxis / sqrt(1 - e*e);
	} else {
		M_t0 = 2*atanh(tan(ret->orbitalPhaseAtStart/2)*sqrt((e-1)/(1+e)));
		M_t0 = M_t0 - e*sinh(M_t0);
		M_t0 -= Pi::game->GetTime()  * ret->velocityAreaPerSecond / ret->semiMajorAxis / ret->semiMajorAxis / sqrt(e*e-1);
	}

	ret->orbitalPhaseAtStart = M_t0;

	return ret;
}
