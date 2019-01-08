#ifndef QUATERNION_BASE_ADDONS_H
#define QUATERNION_BASE_ADDONS_H

#ifndef EIGEN_QUATERNION_PLUGIN
	#define EIGEN_QUATERNIONBASE_PLUGIN "math/eigen_QuaternionBaseAddons.h"
#endif

#else // QUATERNION_BASE_ADDONS_H

Derived& pure(const Vector3 &v) 
{
	w() = 0.0;
	x() = v.x();
	y() = v.y();
	z() = v.z();
	return derived();
}

Vector3 rotate(const Vector3 &v) const
{ 
	double q0 = w();
	Vector3 qv = vec();
	Vector3 Lv = (q0*q0 - qv.squaredNorm())*v + 2 * qv.dot(v)*qv + 2 * q0*qv.cross(v);

	return Lv;
}

Vector3 rotate2(const Vector3 &v) const { return conjugate().rotate(v); }

#endif // QUATERNION_BASE_ADDONS_H
