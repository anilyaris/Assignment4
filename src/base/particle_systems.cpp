#include "particle_systems.hpp"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <numeric>
#include <cstdlib>
#include <ctime>

using namespace std;
using namespace FW;

namespace {

	inline Vec3f fGravity(float mass) {
		return Vec3f(0, -9.8f * mass, 0);
	}

	// force acting on particle at pos1 due to spring attached to pos2 at the other end
	inline Vec3f fSpring(const Vec3f& pos1, const Vec3f& pos2, float k, float rest_length) {
		// YOUR CODE HERE (R2)
		Vec3f disp = pos1 - pos2;
		float d = disp.length();
		return (k * (rest_length - d) / d) * disp;
	}

	inline Vec3f fDrag(const Vec3f& v, float k) {
		// YOUR CODE HERE (R2)
		return -k * v;
	}

	inline Vec3f fWind (float max) {
		srand(time(NULL));
		int mod = 200 * max + 1;
		float x = (rand() % mod) / 100.0f - max;
		float y = (rand() % mod) / 100.0f - max;
		float z = (rand() % mod) / 100.0f - max;
		return Vec3f(x, y, z);
	} 

} // namespace

void SimpleSystem::reset() {
	state_ = State(1, Vec3f(0, radius_, 0));
}

State SimpleSystem::evalF(const State& state) const {
	State f(1, Vec3f(-state[0].y, state[0].x, 0));
	return f;
}

#ifdef EIGEN_SPARSECORE_MODULE_H
// using the implicit Euler method, the simple system should converge towards origin -- as opposed to the explicit Euler, which diverges outwards from the origin.
void SimpleSystem::evalJ(const State&, SparseMatrix& result, bool initial) const {
	if (initial) {
		result.coeffRef(1, 0) = 1.0f;
		result.coeffRef(0, 1) = -1.0f;
	}
}
#endif

Points SimpleSystem::getPoints() {
	return Points(1, state_[0]);
}

Lines SimpleSystem::getLines() {
	static const auto n_lines = 50u;
	auto l = Lines(n_lines * 2);
	const auto angle_incr = 2 * FW_PI / n_lines;
	for (auto i = 0u; i < n_lines; ++i) {
		l[2 * i] = l[2 * i + 1] =
			Vec3f(radius_ * FW::sin(angle_incr * i), radius_ * FW::cos(angle_incr * i), 0);
	}
	rotate(l.begin(), l.begin() + 1, l.end());
	return l;
}

void SpringSystem::reset() {
	const auto start_pos = Vec3f(0.1f, -0.5f, 0.0f);
	const auto spring_k = 30.0f;
	state_ = State(4);
	// YOUR CODE HERE (R2)
	// Set the initial state for a particle system with one particle fixed
	// at origin and another particle hanging off the first one with a spring.
	// Place the second particle initially at start_pos.
	state_[0] = Vec3f();
	state_[1] = Vec3f();
	state_[2] = start_pos;
	state_[3] = Vec3f();

	spring_ = Spring(2, 0, spring_k, start_pos.length());
}

State SpringSystem::evalF(const State& state) const {
	const auto drag_k = 0.5f;
	const auto mass = 1.0f;
	State f(4);
	// YOUR CODE HERE (R2)
	// Return a derivative for the system as if it was in state "state".
	// You can use the fGravity, fDrag and fSpring helper functions for the forces.
	f[0] = Vec3f();
	f[1] = Vec3f();
	f[2] = state_[3];
	f[3] = (fGravity(mass) + fDrag(state_[3], drag_k) + fSpring(state_[spring_.i1], state_[spring_.i2], spring_.k, spring_.rlen)) / mass;
	return f;
}

#ifdef EIGEN_SPARSECORE_MODULE_H

// This is a very useful read for the Jacobians of the spring forces. It deals with spring damping as well, we don't do that -- our drag is simply a linear damping of velocity (that results in some constants in the Jacobian).
// http://blog.mmacklin.com/2012/05/04/implicitsprings/

void SpringSystem::evalJ(const State& state, SparseMatrix& result, bool initial) const {
	const auto drag_k = 0.5f;
	const auto mass = 1.0f;
	// EXTRA: Evaluate the Jacobian into the 'result' matrix here. Only the free end of the spring should have any nonzero values related to it.
}
#endif

Points SpringSystem::getPoints() {
	auto p = Points(2);
	p[0] = state_[0]; p[1] = state_[2];
	return p;
}

Lines SpringSystem::getLines() {
	auto l = Lines(2);
	l[0] = state_[0]; l[1] = state_[2];
	return l;
}

void PendulumSystem::reset() {
	const auto spring_k = 1000.0f;
	const auto start_point = Vec3f(0);
	const auto end_point = Vec3f(0.05, -1.5, 0);
	const auto disp = (end_point - start_point) / (n_ - 1);
	state_ = State(2 * n_);
	// YOUR CODE HERE (R4)
	// Set the initial state for a pendulum system with n_ particles
	// connected with springs into a chain from start_point to end_point.
	for (int i = 0; i < n_; i++) {
		state_[2 * i] = start_point + ((float) i) * disp;
		state_[2 * i + 1] = Vec3f();
		if (i != n_ - 1) springs_.push_back(Spring(i + 1, i, spring_k, disp.length()));
	}	
}

State PendulumSystem::evalF(const State& state) const {
	const auto drag_k = 0.5f;
	const auto mass = 0.5f;
	auto f = State(2 * n_);
	// YOUR CODE HERE (R4)
	// As in R2, return a derivative of the system state "state".
	f[0] = Vec3f();
	f[1] = Vec3f();
	for (int i = 1; i < n_; i++) {
		f[2 * i] = state_[2 * i + 1];
		f[2 * i + 1] = (fGravity(mass) + fDrag(state_[2 * i + 1], drag_k) + fSpring(state_[2 * springs_[i - 1].i1], state_[2 * springs_[i - 1].i2], springs_[i - 1].k, springs_[i - 1].rlen)) / mass;
		if (i != n_ - 1) f[2 * i + 1] -= fSpring(state_[2 * springs_[i].i1], state_[2 * springs_[i].i2], springs_[i].k, springs_[i].rlen);
	}
	return f;
}

#ifdef EIGEN_SPARSECORE_MODULE_H

void PendulumSystem::evalJ(const State& state, SparseMatrix& result, bool initial) const {

	const auto drag_k = 0.5f;
	const auto mass = 0.5f;

	// EXTRA: Evaluate the Jacobian here. Each spring has an effect on four blocks of the matrix -- both of the positions of the endpoints will have an effect on both of the velocities of the endpoints.
}
#endif


Points PendulumSystem::getPoints() {
	auto p = Points(n_);
	for (auto i = 0u; i < n_; ++i) {
		p[i] = state_[i * 2];
	}
	return p;
}

Lines PendulumSystem::getLines() {
	auto l = Lines();
	for (const auto& s : springs_) {
		l.push_back(state_[2 * s.i1]);
		l.push_back(state_[2 * s.i2]);
	}
	return l;
}

void ClothSystem::reset() {
	const auto spring_k = 300.0f;
	const auto width = 1.5f, height = 1.5f; // width and height of the whole grid
	auto w = width / x_, h = height / y_;
	Vec3f origin(-width / 2, 0, 0);
	state_ = State(2 * x_*y_);
	// YOUR CODE HERE (R5)
	// Construct a particle system with a x_ * y_ grid of particles,
	// connected with a variety of springs as described in the handout:
	// structural springs, shear springs and flex springs.
	for (int i = 0; i < x_; i++) {
		for (int j = 0; j < y_; j++) {
			int index = i * y_ + j;
			state_[2 * index] = origin + Vec3f(i * w, 0, -j * h);
			state_[2 * index + 1] = Vec3f();
			if (i < x_ - 1) {
				springs_.push_back(Spring(index + y_, index, spring_k, w));
				if (j > 0) springs_.push_back(Spring(index + y_ - 1, index, spring_k, Vec2f(w, h).length()));
				if (j < y_ - 1) springs_.push_back(Spring(index + y_ + 1, index, spring_k, Vec2f(w, h).length()));
				if (i != x_ - 2) springs_.push_back(Spring(index + 2 * y_, index, spring_k, 2 * w));
			}
			if (j < y_ - 1) {
				springs_.push_back(Spring(index + 1, index, spring_k, h));
				if (j != y_ - 2) springs_.push_back(Spring(index + 2, index, spring_k, 2 * h));
			}
		}		
	}

}

State ClothSystem::evalF(const State& state) const {
	const auto drag_k = 0.08f;
	const auto n = x_ * y_;
	static const auto mass = 0.025f;
	auto f = State(2 * n);
	// YOUR CODE HERE (R5)
	// This will be much like in R2 and R4.
	for (int i = 0; i < n; i++) {
		f[2 * i] = state_[2 * i + 1];
		f[2 * i + 1] = (fGravity(mass) + fDrag(state_[2 * i + 1], drag_k) + fWind(0.7f)) / mass;		
	}

	for (auto s : springs_) {
		auto a = fSpring(state_[2 * s.i1], state_[2 * s.i2], s.k, s.rlen) / mass;
		f[2 * s.i1 + 1] += a;
		f[2 * s.i2 + 1] -= a;
	}
	f[0] = Vec3f();
	f[1] = Vec3f();
	f[2 * (x_ - 1) * y_] = Vec3f();
	f[2 * (x_ - 1) * y_ + 1] = Vec3f();
	return f;
}

#ifdef EIGEN_SPARSECORE_MODULE_H

void ClothSystem::evalJ(const State& state, SparseMatrix& result, bool initial) const {
	const auto drag_k = 0.08f;
	static const auto mass = 0.025f;

	// EXTRA: Evaluate the Jacobian here. The code is more or less the same as for the pendulum.
}

#endif

Points ClothSystem::getPoints() {
	auto n = x_ * y_;
	auto p = Points(n);
	for (auto i = 0u; i < n; ++i) {
		p[i] = state_[2 * i];
	}
	return p;
}

Lines ClothSystem::getLines() {
	float t = 0.05f;
	auto l = Lines();
	for (int s = 0; s < springs_.size(); s++) {
		Spring& spring = springs_[s];
		if (tear && spring.i1 != 0 && spring.i2 != 0 && spring.i1 != (x_ - 1) * y_ && spring.i2 != (x_ - 1) * y_ && (state_[2 * spring.i1] - state_[2 * spring.i2]).length() - spring.rlen > t) {
			springs_.erase(springs_.begin() + s);
		}
		else {
			l.push_back(state_[2 * spring.i1]);
			l.push_back(state_[2 * spring.i2]);
		}
	}
	return l;
}

void GeneratorSystem::reset() {
	state_.clear();
	ages_.clear();
}

State GeneratorSystem::evalF(const State& state) const {
	static const auto mass = 0.025f;
	auto n = ages_.size();
	auto f = State(2 * n);
	// YOUR CODE HERE (R5)
	// This will be much like in R2 and R4.
	for (int i = 0; i < n; i++) {
		if (ages_[i] != t_ - 1) {
			f[2 * i] = state_[2 * i + 1];
			f[2 * i + 1] = (fGravity(mass) + fWind(1.0f)) / mass;			
		}		
	}

	return f;
}

#ifdef EIGEN_SPARSECORE_MODULE_H

void GeneratorSystem::evalJ(const State& state, SparseMatrix& result, bool initial) const {
	const auto drag_k = 0.08f;
	static const auto mass = 0.025f;

	// EXTRA: Evaluate the Jacobian here. The code is more or less the same as for the pendulum.
}

#endif

Points GeneratorSystem::getPoints() {
	for (int i = 0; i < ages_.size(); i++) ages_[i]++;
	while (ages_.size() > 0 && ages_[0] == t_) {
		ages_.erase(ages_.begin());
		state_.erase(state_.begin(), state_.begin() + 2);
	}

	float bound = 1.0f;
	for (int i = 0; i < state_.size(); i = i + 2) {
		Vec3f pos = state_[i], vel = state_[i + 1];
		if (pos.x <= -bound) {
			pos.x = -bound;
			vel.x = -vel.x;
		}
		if (pos.x >= bound) {
			pos.x = bound;
			vel.x = -vel.x;
		}
		if (pos.y <= -bound) {
			pos.y = -bound;
			vel.y = -vel.y;
		}
		if (pos.y >= bound) {
			pos.y = bound;
			vel.y = -vel.y;
		}
		if (pos.z <= -bound) {
			pos.z = -bound;
			vel.z = -vel.z;
		}
		if (pos.z >= bound) {
			pos.z = bound;
			vel.z = -vel.z;
		}
		state_[i] = pos;
		state_[i + 1] = vel;
	}

	const float vlim = 30.0f;
	Vec3f origin(0, 0, 0);
	unsigned numNew = rand() % n_ + 1;
	for (int i = 0; i < numNew; i++) {
		state_.push_back(origin);
		state_.push_back(fWind(vlim));
		ages_.push_back(0);
	}

	auto n = ages_.size();
	auto p = Points(n);
	for (auto i = 0u; i < n; ++i) {
		p[i] = state_[2 * i];
	}
	return p;
}

Lines GeneratorSystem::getLines() {
	return Lines();	
}

