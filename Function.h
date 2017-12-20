#ifndef _FUNCTION_H

#include "DataFile.h"

class Function {
  private:
    int _initial_condition_choice; // Gaussian or Rectangular
    enum {Gaussian, Rectangular};
    int _velocity_choice; // Uniform, Rotational or Sinusoidal
    enum {Uniform, Rotational, Sinusoidal};

    // parameters x0, y0, a, b, c, d, e and f
    double _x0, _y0, _a, _b, _c, _d, _e, _f;

    bool _in_flow;
	public: // Méthodes et opérateurs de la classe
    Function(DataFile* data_file);
    double ExactSolution(const double x, const double y, const double t) const;
    double InitialCondition(const double x, const double y) const;
    double Velocity_x(const double x, const double y, const double t) const;
    double Velocity_y(const double x, const double y, const double t) const;
};

#define _FUNCTION_H
#endif
