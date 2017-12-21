#ifndef _FUNCTION_CPP

#include "Function.h"
#include <cmath>

Function::Function(DataFile* data_file) :
_x0(data_file->Get_param_x0()),_y0(data_file->Get_param_y0()),_a(data_file->Get_param_a()),
_b(data_file->Get_param_b()),_c(data_file->Get_param_c()),_d(data_file->Get_param_d()),
_e(data_file->Get_param_e()),_f(data_file->Get_param_f())
{
  if (data_file->Get_initial_condition_choice() == "gaussian")
  {
    _initial_condition_choice = Gaussian;
  }
  else if (data_file->Get_initial_condition_choice() == "rectangular")
  {
    _initial_condition_choice = Rectangular;
  }
  else
  {
    std::cout << "Invalid initial condition" << std::endl;
    abort();
  }

  if (data_file->Get_velocity_choice() == "uniform")
  {
    _velocity_choice = Uniform;
  }
  else if (data_file->Get_velocity_choice() == "rotational")
  {
    _velocity_choice = Rotational;
  }
  else if (data_file->Get_velocity_choice() == "sinusoidal")
  {
    _velocity_choice = Sinusoidal;
  }
  else
  {
    std::cout << "Invalid velocity" << std::endl;
    abort();
  }
}


double Function::ExactSolution(const double x, const double y, const double t) const
{
  if (_velocity_choice != Uniform)
  {
    std::cout << "The exact solution is not known !" << std::endl;
    abort();
  }
  return InitialCondition(x-_c*t,y-_d*t);
}

double Function::InitialCondition(const double x, const double y) const
{
  switch (_initial_condition_choice)
  {
    case Gaussian:
      return exp( -_a*( pow((x-_x0),2) + pow((y-_y0),2) ) );
      break;
    case Rectangular:
    if (pow((x-_x0),2)+pow((y-_y0),2)<_b)
      return 1.;
    else
      return 0.;
      break;
    default:
      abort();
  }
}

double Function::Velocity_x(const double x, const double y, const double t) const
{
  switch (_velocity_choice)
  {
    case Uniform:
      return _c;
      break;
    case Rotational:
      return _c*y;
      break;
    case Sinusoidal:
      return 0.;
      break;
    default:
      abort();
  }
}

double Function::Velocity_y(const double x, const double y, const double t) const
{
  switch (_velocity_choice)
  {
    case Uniform:
      return _d;
      break;
    case Rotational:
      return _d*x;
      break;
    case Sinusoidal:
      return _e*cos(_f*M_PI*t);
      break;
    default:
      abort();
  }
}


#define _FUNCTION_CPP
#endif
