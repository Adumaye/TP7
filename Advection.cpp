#ifndef _ADVECTION_CPP
#include "Advection.h"
#include <fstream>
#include <iostream>

using namespace std;
using namespace Eigen;

// Constructeur
Advection::Advection(Function* source_fct, DataFile* data_file, Mesh2D* mesh) :
_function(source_fct), _vertices(mesh->GetVertices()),_triangles(mesh->GetTriangles()),
_edges(mesh->GetEdges()), _tri_center(mesh->GetTrianglesCenter()), _tri_area(mesh->GetTrianglesArea()),
_edg_normal(mesh->GetEdgesNormal()), _edg_center(mesh->GetEdgesCenter()), _edg_length(mesh->GetEdgesLength()),
_results(data_file->Get_results())
{
	system(("mkdir -p ./" + _results).c_str());
	system(("rm -f ./" + data_file->Get_results() + "/*.vtk").c_str());
	// Copier le fichier de données dans le dossier résultats
	system(("cp -r ./" + data_file->Get_file_name() + " ./" + _results + "/" + data_file->Get_file_name()).c_str());

	if (data_file->Get_numerical_flux_choice() == "centered")
	_numerical_flux = Centered;
	else
	_numerical_flux = Upwind;
}

// Construit les deux coordonnées du vecteur vitesse au centre des triangles
void Advection::BuildVelocity(const double t)
{
	_V.resize(_triangles.size(),2);
	for (int i=0; i<_triangles.size();i++)
	{
		_V(i,0)=_function->Velocity_x(_tri_center(i,0),_tri_center(i,1),t);
		_V(i,1)=_function->Velocity_y(_tri_center(i,0),_tri_center(i,1),t);
	}
}

// Construit le vecteur f = F(u,t) (EDO : du/dt = F(u,t))
void Advection::BuildF(const double& t, const Eigen::VectorXd& sol)
{
	_f.resize(_triangles.size());
	BuildVelocity(t);
	_f.setZero();
	// boucle sur les arêtes
	for (int i = 0 ; i < _edges.size() ; i++)
	{
		int t1 = _edges[i].GetT1(); int t2= _edges[i].GetT2();
		double scal= _V(t1,0)*_edg_normal(i,0)+_V(t1,1)*_edg_normal(i,1);
		if (t2 == -1)
		{
			if (scal>0.)
				_f(t1)+= ( -1.*_edg_length(i)*scal*sol(t1) ) / _tri_area(t1);
		}
		else
		{ // Arête(2 triangles associés)
			double scal2= _V(t2,0)*-1.*_edg_normal(i,0)+_V(t2,1)*-1.*_edg_normal(i,1);

			if (_numerical_flux == Centered)
			{
				_f(t1) +=  ( -1.*_edg_length(i)*(scal*sol(t1)+scal2*sol(t2))/2. ) / _tri_area(t1);
				_f(t2) +=  ( _edg_length(i)*(scal*sol(t1)+scal2*sol(t2))/2. ) / _tri_area(t2);
			}
			else
			{ // Flux décentré : remplir _f[t1] et _f[t2]
				if (scal>0.)
				{
					_f(t1) +=  ( -1.*_edg_length(i)*scal*sol(t1) ) / _tri_area(t1);
					_f(t2) +=  ( _edg_length(i)*scal*sol(t1) ) / _tri_area(t2);
				}
				else
				{
					_f(t1) +=  ( -1.*_edg_length(i)*scal*sol(t2) ) / _tri_area(t1);
					_f(t2) +=  ( _edg_length(i)*scal*sol(t2) ) / _tri_area(t2);
				}
			}
		}
	}
}

// Construit la condition initiale au centre des triangles
VectorXd Advection::InitialCondition()
{
	VectorXd sol0(_triangles.size());
	for (int i=0; i<_triangles.size();i++)
	{
		sol0(i)=_function->InitialCondition(_tri_center(i,0),_tri_center(i,1));
	}
	return sol0;
}

// Solution exacte au centre des triangles
VectorXd Advection::ExactSolution(const double t)
{
	VectorXd exact_sol(_triangles.size());
	for (int i=0; i<_triangles.size();i++)
	{
		exact_sol(i)=_function->ExactSolution(_tri_center(i,0),_tri_center(i,1),t);
	}
	return exact_sol;
}

// Sauvegarde la solution
void Advection::SaveSol(const Eigen::VectorXd& sol, int n)
{
	string name_file = _results + "/solution_" + std::to_string(n) + ".vtk";

	int nb_vert = _vertices.size();

	assert((sol.size() == _triangles.size()) && "The size of the solution vector is not the same than the number of _triangles !");

	ofstream solution;
	solution.open(name_file, ios::out);
	solution.precision(7);

	solution << "# vtk DataFile Version 3.0 " << endl;
	solution << "2D Unstructured Grid" << endl;
	solution << "ASCII" << endl;
	solution << "DATASET UNSTRUCTURED_GRID" << endl;

	solution << "POINTS " << nb_vert << " float " << endl;
	for (int i = 0 ; i < nb_vert ; ++i)
	{
		solution << ((_vertices[i]).GetCoor())[0] << " " << ((_vertices[i]).GetCoor())[1] << " 0." << endl;
	}
	solution << endl;

	solution << "CELLS " << _triangles.size() << " " << _triangles.size()*4 << endl;
	for (int i = 0 ; i < _triangles.size() ; ++i)
	{
		solution << 3 << " " << ((_triangles[i]).GetVertices())[0] << " " << ((_triangles[i]).GetVertices())[1]
		<< " " << ((_triangles[i]).GetVertices())[2] << endl;
	}
	solution << endl;

	solution << "CELL_TYPES " << _triangles.size() << endl;
	for (int i = 0 ; i < _triangles.size() ; ++i)
	{
		solution << 5 << endl;
	}
	solution << endl;

	solution << "CELL_DATA " << _triangles.size() << endl;
	solution << "SCALARS sol float 1" << endl;
	solution << "LOOKUP_TABLE default" << endl;

	double norm = 0.;
	for (int i = 0 ; i < _triangles.size() ; ++i)
	{
		solution << float(sol[i]) << endl;
		norm+=_tri_area(i)*pow(sol(i),2);
	}
	solution << endl;
	solution.close();
	cout << sqrt(norm)<<endl;
}

#define _ADVECTION_CPP
#endif
