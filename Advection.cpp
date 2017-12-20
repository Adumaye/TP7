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
	// TODO : remplir _V avec les vitesses pour chaque triangle
}

// Construit le vecteur f = F(u,t) (EDO : du/dt = F(u,t))
void Advection::BuildF(const double& t, const Eigen::VectorXd& sol)
{
	_f.resize(_triangles.size());
	BuildVelocity(t);
	// TODO : Evaluer la fonction F(U, t) pour résoudre dt U = F(U,t)
}

// Construit la condition initiale au centre des triangles
VectorXd Advection::InitialCondition()
{
	VectorXd sol0(_triangles.size());
	// TODO : Construire sol0 au centre des triangles
	return sol0;
}

// Solution exacte au centre des triangles
VectorXd Advection::ExactSolution(const double t)
{
	VectorXd exact_sol(_triangles.size());
	//  TODO : calculer la condition initiale au centre des triangles
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
  for (int i = 0 ; i < _triangles.size() ; ++i)
  {
    solution << float(sol[i]) << endl;
  }
  solution << endl;

	solution.close();
}

#define _ADVECTION_CPP
#endif
