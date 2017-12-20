#ifndef _ADVECTION_H
#include <string>
#include "Dense"
#include "Sparse"
#include "Function.h"
#include "Mesh2D.h"

class Advection {
	private:
		// Pointeur de la classe Function (accès à la condition initiale, la solution
		// exacte et à la vitesse)
		const Function* _function;
		// Variables géométriques
		const std::vector<Triangle>& _triangles;
		const std::vector<Vertex>& _vertices;
		const std::vector<Edge>& _edges;
		const Eigen::Matrix<double, Eigen::Dynamic, 2>& _tri_center;
		const Eigen::VectorXd & _tri_area;
		const Eigen::Matrix<double, Eigen::Dynamic, 2>& _edg_normal;
		const Eigen::Matrix<double, Eigen::Dynamic, 2>& _edg_center;
		const Eigen::VectorXd & _edg_length;
		// Dossier qui contient les résultats
    const std::string _results;
		// Le choix du flux
		int _numerical_flux; // Centered ou Upwind
  	enum{Centered, Upwind};

		// Écoulement V aux centres des triangles
		Eigen::Matrix<double, Eigen::Dynamic, 2> _V;

		// Vecteur f = g(u,t) (EDO : du/dt = g(u))
		Eigen::VectorXd _f;

	public:
		// Constructeur
		Advection(Function* source_fct, DataFile* data_file, Mesh2D* mesh);

		// calcul de V au centre des triangles
	  void BuildVelocity(const double t);

		// Construit le vecteur f = F(u,t) (EDO : du/dt = F(u,t))
	  void BuildF(const double& t, const Eigen::VectorXd& sol);

		// Renvoie le vecteur f = F(u,t) (EDO : du/dt = F(u,t))
	  Eigen::VectorXd& GetF() {return _f;};

	  // Condition Initiale au centre des triangles
	  Eigen::VectorXd InitialCondition();

	  // Solution exacte au centre des triangles
	  Eigen::VectorXd ExactSolution(double t);

	  // Sauvegarde la solution
	  void SaveSol(const Eigen::VectorXd& sol, int n);
};

#define _ADVECTION_H
#endif
