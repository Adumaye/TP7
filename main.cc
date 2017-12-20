#include <iostream>
#include <fstream>
#include <chrono>
#include "TimeScheme.h"

using namespace std;
using namespace Eigen;

int main(int argc, char** argv)
{

  if (argc < 2)
  {
    cout << "Please, enter the name of your data file." << endl;
    abort();
  }
  const string data_file_name = argv[1];

  // ----------------------- Fichier de données --------------------------------
  DataFile* data_file = new DataFile(data_file_name);
  data_file->ReadDataFile();
  // ---------------------------------------------------------------------------

  // ------------------Définition du nombre d'itérations------------------------
  int nb_iterations = int(ceil((data_file->Get_tfinal()-data_file->Get_t0())/data_file->Get_dt()));
  data_file->Adapt_dt(data_file->Get_tfinal() / nb_iterations);
  // ---------------------------------------------------------------------------

  // ---------------------------- Résolution  ----------------------------------
  Mesh2D* mesh = new Mesh2D();
  mesh->ReadMesh(data_file->Get_mesh_name());
  Function* source_fct = new Function(data_file);
  Advection* adv = new Advection(source_fct, data_file, mesh);
  TimeScheme* time_scheme = NULL;

  if (data_file->Get_scheme() == "RungeKutta")
    time_scheme = new RungeKuttaScheme();
  else
    time_scheme = new EulerScheme();

  cout << "-------------------------------------------------" << endl;
  cout << "Search u such that : " << endl;
  cout << "dt u + v . grad u = 0,  B" << endl;
  cout << "-------------------------------------------------" << endl;

  // Démarrage du chrono
  auto start = chrono::high_resolution_clock::now();
  time_scheme->Initialize(data_file, adv); // Initialisation
  time_scheme->SaveSolution(0); // Sauvegarde condition initiale

  for (int n = 1; n <= nb_iterations; n++) // Boucle en temps
  {
    time_scheme->Advance();
    time_scheme->SaveSolution(n);
  }
  // Fin du chrono
  auto finish = chrono::high_resolution_clock::now();
  double t = chrono::duration_cast<chrono::milliseconds>(finish-start).count();
  // Affichage du résultat
  cout << "Cela a pris "<< t << " milliseconds" << endl;
  // ---------------------------------------------------------------------------


  //------------------- Si on connait la solution exacte -----------------------
  if (data_file->Get_velocity_choice() == "uniform")
  {
    VectorXd exact_sol = adv->ExactSolution(data_file->Get_tfinal());
    VectorXd approx_sol = time_scheme->GetIterateSolution();
    double error = ((approx_sol-exact_sol).array().abs()).maxCoeff();
    cout << "Erreur = " << error << endl;
    cout << "-------------------------------------------------" << endl;
    ofstream error_file;
    const string error_file_name = data_file->Get_results()+"/error.txt";
    error_file.open(error_file_name, ios::out);
    error_file << error << endl;
    error_file.close();
  }
  // ---------------------------------------------------------------------------

  delete time_scheme;
  delete adv;
  delete data_file;
  delete source_fct;

  return 0;
}
