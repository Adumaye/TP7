#ifndef _MESH_2D_CPP

#include "Mesh2D.h"
#include<iostream>
#include<fstream>
#include <assert.h>

using namespace std;

Vertex::Vertex()
{
  _v_coor(0) = -10000; _v_coor(1) = -10000; _ref = -1;
}

Vertex::Vertex(double x, double y, int ref) : _ref(ref)
{
  _v_coor(0) = x; _v_coor(1) = y;
}

void Vertex::Print() const
{
  cout << "[x, y] = [" << _v_coor(0) << " " << _v_coor(1) << "];" << endl;
  cout << "ref = " << _ref << endl;
}

Edge::Edge()
{
  _v_edge(0) = -1; _v_edge(1) = -1; _ref = -1;
}

Edge::Edge(int vertex1, int vertex2, int ref) : _ref(ref)
{
  // sort
  if (vertex1 > vertex2)
  {
    _v_edge(0) = vertex2;
    _v_edge(1) = vertex1;
  }
  else
  {
    _v_edge(0) = vertex1;
    _v_edge(1) = vertex2;
  }
  _t1 = -1;
  _t2 = -1;
}

void Edge::Print() const
{
  cout << "[pt1, pt2] = [" << _v_edge(0) << " " << _v_edge(1) << "];" << endl;
  cout << "[t1, t2] = [" << _t1 << " " << _t2 << "];" << endl;
  cout << "ref = " << _ref << endl;
}

Triangle::Triangle()
{
  _v_triangle(0) = -1; _v_triangle(1) = -1; _v_triangle(2) = -1; _ref = -1;
}

Triangle::Triangle(int vertex1, int vertex2, int vertex3, int ref) : _ref(ref)
{
  _v_triangle(0) = vertex1; _v_triangle(1) = vertex2; _v_triangle(2) = vertex3;
}

void Triangle::Print() const
{
  cout << "[pt1, pt2, pt3] = [" << _v_triangle(0) << " " << _v_triangle(1) << " " << _v_triangle(2) << "];" << endl;
  cout << "ref = " << _ref << endl;
}

Mesh2D::Mesh2D()
{
}

void Mesh2D::BuildTrianglesCenterAndArea()
{
  _tri_center.resize(_triangles.size(),2);
  _tri_area.resize(_triangles.size());

  Eigen::VectorXd y;
  Eigen::VectorXd x;
  x.resize(3);y.resize(3);

  for (int i=0; i<_triangles.size();i++)
  {
    for (int j=0; j<3; j++)
    {
      int pj = _triangles[i].GetVertices()(j);
       x(j) = _vertices[pj].GetCoor()(0);
       y(j) = _vertices[pj].GetCoor()(1);
    }
    _tri_center(i,0)=(x(0)+x(1)+x(2))/3;
    _tri_center(i,1)=(y(0)+y(1)+y(2))/3;

    _tri_area(i)= abs( (x(1)-x(0))*(y(2)-y(0))-(y(1)-y(0))*(x(2)-x(0)) ) / 2.;
  }
  //  Construire le vecteur contenant les aires de chaque triangle
  // _tri_area(i) correspond à l'aire du triangle i
}

void Mesh2D::BuildEdgesNormalLengthAndCenter()
{
  _edg_center.resize(_edges.size(),2);
  _edg_length.resize(_edges.size());
  _edg_normal.resize(_edges.size(),2);

  Eigen::Vector2d x;
  Eigen::Vector2d y;
  for (int i=0; i<_edges.size();i++)
  {
    for (int j=0; j<2; j++)
    {
      int pj = _edges[i].GetVertices()(j);
      x[j] = _vertices[pj].GetCoor()(0);
      y[j] = _vertices[pj].GetCoor()(1);
    }
    _edg_center(i,0)=(x(0)+x(1))/2.;
    _edg_center(i,1)=(y(0)+y(1))/2.;

    _edg_length(i)=sqrt(pow(x(1)-x(0),2)+pow(y(1)-y(0),2));

    _edg_normal(i,0)= (y(0)-y(1))/_edg_length(i);
    _edg_normal(i,1)= (x(1)-x(0))/_edg_length(i);

    int t1 = _edges[i].GetT1();

    Eigen::Vector2d diff = _edg_center.row(i) - _tri_center.row(t1);

    if (diff(0)*_edg_normal(i,0)+diff(1)*_edg_normal(i,1) < 0)
      _edg_normal.row(i) *= -1;
  }
  // Construire la matrice contenant les coordonnées des normales unitaires de chaque arête
  // _edg_normal(i,0) correspond à la cooordonnée x de la normale unitaires de l'arête i
  // _edg_normal(i,1) correspond à la cooordonnée y de la normale unitaires de l'arête i
}

// methode interne qui rajoute une arete
void Mesh2D::AddSingleEdge(const Edge& edge, int ne, vector<int>& head_minv, vector<int>& next_edge, int& nb_edges)
{
  int n1 = edge.GetVertices()(0);
  int n2 = edge.GetVertices()(1);
  int ref = edge.GetReference();

  bool exist = false;
  // we look at the list of edges leaving from n1
  // if we find the same edge than n1->n2 we add the edge
  for (int e = head_minv[n1]; e != -1; e = next_edge[e])
  {
    if (_edges[e].GetVertices()(1) == n2)
    {
      if (ne >= 0)
      {
        _edges[e].AddTriangle(ne);
      }
      exist = true;
    }
  }

  // if the edge has not been found, we create it
  if (!exist)
  {
    // we initialize the edge
    _edges[nb_edges] = Edge(n1, n2, ref);
    if (ne >= 0)
    {
      _edges[nb_edges].AddTriangle(ne);
    }
    // we update the arrays next_edge and head_minv
    next_edge[nb_edges] = head_minv[n1];
    head_minv[n1] = nb_edges;
    nb_edges++;
  }
}

void Mesh2D::ReadMesh(string name_mesh)
{
  ifstream mesh_file(name_mesh.data());
  if (!mesh_file.is_open())
  {
    cout << "Unable to open file " << name_mesh << endl;
    abort();
  }
  else
  {
    cout << "-------------------------------------------------" << endl;
    cout << "Reading mesh: " << name_mesh << endl;
  }

  string file_line;
  vector<Edge> edges_boundary;
  int dim = 3;

  while (!mesh_file.eof())
  {
    getline(mesh_file, file_line);
    if (file_line.find("Dimension") != std::string::npos)
    {
      mesh_file >> dim;
    }
    else if (file_line.find("Vertices") != std::string::npos)
    {
      int nb_vertices(0);
      mesh_file >> nb_vertices;
      cout << "Number of vertices  (" << nb_vertices << ")" << endl;
      _vertices.resize(nb_vertices);
      for (int i = 0 ; i < nb_vertices ; ++i)
      {
        double x,y,z; int ref;
        mesh_file >> x >> y >> z >> ref;
        _vertices[i] = Vertex(x, y, ref);
      }
    }
    else if (file_line.find("Edges") != std::string::npos)
    {
      int nb_edges(0);
      mesh_file >> nb_edges;
      cout << "Number of edges (" << nb_edges << ")" << endl;
      edges_boundary.resize(nb_edges);
      int n1, n2, ref;
      for (int i = 0 ; i < nb_edges ; ++i)
      {
        mesh_file >> n1 >> n2 >> ref;
        n1--; n2--;
        edges_boundary[i] = Edge(n1, n2, ref);
      }
    }
    else if (file_line.find("Triangles") != std::string::npos)
    {
      int nb_triangles(0);
      mesh_file >> nb_triangles;
      cout << "Number of triangles (" << nb_triangles << ")" << endl;
      _triangles.resize(nb_triangles);
      for (int i = 0 ; i < nb_triangles ; ++i)
      {
        int vertex1, vertex2, vertex3, ref;
        mesh_file >> vertex1 >> vertex2 >> vertex3 >> ref;
        vertex1--; vertex2--; vertex3--;
        _triangles[i] = Triangle(vertex1, vertex2, vertex3, ref);
      }
    }
  }

  cout << "---------Edges and Associated Triangles----------" << endl;
  // Toutes les aretes exterieures du maillage sont presentes
  int nb_edges = (3*_triangles.size() + edges_boundary.size())/2;
  _edges.resize(nb_edges);

  int nb_vertices = _vertices.size();
  vector<int> head_minv(nb_vertices, -1);
  vector<int> next_edge(nb_edges, -1);

  // on rajoute d'abord les aretes du bord
  nb_edges = 0;
  for (int i = 0; i < edges_boundary.size(); i++)
  {
    this->AddSingleEdge(edges_boundary[i], -1, head_minv, next_edge, nb_edges);
  }

  // ensuite les aretes interieures
  for (int i = 0; i < _triangles.size(); i++)
  {
    const Eigen::Vector3i& nv = _triangles[i].GetVertices();
    for (int j = 0; j < 3; j++)
    {
      Edge edge(nv(j), nv((j+1)%3), 0);
      AddSingleEdge(edge, i, head_minv, next_edge, nb_edges);
    }
  }

  cout << "-----------Triangles center and area-------------" << endl;
  BuildTrianglesCenterAndArea();

  cout << "------------Edges Normal --------------" << endl;
  BuildEdgesNormalLengthAndCenter();

  cout << "-------------------------------------------------" << endl;
}

#define _MESH_2D_CPP
#endif
