// LiquidFlow.cpp : Defines the entry point for the console application.
/*############################################################################

       Copyright Alexander Epifanov, Inc. 2014 - 2014

       The copyright notice above does not evidence any
       actual or intended publication of such source code.
       The code contains
       Alexander Epifanov's Confidential Restricted Information.

       FILE NAME:              LiquidFlow.cpp

       LOCATIONS:              
		
       OWNER:                  Alexander Epifanov

       DATE CREATED:           06/04/2017

       DESCRIPTION:            LiquidFlow. FEM

############################################################################

       MODIFICATION HISTORY:

 Ver Date     Engineer DR         Change
 --- -------- -------- ---------- ------------------------------------------
 1.0 03/18/14 Alex		1		initial version
############################################################################
 */

#include <time.h>
#include <iostream>
#include "Common.h"
#include "Mesh.h"
#include "Fem.h"
#include "Problem.h"

using namespace std;

int main()
{
	CConfigFileParser cfp(CONFIG_FILENAME);

	int Re = 0,
		Iterations = 0,
		SolverType = 0,
		MeshType = 0,
		Dimension = 0,
		Problem = 0,
		MeshGenerator = 0;
	string MeshFileName;
	real_t Tau = 0;

	Dimension = std::stoi(cfp.getParameter("^Dimension=(\\S+)"));
	Problem = std::stoi(cfp.getParameter("^Problem=(\\S+)"));
	MeshGenerator = std::stoi(cfp.getParameter("^MeshGenerator=(\\S+)"));
	Re = std::stoi(cfp.getParameter("^Re=(\\S+)"));
	Tau = std::stof(cfp.getParameter("^Tau=(\\S+)"));
	Iterations = std::stoi(cfp.getParameter("^Iterations=(\\S+)"));
	MeshFileName = cfp.getParameter("^MeshFile=(\\S+)");

	CMesh *mesh = nullptr;
	switch(MeshGenerator) {
		case (int)MeshGeneratorType::Salome:
			mesh = new CSalomeMesh(MeshFileName);
			break;
		case (int)MeshGeneratorType::FreeFem:
			mesh = new CFreeFemMesh(MeshFileName);
			break;
	};

	mesh->Init((MeshGeometryType)Dimension);
	CFem *fem = nullptr;
	CProblem *pr = NULL;

	if (Dimension == (int)MeshGeometryType::G2D) {
		if (mesh->getPointsNumberPerElement() == 3) {
			fem = new CFemLocalLinear2D(mesh);
			cout << "Linear" << endl;
		}
		if (mesh->getPointsNumberPerElement() == 6) {
			fem = new CFemLocalQuad2D(mesh);
			cout << "Quadratic" << endl;
		}
	}
	else {
		if (mesh->getPointsNumberPerElement() == 4) {
			fem = new CFemLocalLinear3D(mesh);
			cout << "Linear" << endl;
		}
		if (mesh->getPointsNumberPerElement() == 10) {
			fem = new CFemLocalQuad3D(mesh);
			cout << "Quadratic" << endl;
		}
	}

	switch(Problem) {
		case (int)ProblemType::P2DCircle:
			pr = new CProblem2DCircle(mesh);
			break;
		case (int)ProblemType::P2DMixer:
			pr = new CProblem2DMixer(mesh);
			break;
		case (int)ProblemType::P3DPipe:
			pr = new CProblem3DPipe(mesh);
			break;
		case (int)ProblemType::P3DTest4:
			pr = new CProblem3DCubeTest4(mesh);
			break;
	}
	pr->setRe(1.0);
	pr->setTau(Tau);
	fem->init(pr);
	fem->perform(Iterations);

	delete pr;
	delete mesh;
	delete fem;


	return 0;
}
