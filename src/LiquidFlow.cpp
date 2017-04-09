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

int main()
{
	CConfigFileParser cfp(CONFIG_FILENAME);

	int Re = 0,
		Iterations = 0,
		SolverType = 0,
		MeshType = 0;

	Re = std::stoi(cfp.getParameter("^Re=(\\S+)"));
	Iterations = std::stoi(cfp.getParameter("^Iterations=(\\S+)"));
	SolverType = std::stoi(cfp.getParameter("^SolverType=(\\S+)"));
	MeshType = std::stoi(cfp.getParameter("^MeshType=(\\S+)"));

	CMesh *mesh = new CSalomeMesh("../mesh/Mesh_box.dat");
	mesh->Init();

	delete (mesh);
	return 0;
}