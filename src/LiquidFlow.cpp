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

       DATE CREATED:           04/27/2013

       DESCRIPTION:            LiquidFlow. FEM

############################################################################

       MODIFICATION HISTORY:

 Ver Date     Engineer DR         Change
 --- -------- -------- ---------- ------------------------------------------
 1.0 03/18/14 Alex		1		initial version
############################################################################
 */

#define WIN

#include <time.h>
#include <stdio.h>
#include "Mesh.h"
#include "Point.h"
//#include "CommonSolver.h"
#include "Common.h"

void meshTest1(CExtGenMesh *m_mb);

int main()
{

    //CBaseFEMSolver *m_bs = NULL;
	CExtGenMesh *m_mb = NULL;
	CConfigFileParser cfp(CONFIG_FILENAME);

	int Re = 0,
		H = 0,
		n = 0,
		nz = 0,
		Iterations = 0,
		SolverType = 0,
		MeshType = 0,
		MeshGenType = 0,
		res = OK;
	real_t R = 0,
		   Tau = 0,
		   StartVz = 10000;
	char *fn;
	char mesh_folder[256];
	char mesh_mask[256];
	char *cfg_str;

	Re = atoi(cfp.GetParameter("^Re=(\\S+)"));
	Tau = atof(cfp.GetParameter("^Tau=(\\S+)"));
	Iterations = atoi(cfp.GetParameter("^Iterations=(\\S+)"));
	SolverType = atoi(cfp.GetParameter("^SolverType=(\\S+)"));
	MeshType = atoi(cfp.GetParameter("^MeshType=(\\S+)"));


#ifdef WIN
	clock_t start = clock();
	clock_t stop;
#endif
#ifdef LINUX
	struct timespec start;
	struct timespec stop;
	clock_gettime(CLOCK_REALTIME, &start);
#endif
	printf ("MeshType %d\n", MeshType);
	switch(MeshType)
	{
		case 1:
			//strcpy(fn, cfp.GetParameter("^MeshFile=(\\S+)"));
			fn = strdup(cfp.GetParameter("^MeshFile=(\\S+)"));			
			printf("DEBUG: CExtGenMesh test1 mesh %s\n", fn);
			m_mb = new CExtGenMesh(fn, MG_SALOME, SCHEME_ORDER_1);			
			if (!m_mb) 
				printf("new CExtGenMesh Error!\n");
			m_mb->SalomeDat3DParser();
			meshTest1(m_mb);
			break;
		case 2:
			//group of files
			strcpy(mesh_folder, cfp.GetParameter("^MeshFolder=(\\S+)"));
			strcpy(mesh_mask, cfp.GetParameter("^StartMeshMask=(\\S+)"));
			strcpy(&fn[0], get_nextfile(mesh_folder, mesh_mask));
			//fn = get_nextfile(mesh_folder, mesh_mask);
			m_mb = new CExtGenMesh(fn, MG_SALOME, SCHEME_ORDER_1);
			if (!m_mb) printf("Error!\n");
			printf("DEBUG: CExtGenMesh\n");
			break;
		default: break;
	}



#ifdef WIN
	stop = clock();
	printf("Mesh Generate time seconds: %ul\n", (stop - start)/*CLOCKS_PER_SEC*/);
	start = clock();
#endif

#ifdef LINUX
	clock_t stop = clock();
	clock_gettime(CLOCK_REALTIME, &stop);
	printf("Mesh Generate time seconds: %ul\n", (stop.tv_sec - start.tv_sec) + (float)(stop.tv_nsec - start.tv_nsec)/1000000000);
	clock_gettime(CLOCK_REALTIME, &start);
#endif
	printf ("SolverType %d\n", SolverType);
	// switch(SolverType)
	// 		{
	// 			case 1: //PipeSolver
	// 				m_bs = new_CSimpleTriangleSolver(m_mb, Re, Tau);
	// 				printf("DEBUG: Create PipeSolver\n");
	// 				Initi_SimpleTriangleSolver((CSimpleTriangleSolver*)m_bs);
	// 				printf("DEBUG: Init PipeSolver\n");
	// 				for (int i = 1; i < Iterations; i++)
	// 				{
	// 					Run_SimpleTriangleSolver((CSimpleTriangleSolver*)m_bs, Iterations - 1, i);
	// 					printf ("Iterations %d\n", i);
	// 				}
	// 				printf("DEBUG: Done PipeSolver\n");
	// 				break;
	// 			case 2: //TsunamiSolver
	// 				m_bs = new_CSimpleTriangleTsunamiSolver(m_mb, Re, Tau, StartVz);
	// 				printf("DEBUG: Create TsunamiSolver\n");
	// 				Init_SimpleTriangleTsunamiSolver((CSimpleTriangleTsunamiSolver*)m_bs);
	// 				printf("DEBUG: Init TsunamiSolver\n");
	// 				for (int i = 1; i < Iterations; i++)
	// 				{
	// 					Run_SimpleTriangleTsunamiSolver((CSimpleTriangleTsunamiSolver*)m_bs, Iterations - 1, i);
	// 				}
	// 				printf("DEBUG: Done TsunamiSolver\n");
	// 				break;

	// 			case 3: //PipeSolver2D
	// 				m_bs = new_CSimpleTriangleSolver2D(m_mb, Re, Tau);
	// 				printf("DEBUG: CSimpleTriangleSolver2D\n");
	// 				res = Init_SimpleTriangleSolver2D((CSimpleTriangleSolver2D*)m_bs);
	// 				if (OK == res)
	// 				{
	// 					printf("DEBUG: Init CSimpleTriangleSolver2D\n");
	// 					for (int i = 1; i < Iterations; i++)
	// 					{
	// 						res = Run_SimpleTriangleSolver2D((CSimpleTriangleSolver2D*)m_bs, Iterations - 1, i);
	// 						if (OK != res) break;
	// 						printf("iter %d\n", i);
	// 					}
	// 					printf("DEBUG: Done CSimpleTriangleSolver2D\n");
	// 				}
	// 				break;

	// 			case 4: //MixerSolver2D
	// 				m_bs = new_CSimpleMixerSolver2D(m_mb, Re, Tau);
	// 				printf("DEBUG: Create MixerSolver2D\n");
	// 				res = Init_SimpleMixerSolver2D((CSimpleMixerSolver2D*)m_bs);
	// 				if (OK == res)
	// 				{
	// 					printf("DEBUG: Init MixerSolver2D\n");
	// 					for (int i = 1; i < Iterations; i++)
	// 					{
	// 						res = Run_SimpleMixerSolver2D((CSimpleMixerSolver2D*)m_bs, Iterations - 1, i);
	// 						if (OK != res) break;
	// 						printf("iter %d\n", i);
	// 					}
	// 					printf("DEBUG: Done MixerSolver2D\n");
	// 				}
	// 				break;

	// 			case 5: //Submarine Solver
	// 				m_bs = new_CSimpleSubmarineSolver3D(m_mb, Re, Tau);
	// 				printf("DEBUG: Create new_CSimpleSubmarineSolver3D\n");
	// 				res = Init_SimpleSubmarineSolver3D((CSimpleSubmarineSolver3D*)m_bs);
	// 				if (OK == res)
	// 				{
	// 					printf("DEBUG: Init_SimpleSubmarineSolver3D\n");
	// 					for (int i = 1; i < Iterations; i++)
	// 					{
	// 						res = Run_SimpleSubmarineSolver3D((CSimpleSubmarineSolver3D*)m_bs, Iterations - 1, i);
	// 						if (OK != res) break;
	// 						printf("iter %d\n", i);
	// 					}
	// 					printf("DEBUG: Done SimpleSubmarineSolver3D\n");
	// 				}
	// 				break;

	// 			case 6: //SuperMixer Solver
	// 				m_bs = new_CSimpleSuperMixerSolver(m_mb, Re, Tau);
	// 				printf("DEBUG: Create new_CSimpleSuperMixerSolver\n");
	// 				res = Init_SimpleSuperMixerSolver((CSimpleSuperMixerSolver*)m_bs);
	// 				if (OK == res)
	// 				{
	// 					printf("DEBUG: Init_SimpleSuperMixerSolver\n");
	// 					for (int i = 1; i < Iterations; i++)
	// 					{
	// 						res = Run_SimpleSuperMixerSolver((CSimpleSuperMixerSolver*)m_bs, Iterations - 1, i);
	// 						if (OK != res) break;
	// 						printf("iter %d\n", i);
	// 					}
	// 					printf("DEBUG: Done SimpleSuperMixerSolver\n");
	// 				}
	// 				break;

	// 			case 7: //SuperMixer Solver
	// 				m_bs = new_CSimpleSuperMixerSolver2D(m_mb, Re, Tau);
	// 				printf("DEBUG: Create new_CSimpleSuperMixerSolver2D\n");
	// 				res = Init_SimpleSuperMixerSolver2D((CSimpleSuperMixerSolver2D*)m_bs);
	// 				if (OK == res)
	// 				{
	// 					printf("DEBUG: Init_SimpleSuperMixerSolver2D\n");
	// 					for (int i = 1; i < Iterations; i++)
	// 					{
	// 						res = Run_SimpleSuperMixerSolver2D((CSimpleSuperMixerSolver2D*)m_bs, Iterations - 1, i);
	// 						if (OK != res) break;
	// 						printf("iter %d\n", i);
	// 					}
	// 					printf("DEBUG: Done SimpleSuperMixerSolver2D\n");
	// 				}
	// 				break;

	// 			case 8: //Poisson Solver
	// 				m_bs = new_SimplePoissonSolver2D(m_mb, Re, Tau);
	// 				printf("DEBUG: Create new_SimplePoissonSolver2D\n");
	// 				res = Init_SimplePoissonSolver2D((CSimplePoissonSolver2D*)m_bs);
	// 				if (OK == res)
	// 				{
	// 					printf("DEBUG: Run_SimplePoissonSolver2D\n");
	// 					for (int i = 1; i < Iterations; i++)
	// 					{
	// 						res = Run_SimplePoissonSolver2D((CSimplePoissonSolver2D*)m_bs, Iterations - 1, i);
	// 						if (OK != res) break;
	// 						printf("iter %d\n", i);
	// 					}
	// 					printf("DEBUG: Done SimplePoissonSolver2D\n");
	// 				}
	// 				break;
	// 			#ifdef LIBMAGMASPARSE
	// 			case 9: 
	// 				m_bs = new_CSimpleTrinagleSparseSolver3D(m_mb, Re, Tau);
	// 				printf("DEBUG: Create new_CSimpleTrinagleSparseSolver3D\n");
	// 				res = Init_CSimpleTrinagleSparseSolver3D((CSimpleTrinagleSparseSolver3D*)m_bs);
	// 				if (OK == res)
	// 				{
	// 					printf("DEBUG: Run_CSimpleTrinagleSparseSolver3D\n");
	// 					// for (int i = 1; i < Iterations; i++)
	// 					// {
	// 					// 	res = Run_CSimpleTrinagleSparseSolver3D((CSimpleTrinagleSparseSolver3D*)m_bs, Iterations - 1, i);
	// 					// 	if (OK != res) break;
	// 					// 	printf("iter %d\n", i);
	// 					// }
	// 					// printf("DEBUG: Done CSimpleTrinagleSparseSolver3D\n");
	// 				}
	// 				break;					
	// 			#endif

	// 		}
#ifdef WIN
	stop = clock();
	printf("Calculations time seconds: %ul\n", (stop - start)/*/CLOCKS_PER_SEC*/);
#endif

#ifdef LINUX
	clock_gettime(CLOCK_REALTIME, &stop);
	printf("Calculations time seconds: %f\n", (stop.tv_sec - start.tv_sec) + (float)(stop.tv_nsec - start.tv_nsec)/1000000000);
#endif
	printf ("finish\n");	
	//delete m_bs;
	//if (fn) free(fn);
	if (m_mb) delete m_mb;
	

	return 0;
}

void meshTest1(CExtGenMesh *m_mb)
{
	printf ("meshTest1\n");
	printf ("GetElnumber() %d\n", m_mb->GetElnumber());
	printf ("GetPtnumber() %d\n", m_mb->GetPtnumber());
	printf ("GetTriangleNumber() %d\n", m_mb->GetTriangleNumber());
	m_mb->PrintElements("elements_test.txt");
}
