#ifndef CONSTANTS_H_
#define CONSTANTS_H_
#include <string>

enum class MeshGeometryType {G2D = 0, G3D = 1};
enum class ProblemType {P2DCircle = 0, P3DPipe = 1, P2DMixer = 2, P3DTest4 = 3};
enum class MeshGeneratorType {Salome = 0, FreeFem = 1};

namespace MeshElementsType {
	namespace Salome {
		namespace Tetrahedron {
			const std::string sLin1D = "102",
							  sLin2D = "203",
							  sLin3D = "304";
			const std::string sQuad1D = "103",
							  sQuad2D = "206",
							  sQuad3D = "310";
//			enum class Linear {
//				T1D = 102,
//				T2D = 203,
//				T3D = 304
//			};
//			enum class Quad {
//				T1D = 103,
//				T2D = 206,
//				T3D = 310
//			};
		};
	}
};

#endif //CONSTANTS_H_
