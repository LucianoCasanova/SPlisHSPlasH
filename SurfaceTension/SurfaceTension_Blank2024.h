#ifndef __SurfaceTension_Blank2024_h__
#define __SurfaceTension_Blank2024_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/FluidModel.h"
#include "SurfaceTensionBase.h"
#include "SPlishSPlasH/NeighborhoodSearch.h"

namespace SPH
{
	/** \brief This class implements the surface tension method introduced
	* by TODO
	*
	* References:
	* TODO
	*/
	class SurfaceTension_Blank2024 : public SurfaceTensionBase
	{
	protected:
		std::vector<Vector3r> m_normals;
		std::vector<Vector3r> m_betterNormals;
		std::vector<Real> m_curvatures;

		void computeRawNormals();
		void smoothAndCapNormals();
		Matrix3r computeGradCorrectionMatrix(const unsigned int i);

	public:
		SurfaceTension_Blank2024(FluidModel* model);
		virtual ~SurfaceTension_Blank2024(void);

		static NonPressureForceBase* creator(FluidModel* model) { return new SurfaceTension_Blank2024(model); }

		virtual void step();
		virtual void reset();
	};
}

#endif
