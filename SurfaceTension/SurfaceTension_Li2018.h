#ifndef __SurfaceTension_Li2018_h__
#define __SurfaceTension_Li2018_h__

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
	class SurfaceTension_Li2018 : public SurfaceTensionBase
	{
	protected:

		// Method

		Real m_a1, m_a2;
		Real m_h1, m_h2;
		Gaussian3 m_kernel1, m_kernel2;
		WendlandC4 m_boundaryKernel;

		// Needed for tuned neighborhoods

		NeighborhoodSearch* m_neighborhoodSearch = nullptr;

		// Virial Pressure

		Real m_virialPG;			  // Global
		std::vector<Real> m_virialPL; // Local
		std::vector<Matrix3r> m_tensors; // Tension tensor 

	public:
		SurfaceTension_Li2018(FluidModel* model);
		virtual ~SurfaceTension_Li2018(void);

		static NonPressureForceBase* creator(FluidModel* model) { return new SurfaceTension_Li2018(model); }

		virtual void step();
		virtual void reset();
		virtual void deferredInit();
		virtual void initParameters();

		static int A1;
		static int A2;
		static int H1;
		static int H2;

		FORCE_INLINE const Real getVirialP(const unsigned int i) const
		{
			return m_virialPL[i];
		}

		FORCE_INLINE Real& getVirialP(const unsigned int i)
		{
			return m_virialPL[i];
		}

		FORCE_INLINE void setVirialP(const unsigned int i, const Real p)
		{
			m_virialPL[i] = p;
		}

		FORCE_INLINE const Matrix3r getTensor(const unsigned int i) const
		{
			return m_tensors[i];
		}

		FORCE_INLINE Matrix3r& getTensor(const unsigned int i)
		{
			return m_tensors[i];
		}

		FORCE_INLINE void setTensor(const unsigned int i, const Matrix3r t)
		{
			m_tensors[i] = t;
		}
	};
}

#endif
