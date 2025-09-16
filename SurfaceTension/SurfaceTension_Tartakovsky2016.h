#ifndef __SurfaceTension_Tartakovsky2016_h__
#define __SurfaceTension_Tartakovsky2016_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/FluidModel.h"
#include "SurfaceTensionBase.h"
#include "SPlishSPlasH/NeighborhoodSearch.h" //!!!

namespace SPH
{
	/** \brief This class implements the surface tension method introduced
	* by TODO
	*
	* References: 
	* - TODO
	*/
	class SurfaceTension_Tartakovsky2016 : public SurfaceTensionBase
	{
	protected:

		// Needed for tuned neighborhoods

		NeighborhoodSearch* m_neighborhoodSearch = nullptr;
		Real m_supportMultiplier = static_cast<Real>(1.0);

		// Virial Pressure

		Real m_virialPG;			  // Global
		std::vector<Real> m_virialPL; // Local
		std::vector<Matrix3r> m_tensors; // Tension tensor 

	public:
		SurfaceTension_Tartakovsky2016(FluidModel *model);
		virtual ~SurfaceTension_Tartakovsky2016(void);

		static NonPressureForceBase* creator(FluidModel* model) { return new SurfaceTension_Tartakovsky2016(model); }

		virtual void step();
		virtual void reset();
		virtual void deferredInit();
		virtual void initParameters();

		static int SUPPORT_MULTIPLIER;

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
