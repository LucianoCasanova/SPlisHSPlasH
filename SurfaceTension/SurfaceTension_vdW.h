#ifndef __SurfaceTension_vdW_h__
#define __SurfaceTension_vdW_h__

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
	class SurfaceTension_vdW : public SurfaceTensionBase
	{
	protected:
		NeighborhoodSearch* m_neighborhoodSearch = nullptr;
		unsigned int getNeighbor(const unsigned int pointSetIndex, const unsigned int neighborPointSetIndex, const unsigned int index, const unsigned int k) const;
		unsigned int numberOfNeighbors(const unsigned int pointSetIndex, const unsigned int neighborPointSetIndex, const unsigned int index) const;

	public:
		SurfaceTension_vdW(FluidModel* model);
		virtual ~SurfaceTension_vdW(void);

		static NonPressureForceBase* creator(FluidModel* model) { return new SurfaceTension_vdW(model); }

		virtual void step();
		virtual void reset();
	};
}

#endif
