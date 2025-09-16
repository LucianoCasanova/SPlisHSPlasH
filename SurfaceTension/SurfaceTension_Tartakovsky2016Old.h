#ifndef __SurfaceTension_Tartakovsky2016_h__
#define __SurfaceTension_Tartakovsky2016_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/FluidModel.h"
#include "SurfaceTensionBase.h"

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
	public:
		SurfaceTension_Tartakovsky2016(FluidModel *model);
		virtual ~SurfaceTension_Tartakovsky2016(void);

		static NonPressureForceBase* creator(FluidModel* model) { return new SurfaceTension_Tartakovsky2016(model); }

		virtual void step();
		virtual void reset();
	};
}

#endif
