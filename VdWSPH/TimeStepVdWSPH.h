#ifndef __TimeStepVdWSPH_h__
#define __TimeStepVdWSPH_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/TimeStep.h"
#include "SimulationDataVdWSPH.h"
#include "SPlisHSPlasH/SPHKernels.h"

namespace SPH
{
	class SimulationDataVdWSPH;

	/** \brief This class implements the van der Waals model
	*
	* References:
	*/
	class TimeStepVdWSPH : public TimeStep
	{
	protected:
		Real m_bParameter;
		Real m_kParameter;
		Real m_T;

		SimulationDataVdWSPH m_simulationData;

		/** Determine the pressure accelerations when the pressure is already known. */
		void computePressureAccels(const unsigned int fluidModelIndex);
		void advectParticles();

		virtual void performNeighborhoodSearchSort();

		virtual void emittedParticles(FluidModel* model, const unsigned int startIndex);
		virtual void initParameters();

	public:
		static int B_PARAM;
		static int K_PARAM;
		static int T;

		TimeStepVdWSPH();
		virtual ~TimeStepVdWSPH(void);

		virtual void step();
		virtual void reset();
		virtual void resize();
	};
}

#endif
