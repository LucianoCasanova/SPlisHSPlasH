#ifndef __TimeStepISPH_h__
#define __TimeStepISPH_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/TimeStep.h"
#include "SimulationDataISPH.h"
#include "SPlisHSPlasH/SPHKernels.h"

namespace SPH
{
	class SimulationDataISPH;

	/** \brief This class implements the Incompressible SPH approach introduced
	 * by Cummings et al. [CR+99]. A correction for free surface by [NT+14] is also available.
	*
	* References:
	* - [CR+99] Cummins, S. J., & Rudman, M. (1999). An SPH projection method. Journal of Computational Physics, 152(2), 584–607. https://doi.org/10.1006/jcph.1999.6246
	* - [NT+14] Nair, P., & Tomar, G. (2014). An improved free surface modeling for incompressible SPH. Computers & Fluids, 102, 304–314. https://doi.org/10.1016/j.compfluid.2014.07.006
	*/
	class TimeStepISPH : public TimeStep
	{
	protected:
		SimulationDataISPH m_simulationData;

		bool m_useSurfaceCorrection;
		Real m_omega;

		bool m_kObtained;

		void calculateK(const unsigned int fluidModelIndex);
		void predictAdvection(const unsigned int fluidModelIndex);
		void pressureSolve();
		void pressureSolveIteration(const unsigned int fluidModelIndex, Real& avg_density_err);
		void integration(const unsigned int fluidModelIndex);

		/** Determine the pressure accelerations when the pressure is already known. */
		void computePressureAccels(const unsigned int fluidModelIndex);

		virtual void performNeighborhoodSearchSort();

		virtual void emittedParticles(FluidModel* model, const unsigned int startIndex);
		virtual void initParameters();

	public:
		static int P_AMB;
		static int USE_SURFACE_CORRECTION;
		static int OMEGA;

		TimeStepISPH();
		virtual ~TimeStepISPH(void);

		virtual void step();
		virtual void reset();
		virtual void resize();

		SimulationDataISPH& getSimulationData() { return m_simulationData; };
	};
}

#endif
