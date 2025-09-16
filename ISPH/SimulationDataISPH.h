#ifndef __SimulationDataISPH_h__
#define __SimulationDataISPH_h__

#include "SPlisHSPlasH/Common.h"
#include <vector>
#include "SPlisHSPlasH/FluidModel.h"

namespace SPH
{
	/** \brief Simulation data which is required by the method Incompressible SPH introduced
	* by Cummings et al. [CR+99]. A correction for free surface by [NT+14] is also available.
	*
	* References:
	* - [CR+99] Cummins, S. J., & Rudman, M. (1999). An SPH projection method. Journal of Computational Physics, 152(2), 584–607. https://doi.org/10.1006/jcph.1999.6246
	* - [NT+14] Nair, P., & Tomar, G. (2014). An improved free surface modeling for incompressible SPH. Computers & Fluids, 102, 304–314. https://doi.org/10.1016/j.compfluid.2014.07.006
	*/
	class SimulationDataISPH
	{
	public:
		SimulationDataISPH();
		virtual ~SimulationDataISPH();

	protected:
		// Surface Correction
		std::vector<std::vector<Real>> m_Shi;
		std::vector<std::vector<Real>> m_corrected;
		std::vector<std::vector<Real>> m_predRho;
		std::vector<std::vector<Real>> m_errPredRho;
		std::vector<Real> m_K;
		std::vector<std::vector<Real>> m_pZero;
		Real m_pAmb;
		// Cummings
		std::vector<std::vector<Real>> m_aii;
		std::vector<std::vector<Real>> m_pressure;
		std::vector<std::vector<Real>> m_lastPressure;
		std::vector<std::vector<Vector3r>> m_pressureAccel;
		std::vector<std::vector<Vector3r>> m_prevPos;
		std::vector<std::vector<Real>> m_si;
		std::vector<std::vector<Real>> m_aij_pj;

	public:
		/** Initialize the arrays containing the particle data.
		*/
		virtual void init();

		/** Release the arrays containing the particle data.
		*/
		virtual void cleanup();

		/** Reset the particle data.
		*/
		virtual void reset();

		/** Important: First call m_model->performNeighborhoodSearchSort()
		 * to call the z_sort of the neighborhood search.
		 */
		void performNeighborhoodSearchSort();

		void emittedParticles(FluidModel* model, const unsigned int startIndex);

		FORCE_INLINE const Real getShi(const unsigned int fluidIndex, const unsigned int i) const
		{
			return m_Shi[fluidIndex][i];
		}

		FORCE_INLINE Real& getShi(const unsigned int fluidIndex, const unsigned int i)
		{
			return m_Shi[fluidIndex][i];
		}

		FORCE_INLINE void setShi(const unsigned int fluidIndex, const unsigned int i, const Real sh_i)
		{
			m_Shi[fluidIndex][i] = sh_i;
		}

		FORCE_INLINE const Real getCorrected(const unsigned int fluidIndex, const unsigned int i) const
		{
			return m_corrected[fluidIndex][i];
		}

		FORCE_INLINE Real& getCorrected(const unsigned int fluidIndex, const unsigned int i)
		{
			return m_corrected[fluidIndex][i];
		}

		FORCE_INLINE void setCorrected(const unsigned int fluidIndex, const unsigned int i, const Real correction)
		{
			m_corrected[fluidIndex][i] = correction;
		}

		FORCE_INLINE const Real getPredictedDensity(const unsigned int fluidIndex, const unsigned int i) const
		{
			return m_predRho[fluidIndex][i];
		}

		FORCE_INLINE Real& getPredictedDensity(const unsigned int fluidIndex, const unsigned int i)
		{
			return m_predRho[fluidIndex][i];
		}

		FORCE_INLINE void setPredictedDensity(const unsigned int fluidIndex, const unsigned int i, const Real prediction)
		{
			m_predRho[fluidIndex][i] = prediction;
		}

		FORCE_INLINE const Real getPredictedDensityError(const unsigned int fluidIndex, const unsigned int i) const
		{
			return m_errPredRho[fluidIndex][i];
		}

		FORCE_INLINE Real& getPredictedDensityError(const unsigned int fluidIndex, const unsigned int i)
		{
			return m_errPredRho[fluidIndex][i];
		}

		FORCE_INLINE void setPredictedDensityError(const unsigned int fluidIndex, const unsigned int i, const Real error)
		{
			m_errPredRho[fluidIndex][i] = error;
		}

		FORCE_INLINE const Real getK(const unsigned int fluidIndex) const
		{
			return m_K[fluidIndex];
		}

		FORCE_INLINE Real& getK(const unsigned int fluidIndex)
		{
			return m_K[fluidIndex];
		}

		FORCE_INLINE void setK(const unsigned int fluidIndex, const Real k)
		{
			m_K[fluidIndex] = k;
		}

		FORCE_INLINE const Real getPZero(const unsigned int fluidIndex, const unsigned int i) const
		{
			return m_pZero[fluidIndex][i];
		}

		FORCE_INLINE Real& getPZero(const unsigned int fluidIndex, const unsigned int i)
		{
			return m_pZero[fluidIndex][i];
		}

		FORCE_INLINE void setPZero(const unsigned int fluidIndex, const unsigned int i, const Real p)
		{
			m_pZero[fluidIndex][i] = p;
		}

		FORCE_INLINE const Real getPamb() const
		{
			return m_pAmb;
		}

		FORCE_INLINE Real& getPamb()
		{
			return m_pAmb;
		}

		FORCE_INLINE void setPamb(const Real p)
		{
			m_pAmb = p;
		}

		FORCE_INLINE const Real getAii(const unsigned int fluidIndex, const unsigned int i) const
		{
			return m_aii[fluidIndex][i];
		}

		FORCE_INLINE Real& getAii(const unsigned int fluidIndex, const unsigned int i)
		{
			return m_aii[fluidIndex][i];
		}

		FORCE_INLINE void setAii(const unsigned int fluidIndex, const unsigned int i, const Real aii)
		{
			m_aii[fluidIndex][i] = aii;
		}

		FORCE_INLINE const Real getPressure(const unsigned int fluidIndex, const unsigned int i) const
		{
			return m_pressure[fluidIndex][i];
		}

		FORCE_INLINE Real& getPressure(const unsigned int fluidIndex, const unsigned int i)
		{
			return m_pressure[fluidIndex][i];
		}

		FORCE_INLINE void setPressure(const unsigned int fluidIndex, const unsigned int i, const Real p)
		{
			m_pressure[fluidIndex][i] = p;
		}

		FORCE_INLINE const Real getLastPressure(const unsigned int fluidIndex, const unsigned int i) const
		{
			return m_lastPressure[fluidIndex][i];
		}

		FORCE_INLINE Real& getLastPressure(const unsigned int fluidIndex, const unsigned int i)
		{
			return m_lastPressure[fluidIndex][i];
		}

		FORCE_INLINE void setLastPressure(const unsigned int fluidIndex, const unsigned int i, const Real p)
		{
			m_lastPressure[fluidIndex][i] = p;
		}

		FORCE_INLINE Vector3r& getPressureAccel(const unsigned int fluidIndex, const unsigned int i)
		{
			return m_pressureAccel[fluidIndex][i];
		}

		FORCE_INLINE const Vector3r& getPressureAccel(const unsigned int fluidIndex, const unsigned int i) const
		{
			return m_pressureAccel[fluidIndex][i];
		}

		FORCE_INLINE void setPressureAccel(const unsigned int fluidIndex, const unsigned int i, const Vector3r& val)
		{
			m_pressureAccel[fluidIndex][i] = val;
		}

		FORCE_INLINE Vector3r& getPrevPos(const unsigned int fluidIndex, const unsigned int i)
		{
			return m_prevPos[fluidIndex][i];
		}

		FORCE_INLINE const Vector3r& getPrevPos(const unsigned int fluidIndex, const unsigned int i) const
		{
			return m_prevPos[fluidIndex][i];
		}

		FORCE_INLINE void setPrevPos(const unsigned int fluidIndex, const unsigned int i, const Vector3r& val)
		{
			m_prevPos[fluidIndex][i] = val;
		}

		FORCE_INLINE const Real getSi(const unsigned int fluidIndex, const unsigned int i) const
		{
			return m_si[fluidIndex][i];
		}

		FORCE_INLINE Real& getSi(const unsigned int fluidIndex, const unsigned int i)
		{
			return m_si[fluidIndex][i];
		}

		FORCE_INLINE void setSi(const unsigned int fluidIndex, const unsigned int i, const Real si)
		{
			m_si[fluidIndex][i] = si;
		}

		FORCE_INLINE const Real getAijPj(const unsigned int fluidIndex, const unsigned int i) const
		{
			return m_aij_pj[fluidIndex][i];
		}

		FORCE_INLINE Real& getAijPj(const unsigned int fluidIndex, const unsigned int i)
		{
			return m_aij_pj[fluidIndex][i];
		}

		FORCE_INLINE void setAijPj(const unsigned int fluidIndex, const unsigned int i, const Real aij_pj)
		{
			m_aij_pj[fluidIndex][i] = aij_pj;
		}

	};
}

#endif