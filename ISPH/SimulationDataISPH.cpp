#include "SimulationDataISPH.h"
#include "SPlisHSPlasH/SPHKernels.h"
#include "SPlisHSPlasH/Simulation.h"

using namespace SPH;

SimulationDataISPH::SimulationDataISPH()
{
}

SimulationDataISPH::~SimulationDataISPH(void)
{
	cleanup();
}

void SimulationDataISPH::init()
{
	Simulation* sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();

	m_pAmb = 0.0;
	m_Shi.resize(nModels);
	m_corrected.resize(nModels);
	m_predRho.resize(nModels);
	m_errPredRho.resize(nModels);
	m_K.resize(nModels);
	m_pZero.resize(nModels);
	m_aii.resize(nModels);
	m_pressure.resize(nModels);
	m_lastPressure.resize(nModels);
	m_pressureAccel.resize(nModels);
	m_prevPos.resize(nModels);
	m_si.resize(nModels);
	m_aij_pj.resize(nModels);
	for (unsigned int i = 0; i < nModels; i++)
	{
		FluidModel* fm = sim->getFluidModel(i);
		m_Shi[i].resize(fm->numParticles(), 0.0);
		m_corrected[i].resize(fm->numParticles(), 0.0);
		m_predRho[i].resize(fm->numParticles(), 0.0);
		m_errPredRho[i].resize(fm->numParticles(), 0.0);
		m_pZero[i].resize(fm->numParticles(), m_pAmb);
		m_aii[i].resize(fm->numParticles(), 0.0);
		m_pressure[i].resize(fm->numParticles(), 0.0);
		m_lastPressure[i].resize(fm->numParticles(), 0.0);
		m_pressureAccel[i].resize(fm->numParticles(), Vector3r::Zero());
		m_prevPos[i].resize(fm->numParticles(), Vector3r::Zero());
		m_si[i].resize(fm->numParticles(), 0.0);
		m_aij_pj[i].resize(fm->numParticles(), 0.0);
	}
}

void SimulationDataISPH::cleanup()
{
	Simulation* sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();

	for (unsigned int i = 0; i < nModels; i++)
	{
		m_Shi[i].clear();
		m_corrected[i].clear();
		m_predRho[i].clear();
		m_errPredRho[i].clear();
		m_pZero[i].clear();
		m_aii[i].clear();
		m_pressure[i].clear();
		m_lastPressure[i].clear();
		m_pressureAccel[i].clear();
		m_prevPos[i].clear();
		m_si[i].clear();
		m_aij_pj[i].clear();
	}
	m_Shi.clear();
	m_corrected.clear();
	m_predRho.clear();
	m_errPredRho.clear();
	m_K.clear();
	m_pZero.clear();
	m_aii.clear();
	m_pressure.clear();
	m_lastPressure.clear();
	m_pressureAccel.clear();
	m_prevPos.clear();
	m_si.clear();
	m_aij_pj.clear();
}

void SimulationDataISPH::reset()
{
	Simulation* sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();

	for (unsigned int i = 0; i < nModels; i++)
	{
		FluidModel* fm = sim->getFluidModel(i);
		for (unsigned int j = 0; j < fm->numParticles(); j++)
		{
			m_lastPressure[i][j] = 0.0;
			m_pZero[i][j] = m_pAmb;
		}
	}
}

void SimulationDataISPH::performNeighborhoodSearchSort()
{
	Simulation* sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();

	for (unsigned int i = 0; i < nModels; i++)
	{
		FluidModel* fm = sim->getFluidModel(i);
		const unsigned int numPart = fm->numActiveParticles();
		if (numPart != 0)
		{
			auto const& d = sim->getNeighborhoodSearch()->point_set(fm->getPointSetIndex());
			d.sort_field(&m_Shi[i][0]);
			d.sort_field(&m_corrected[i][0]);
			d.sort_field(&m_predRho[i][0]);
			d.sort_field(&m_errPredRho[i][0]);
			d.sort_field(&m_pZero[i][0]);
			d.sort_field(&m_aii[i][0]);
			d.sort_field(&m_pressure[i][0]);
			d.sort_field(&m_lastPressure[i][0]);
			d.sort_field(&m_pressureAccel[i][0]);
			d.sort_field(&m_prevPos[i][0]);
			d.sort_field(&m_si[i][0]);
			d.sort_field(&m_aij_pj[i][0]);
		}
	}
}

void SimulationDataISPH::emittedParticles(FluidModel* model, const unsigned int startIndex)
{
	// initialize last pressure values for new particles
	const unsigned int fluidModelIndex = model->getPointSetIndex();
	for (unsigned int j = startIndex; j < model->numActiveParticles(); j++)
	{
		m_lastPressure[fluidModelIndex][j] = 0.0;
	}
}
