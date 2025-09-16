#include "TimeStepVdWSPH.h"
#include "SPlisHSPlasH/TimeManager.h"
#include "SPlisHSPlasH/SPHKernels.h"
#include "SimulationDataVdWSPH.h"
#include <iostream>
#include "Utilities/Timing.h"
#include "SPlisHSPlasH/Simulation.h"
#include "SPlisHSPlasH/BoundaryModel_Akinci2012.h"
#include "SPlisHSPlasH/BoundaryModel_Koschier2017.h"
#include "SPlisHSPlasH/BoundaryModel_Bender2019.h"


using namespace SPH;
using namespace std;
using namespace GenParam;

//GUI params
int TimeStepVdWSPH::B_PARAM = -1;
int TimeStepVdWSPH::K_PARAM = -1;
int TimeStepVdWSPH::T = -1;


TimeStepVdWSPH::TimeStepVdWSPH() :
	TimeStep()
{
	m_simulationData.init();

	m_bParameter = 0.5;
	m_kParameter = 1.0;
	m_T = 0.4;

	Simulation* sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
	{
		FluidModel* model = sim->getFluidModel(fluidModelIndex);
		model->addField({ "pressure", FieldType::Scalar, [this, fluidModelIndex](const unsigned int i) -> Real* { return &m_simulationData.getPressure(fluidModelIndex, i); } });
		model->addField({ "pressure acceleration", FieldType::Vector3, [this, fluidModelIndex](const unsigned int i) -> Real* { return &m_simulationData.getPressureAccel(fluidModelIndex, i)[0]; } });
	}
}

TimeStepVdWSPH::~TimeStepVdWSPH(void)
{
	Simulation* sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
	{
		FluidModel* model = sim->getFluidModel(fluidModelIndex);
		model->removeFieldByName("pressure");
		model->removeFieldByName("pressure acceleration");
	}
}

void TimeStepVdWSPH::initParameters()
{
	TimeStep::initParameters();

	B_PARAM = createNumericParameter("size", "b_bar", &m_bParameter);
	setGroup(B_PARAM, "Simulation|VdWSPH");
	setDescription(B_PARAM, "Size of the Particle.");
	static_cast<RealParameter*>(getParameter(B_PARAM))->setMinValue(static_cast<Real>(1e-6));

	K_PARAM = createNumericParameter("boltzmann", "Boltzmann", &m_kParameter);
	setGroup(K_PARAM, "Simulation|VdWSPH");
	setDescription(K_PARAM, "Reduced Boltzmann Constant.");
	static_cast<RealParameter*>(getParameter(K_PARAM))->setMinValue(static_cast<Real>(1e-6));

	T = createNumericParameter("temperature", "Temperature", &m_T);
	setGroup(T, "Simulation|VdWSPH");
	setDescription(T, "Reduced Temperature.");
	static_cast<RealParameter*>(getParameter(T))->setMinValue(static_cast<Real>(1e-6));
}

void TimeStepVdWSPH::step()
{
	Simulation* sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();
	TimeManager* tm = TimeManager::getCurrent();
	const Real h = tm->getTimeStepSize();

	sim->performNeighborhoodSearch();

	// I have no idea what this does
#ifdef USE_PERFORMANCE_OPTIMIZATION
	precomputeValues();
#endif

	if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
		computeVolumeAndBoundaryX();
	else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
		computeDensityAndGradient();

	// Compute accelerations: a(t)
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
	{
		clearAccelerations(fluidModelIndex);
		computeDensities(fluidModelIndex);
	}
	sim->computeNonPressureForces();


	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
	{
		FluidModel* model = sim->getFluidModel(fluidModelIndex);
		//const Real density0 = model->getDensity0();
#pragma omp parallel default(shared)
		{
#pragma omp for schedule(static)  
			for (int i = 0; i < (int)model->numActiveParticles(); i++)
			{
				Real& density = model->getDensity(i);
				//density = max(density, density0);
				m_simulationData.getPressure(fluidModelIndex, i) = density * m_kParameter * m_T / (1.0 - density * m_bParameter);
			}
		}

		computePressureAccels(fluidModelIndex);
	}

	sim->updateTimeStepSize();

	advectParticles();

	sim->emitParticles();
	sim->animateParticles();

	// Compute new time	
	tm->setTime(tm->getTime() + h);
}


void TimeStepVdWSPH::reset()
{
	TimeStep::reset();
	m_simulationData.reset();
}

void TimeStepVdWSPH::advectParticles()
{
	Simulation* sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();
	TimeManager* tm = TimeManager::getCurrent();
	const Real h = tm->getTimeStepSize();

	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
	{
		FluidModel* model = sim->getFluidModel(fluidModelIndex);
#pragma omp parallel default(shared)
		{
#pragma omp for schedule(static) 
			for (int i = 0; i < (int)model->numActiveParticles(); i++)
			{
				if (model->getParticleState(i) == ParticleState::Active)
				{
					Vector3r& pos = model->getPosition(i);
					Vector3r& vel = model->getVelocity(i);
					Vector3r& accel = model->getAcceleration(i);
					accel += m_simulationData.getPressureAccel(fluidModelIndex, i);
					vel += accel * h;
					pos += vel * h;
				}
			}
		}
	}
}

void TimeStepVdWSPH::computePressureAccels(const unsigned int fluidModelIndex)
{
	Simulation* sim = Simulation::getCurrent();
	FluidModel* model = sim->getFluidModel(fluidModelIndex);
	const Real density0 = model->getDensity0();
	const unsigned int numParticles = model->numActiveParticles();
	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int nBoundaries = sim->numberOfBoundaryModels();

	// Compute pressure forces
#pragma omp parallel default(shared)
	{
#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			const Vector3r& xi = model->getPosition(i);
			const Real density_i = model->getDensity(i);

			Vector3r& ai = m_simulationData.getPressureAccel(fluidModelIndex, i);
			ai.setZero();

			const Real dpi = m_simulationData.getPressure(fluidModelIndex, i) / (density_i * density_i);
			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			forall_fluid_neighbors(
				// Pressure 
				const Real density_j = fm_neighbor->getDensity(neighborIndex);
				const Real dpj = m_simulationData.getPressure(pid, neighborIndex) / (density_j * density_j);
				ai -= fm_neighbor->getMass(neighborIndex) * (dpi + dpj) * sim->gradW(xi - xj);
			);

			//////////////////////////////////////////////////////////////////////////
			// Boundary
			//////////////////////////////////////////////////////////////////////////
			const Real dpj = m_simulationData.getPressure(fluidModelIndex, i) / (density0 * density0);
			if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
			{
				forall_boundary_neighbors(
					const Vector3r a = density0 * bm_neighbor->getVolume(neighborIndex) * (dpi + dpj) * sim->gradW(xi - xj);
				ai -= a;
				bm_neighbor->addForce(xj, model->getMass(i) * a);
					);
			}
			else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
			{
				forall_density_maps(
					const Vector3r a = -density0 * (dpi + dpj) * gradRho;
				ai -= a;
				bm_neighbor->addForce(xj, model->getMass(i) * a);
					);
			}
			else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
			{
				forall_volume_maps(
					const Vector3r a = density0 * Vj * (dpi + dpj) * sim->gradW(xi - xj);
				ai -= a;
				bm_neighbor->addForce(xj, model->getMass(i) * a);
					);
			}
		}
	}
}

void TimeStepVdWSPH::performNeighborhoodSearchSort()
{
	m_simulationData.performNeighborhoodSearchSort();
}

void TimeStepVdWSPH::emittedParticles(FluidModel* model, const unsigned int startIndex)
{
	m_simulationData.emittedParticles(model, startIndex);
}

void TimeStepVdWSPH::resize()
{
	m_simulationData.init();
}

