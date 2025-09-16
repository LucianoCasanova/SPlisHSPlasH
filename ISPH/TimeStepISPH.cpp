#include "TimeStepISPH.h"
#include "SPlisHSPlasH/TimeManager.h"
#include "SPlisHSPlasH/SPHKernels.h"
#include "SimulationDataISPH.h"
#include <iostream>
#include "Utilities/Timing.h"
#include "SPlisHSPlasH/Simulation.h"
#include "Utilities/Counting.h"
#include "SPlisHSPlasH/BoundaryModel_Akinci2012.h"
#include "SPlisHSPlasH/BoundaryModel_Koschier2017.h"
#include "SPlisHSPlasH/BoundaryModel_Bender2019.h"

//std::ofstream out("convergencia.csv");

using namespace SPH;
using namespace std;

int TimeStepISPH::USE_SURFACE_CORRECTION = -1;
int TimeStepISPH::P_AMB = -1;
int TimeStepISPH::OMEGA = -1;

TimeStepISPH::TimeStepISPH() :
	TimeStep()
{
	m_simulationData.init();
	m_useSurfaceCorrection = true;
	m_omega = 0.5;
	m_kObtained = false;

	Simulation* sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
	{
		FluidModel* model = sim->getFluidModel(fluidModelIndex);
		model->addField({ "Sh_i", FieldType::Scalar, [this, fluidModelIndex](const unsigned int i) -> Real* { return &m_simulationData.getShi(fluidModelIndex, i); } });
		model->addField({ "Indicator", FieldType::Scalar, [this, fluidModelIndex](const unsigned int i) -> Real* { return &m_simulationData.getCorrected(fluidModelIndex, i); } });
		model->addField({ "Predicted Density", FieldType::Scalar, [this, fluidModelIndex](const unsigned int i) -> Real* { return &m_simulationData.getPredictedDensity(fluidModelIndex, i); } });
		model->addField({ "Predicted Density Error", FieldType::Scalar, [this, fluidModelIndex](const unsigned int i) -> Real* { return &m_simulationData.getPredictedDensityError(fluidModelIndex, i); } });
		model->addField({ "Young-Laplace Pressure Jump", FieldType::Scalar, [this, fluidModelIndex](const unsigned int i) -> Real* { return &m_simulationData.getPZero(fluidModelIndex, i); } });
		model->addField({ "a_ii", FieldType::Scalar, [this, fluidModelIndex](const unsigned int i) -> Real* { return &m_simulationData.getAii(fluidModelIndex, i); } });
		model->addField({ "a_ij_pj", FieldType::Scalar, [this, fluidModelIndex](const unsigned int i) -> Real* { return &m_simulationData.getAijPj(fluidModelIndex, i); } });
		model->addField({ "pressure", FieldType::Scalar, [this, fluidModelIndex](const unsigned int i) -> Real* { return &m_simulationData.getPressure(fluidModelIndex, i); }, true });
		model->addField({ "last pressure", FieldType::Scalar, [this, fluidModelIndex](const unsigned int i) -> Real* { return &m_simulationData.getLastPressure(fluidModelIndex, i); } });
		model->addField({ "pressure acceleration", FieldType::Vector3, [this, fluidModelIndex](const unsigned int i) -> Real* { return &m_simulationData.getPressureAccel(fluidModelIndex, i)[0]; } });
		model->addField({ "s_i", FieldType::Scalar, [this, fluidModelIndex](const unsigned int i) -> Real* { return &m_simulationData.getSi(fluidModelIndex, i); }, true });
	}
}

TimeStepISPH::~TimeStepISPH(void)
{
	Simulation* sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
	{
		FluidModel* model = sim->getFluidModel(fluidModelIndex);
		model->removeFieldByName("Sh_i");
		model->removeFieldByName("Predicted Density");
		model->removeFieldByName("Predicted Density Error");
		model->removeFieldByName("Young-Laplace Pressure Jump");
		model->removeFieldByName("Indicator");
		model->removeFieldByName("a_ii");
		model->removeFieldByName("a_ij_pj");
		model->removeFieldByName("pressure");
		model->removeFieldByName("last pressure");
		model->removeFieldByName("pressure acceleration");
		model->removeFieldByName("s_i");
	}
	//out.close();
}

void TimeStepISPH::initParameters()
{
	TimeStep::initParameters();

	USE_SURFACE_CORRECTION = createBoolParameter("useSurfaceCorrection", "Enable surface correction (Nair2014)", &m_useSurfaceCorrection);
	setGroup(USE_SURFACE_CORRECTION, "Simulation|ISPH");
	setDescription(USE_SURFACE_CORRECTION, "Turn surface correction on/off.");

	OMEGA = createNumericParameter("omega", "Jacobi Relaxation", &m_omega);
	setGroup(OMEGA, "Simulation|ISPH");
	setDescription(OMEGA, "Set Jacobi relaxation coefficient.");

	P_AMB = createNumericParameter("pAmb", "External Pressure", &m_simulationData.getPamb());
	setGroup(P_AMB, "Simulation|ISPH");
	setDescription(P_AMB, "Set external pressure.");
}

void TimeStepISPH::step()
{
	Simulation* sim = Simulation::getCurrent();
	TimeManager* tm = TimeManager::getCurrent();
	const Real h = tm->getTimeStepSize();
	const unsigned int nModels = sim->numberOfFluidModels();

	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
		clearAccelerations(fluidModelIndex);

#ifdef USE_PERFORMANCE_OPTIMIZATION
	precomputeValues(); 
	// I am not sure if this can interfere, gonna activate once it is working
#endif

	// Using densities at x_i

	sim->performNeighborhoodSearch();

	if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
		computeVolumeAndBoundaryX();
	else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
		computeDensityAndGradient();

	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
	{
		computeDensities(fluidModelIndex);
	}

	// CFL condition at this step, because then nonPressureForces are calculated at r_adv = r_i + h * v_i
	// Anyhow this seems to change it for i+1 (nope, it is in pressureSolveIteration, why?)

	sim->updateTimeStepSize();

	// Calculates k for a bulk particle
	if (!m_kObtained && m_useSurfaceCorrection)
	{
		for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
			calculateK(fluidModelIndex);
	}

	// Calculate r_adv, v_adv, a_ii and s_i. Saves r_prev (for kernel evaluation)

	START_TIMING("predictAdvection");
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
		predictAdvection(fluidModelIndex);
	STOP_TIMING_AVG;

	// Solve PPE

	START_TIMING("pressureSolve");
	pressureSolve();
	STOP_TIMING_AVG;

	// Predictor-Corrector integrator

	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
		integration(fluidModelIndex);

	sim->emitParticles();
	sim->animateParticles();

	// Compute new time	
	tm->setTime(tm->getTime() + h);
}


void TimeStepISPH::reset()
{
	TimeStep::reset();
	m_simulationData.reset();
}

void TimeStepISPH::calculateK(const unsigned int fluidModelIndex)
{
	Simulation* sim = Simulation::getCurrent();
	FluidModel* model = sim->getFluidModel(fluidModelIndex);
	const unsigned int numParticles = model->numActiveParticles();
	if (numParticles == 0)
		return;

	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int nBoundaries = sim->numberOfBoundaryModels();

	Real& K = m_simulationData.getK(fluidModelIndex);
	K = 0.0;

	for (int i = 0; i < (int)numParticles && !m_kObtained; i++)
	{
		Real Shi = 0.0;
		const Vector3r& xi = model->getPosition(i);
		const Real& mi = model->getMass(i);
		const Real& rho_i = model->getDensity(i);
		Shi += mi / rho_i * sim->W_zero();

		forall_fluid_neighbors(
			const Real & mj = fm_neighbor->getMass(neighborIndex);
			const Real & rho_j = fm_neighbor->getDensity(neighborIndex);

			Shi += mj / rho_j * sim->W(xi - xj);
		);

		if (Shi > 0.95)
		{
			const Real& h_smooth = sim->getSupportRadius() / 2.0;
			forall_fluid_neighbors(
				const Real & mj = fm_neighbor->getMass(neighborIndex);
				const Real & rho_j = fm_neighbor->getDensity(neighborIndex);

				const Real Fij = (xi - xj).dot(sim->gradW(xi - xj)) / ((xi - xj).squaredNorm() + (0.01 * h_smooth)* (0.01 * h_smooth));
				K += (4.0 * mj * Fij) / (rho_j * (rho_i + rho_j));
			);
			m_kObtained = true;
		}
	}
}

void TimeStepISPH::predictAdvection(const unsigned int fluidModelIndex)
{
	Simulation* sim = Simulation::getCurrent();
	FluidModel* model = sim->getFluidModel(fluidModelIndex);
	const unsigned int numParticles = model->numActiveParticles();
	if (numParticles == 0)
		return;

	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int nBoundaries = sim->numberOfBoundaryModels();
	const Real density0 = model->getDensity0();

	const Real h = TimeManager::getCurrent()->getTimeStepSize();

	// Compute r_adv 

#pragma omp parallel default(shared)
	{
#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			Vector3r& xi = model->getPosition(i);
			const Vector3r& vi = model->getVelocity(i);

			Vector3r& xi_prev = m_simulationData.getPrevPos(fluidModelIndex, i);
			xi_prev = xi;

			if (model->getParticleState(i) == ParticleState::Active)
				xi += h * vi;
		}
	}

	// Initialize p_Zero with p_Amb

#pragma omp parallel default(shared)
	{
#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			Real& pZero_i = m_simulationData.getPZero(fluidModelIndex, i);
			pZero_i = m_simulationData.getPamb();
		}
	}

	// Compute v_adv, go back to PrevPos

	sim->performNeighborhoodSearch();
	sim->computeNonPressureForces();

#pragma omp parallel default(shared)
	{
#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			Vector3r& xi = model->getPosition(i);
			Vector3r& xi_prev = m_simulationData.getPrevPos(fluidModelIndex, i);

			xi = xi_prev;

			Vector3r& vi = model->getVelocity(i);
			const Vector3r& ai = model->getAcceleration(i);

			if (model->getParticleState(i) == ParticleState::Active)
				vi += h * ai;
		}
	}
	sim->performNeighborhoodSearch();

	// Compute Sh_i

	if (m_useSurfaceCorrection)
	{
#pragma omp parallel default(shared)
		{
#pragma omp for schedule(static)
			for (int i = 0; i < (int)numParticles; i++)
			{
				Real& Shi = m_simulationData.getShi(fluidModelIndex, i);
				Shi = 0.0;

				const Vector3r& xi = model->getPosition(i);
				const Real& mi = model->getMass(i);
				const Real& rho_i = model->getDensity(i);

				Shi += mi / rho_i * sim->W_zero();

				forall_fluid_neighbors(
					const Real & mj = fm_neighbor->getMass(neighborIndex);
					const Real & rho_j = fm_neighbor->getDensity(neighborIndex);

					Shi += mj / rho_j * sim->W(xi - xj);
					);
			}
		}
	}

	// Compute a_ii and s_i (with surface corrections), first guess of Pi

#pragma omp parallel default(shared)
	{
#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			Real& aii = m_simulationData.getAii(fluidModelIndex, i);
			Real& si = m_simulationData.getSi(fluidModelIndex, i);
			Real& corrected_i = m_simulationData.getCorrected(fluidModelIndex, i);

			aii = 0.0;
			si = 0.0;
			corrected_i = 0.0;

			const Vector3r& xi = model->getPosition(i);
			const Vector3r& vi = model->getVelocity(i);
			const Real& rho_i = model->getDensity(i);
			const Real& h_smooth = sim->getSupportRadius() / 2.0;

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			forall_fluid_neighbors(
				const Vector3r & vj = fm_neighbor->getVelocity(neighborIndex);
				const Real& mj = fm_neighbor->getMass(neighborIndex);
				const Real& rho_j = fm_neighbor->getDensity(neighborIndex);

				const Real Fij = (xi - xj).dot(sim->gradW(xi - xj)) / ((xi - xj).squaredNorm() + (0.01 * h_smooth)* (0.01 * h_smooth));
				aii += (4.0 * mj * Fij) / (rho_j * (rho_i + rho_j));
				si -= (mj * (vi - vj).dot(sim->gradW(xi - xj))) / (rho_j * h);
			);

			// Surface correction

			if (m_useSurfaceCorrection)
			{
				Real& Shi = m_simulationData.getShi(fluidModelIndex, i);

				if (Shi < 0.95)
				{
					corrected_i = 100.0;
					const Real& K = m_simulationData.getK(fluidModelIndex);
					const Real& pZero_i = m_simulationData.getPZero(fluidModelIndex, i);

					si += pZero_i * (K - aii);
					aii = K;
				}
			}

			//////////////////////////////////////////////////////////////////////////
			// Boundary TODO
			//////////////////////////////////////////////////////////////////////////
			/*if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
			{
				forall_boundary_neighbors(
					const Vector3r & vj = bm_neighbor->getVelocity(neighborIndex);
				densityAdv += h * bm_neighbor->getVolume(neighborIndex) * (vi - vj).dot(sim->gradW(xi - xj));
					);
			}
			else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
			{
				forall_density_maps(
					Vector3r vj;
				bm_neighbor->getPointVelocity(xi, vj);
				densityAdv -= h * (vi - vj).dot(gradRho);
					);
			}
			else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
			{
				forall_volume_maps(
					Vector3r vj;
				bm_neighbor->getPointVelocity(xj, vj);
				densityAdv += h * Vj * (vi - vj).dot(sim->gradW(xi - xj));
					);
			}*/

			// First guess for P^{0}_{n+1} = 1/2 P^{last}_{n}
			const Real& pressure = m_simulationData.getPressure(fluidModelIndex, i);
			Real& lastPressure = m_simulationData.getLastPressure(fluidModelIndex, i);
			lastPressure = static_cast<Real>(0.5) * pressure;
		}
	}
}

void TimeStepISPH::pressureSolve()
{
	Simulation* sim = Simulation::getCurrent();
	const unsigned int nFluids = sim->numberOfFluidModels();

	// Compute error in density prediction

	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++)
	{
		FluidModel* model = sim->getFluidModel(fluidModelIndex);
		const unsigned int numParticles = model->numActiveParticles();

#pragma omp parallel default(shared)
		{
#pragma omp for schedule(static)  
			for (int i = 0; i < (int)numParticles; i++)
			{
				const Real& rho_i = model->getDensity(i);
				const Real& rhoPred_i = m_simulationData.getPredictedDensity(fluidModelIndex, i);

				Real& rhoPredErr = m_simulationData.getPredictedDensityError(fluidModelIndex, i);

				rhoPredErr = abs(rho_i - rhoPred_i);
			}
		}
	}

	// Pressure Solve

	Real avg_density_err = 0;
	m_iterations = 0;
	bool chk = false;

	while ((!chk || (m_iterations < m_minIterations)) && (m_iterations < m_maxIterations))
	{
		chk = true;
		for (unsigned int i = 0; i < nFluids; i++)
		{
			FluidModel* model = sim->getFluidModel(i);
			const Real density0 = model->getDensity0();

			avg_density_err = 0.0;
			pressureSolveIteration(i, avg_density_err);

			// Maximal allowed density fluctuation
			const Real eta = m_maxError * static_cast<Real>(0.01) * density0;  // maxError is given in percent
			chk = chk && (fabs(avg_density_err) <= eta);
			bool condition = (!chk || (m_iterations < m_minIterations - 1)) && (m_iterations < m_maxIterations - 1);
			if (!condition)
			{
				LOG_INFO << TimeManager::getCurrent()->getTime() << "\t" << avg_density_err << "\t" << m_iterations+1;
			}	
		}

		m_iterations++;
	}
	INCREASE_COUNTER("ISPH - iterations", static_cast<Real>(m_iterations));
}

void TimeStepISPH::pressureSolveIteration(const unsigned int fluidModelIndex, Real& avg_density_err)
{
	Simulation* sim = Simulation::getCurrent();
	FluidModel* model = sim->getFluidModel(fluidModelIndex);
	const unsigned int numParticles = model->numActiveParticles();
	if (numParticles == 0)
		return;

	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int nBoundaries = sim->numberOfBoundaryModels();

	const Real h = TimeManager::getCurrent()->getTimeStepSize();

	// Compute aij_pj
#pragma omp parallel default(shared)
	{
#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			Real& aij_pj = m_simulationData.getAijPj(fluidModelIndex, i);
			aij_pj = 0.0;

			const Vector3r& xi = model->getPosition(i);
			const Real& rho_i = model->getDensity(i);
			const Real& h_smooth = sim->getSupportRadius() / 2.0;

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			forall_fluid_neighbors(
				const Real& mj = fm_neighbor->getMass(neighborIndex);
				const Real& rho_j = fm_neighbor->getDensity(neighborIndex);
				const Real Fij = (xi - xj).dot(sim->gradW(xi - xj)) / ((xi - xj).squaredNorm() + (0.01 * h_smooth)* (0.01 * h_smooth));
				aij_pj -= (4.0 * mj * Fij * m_simulationData.getLastPressure(fluidModelIndex, neighborIndex)) / (rho_j * (rho_i + rho_j));
			);
		}
	}

	// Compute new pressure
#pragma omp parallel default(shared)
	{
#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			Real& pi = m_simulationData.getPressure(fluidModelIndex, i);
			pi = 0.0;

			const Real& aii = m_simulationData.getAii(fluidModelIndex, i);
			const Real& si = m_simulationData.getSi(fluidModelIndex, i);
			const Real& lastPi = m_simulationData.getLastPressure(fluidModelIndex, i);
			const Real& aij_pj = m_simulationData.getAijPj(fluidModelIndex, i);
			const Real& rho_0 = model->getDensity0();
			const Real& rho_i = model->getDensity(i);
			Real& rhoPred_i = m_simulationData.getPredictedDensity(fluidModelIndex,i);

			if (fabs(aii) > 1.0e-9)
			{
				const Real pi_candidate = (static_cast<Real>(1.0) - m_omega) * lastPi + m_omega / aii * (si - aij_pj);
				pi = max(pi_candidate, static_cast<Real>(0.0));
			}
			else
				pi = 0.0;

			if (pi != 0.0)
			{
				
				/*const Real newRho = h * h * rho_i * (aii * pi + aij_pj - si) + rho_i;
				rhoPred_i = newRho;*/

//#pragma omp atomic
				//avg_density_err += newRho - rho_0;

				avg_density_err += (aii * pi + aij_pj - si) * (aii * pi + aij_pj - si);
			}
		}
	}

	for (int i = 0; i < (int)numParticles; i++)
	{
		const Real& pi = m_simulationData.getPressure(fluidModelIndex, i);
		Real& lastPi = m_simulationData.getLastPressure(fluidModelIndex, i);
		lastPi = pi;
	}

	//avg_density_err /= numParticles;
	//out << TimeManager::getCurrent()->getTime() << "," << m_iterations << "," << avg_density_err << "\n";
}

void TimeStepISPH::integration(const unsigned int fluidModelIndex)
{
	Simulation* sim = Simulation::getCurrent();
	FluidModel* model = sim->getFluidModel(fluidModelIndex);
	const unsigned int numParticles = model->numActiveParticles();
	if (numParticles == 0)
		return;

	// Compute pressure forces
	computePressureAccels(fluidModelIndex);

	Real h = TimeManager::getCurrent()->getTimeStepSize();

#pragma omp parallel default(shared)
	{
#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			if (model->getParticleState(i) == ParticleState::Active)
			{
				Vector3r& pos = model->getPosition(i);
				Vector3r& vel = model->getVelocity(i);
				Vector3r vel_prev = vel;
				vel += m_simulationData.getPressureAccel(fluidModelIndex, i) * h;
				pos += 0.5 * (vel_prev + vel) * h;
			}
		}
	}
}

void TimeStepISPH::computePressureAccels(const unsigned int fluidModelIndex)
{
	Simulation* sim = Simulation::getCurrent();
	FluidModel* model = sim->getFluidModel(fluidModelIndex);
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
			const Real& rho_i = model->getDensity(i);

			Vector3r& ai = m_simulationData.getPressureAccel(fluidModelIndex, i);
			ai.setZero();

			const Real rho_i2 = rho_i * rho_i;
			const Real& pi = m_simulationData.getPressure(fluidModelIndex, i);
			const Real dpi = pi / rho_i2;

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			if (m_simulationData.getShi(fluidModelIndex, i) >= 0.95)
			{
				forall_fluid_neighbors(
					const Real & rho_j = fm_neighbor->getDensity(neighborIndex);
					const Real rho_j2 = rho_j * rho_j;
					const Real & pj = m_simulationData.getPressure(pid, neighborIndex);
					const Real dpj = pj / rho_j2;
					const Real & mj = fm_neighbor->getMass(neighborIndex);
					ai -= mj * (dpi + dpj) * sim->gradW(xi - xj);
				);
			}
			else
			{
				forall_fluid_neighbors(
					const Real & rho_j = fm_neighbor->getDensity(neighborIndex);
					const Real & pj = m_simulationData.getPressure(pid, neighborIndex);
					const Real & mj = fm_neighbor->getMass(neighborIndex);
					const Real& pZero_i = m_simulationData.getPZero(fluidModelIndex, i);

					ai += mj * (pZero_i - pj) / (rho_j * rho_j) * sim->gradW(xi - xj);
				);
			}

			//////////////////////////////////////////////////////////////////////////
			// Boundary TODO
			//////////////////////////////////////////////////////////////////////////
			/*if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
			{
				forall_boundary_neighbors(
					const Vector3r a = bm_neighbor->getVolume(neighborIndex) * dpi * sim->gradW(xi - xj);
				ai -= a;
				bm_neighbor->addForce(xj, model->getMass(i) * a);
					);
			}
			else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
			{
				forall_density_maps(
					const Vector3r a = -dpi * gradRho;
				ai -= a;
				bm_neighbor->addForce(xj, model->getMass(i) * a);
					);
			}
			else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
			{
				forall_volume_maps(
					const Vector3r a = Vj * dpi * sim->gradW(xi - xj);
				ai -= a;
				bm_neighbor->addForce(xj, model->getMass(i) * a);
					);
			}*/
		}
	}
}

void TimeStepISPH::performNeighborhoodSearchSort()
{
	m_simulationData.performNeighborhoodSearchSort();
}

void TimeStepISPH::emittedParticles(FluidModel* model, const unsigned int startIndex)
{
	m_simulationData.emittedParticles(model, startIndex);
}

void TimeStepISPH::resize()
{
	m_simulationData.init();
}

