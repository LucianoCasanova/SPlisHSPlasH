#include "SurfaceTension_Blank2024.h"
#include "SPlisHSPlasH/ISPH/TimeStepISPH.h"
#include <iostream>
#include "../Simulation.h"
#include "SPlisHSPlasH/BoundaryModel_Akinci2012.h"
#include "SPlisHSPlasH/BoundaryModel_Koschier2017.h"
#include "SPlisHSPlasH/BoundaryModel_Bender2019.h"
#include "Utilities/Logger.h"


using namespace SPH;

SurfaceTension_Blank2024::SurfaceTension_Blank2024(FluidModel* model) :
	SurfaceTensionBase(model)
{
	m_normals.resize(model->numParticles(), Vector3r::Zero());
	m_betterNormals.resize(model->numParticles(), Vector3r::Zero());
	m_curvatures.resize(model->numParticles(), 0.0);
	model->addField({ "normal", FieldType::Vector3, [&](const unsigned int i) -> Real* { return &m_normals[i][0]; }, false });
	model->addField({ "~^normal", FieldType::Vector3, [&](const unsigned int i) -> Real* { return &m_betterNormals[i][0]; }, false });
	model->addField({ "curvature", FieldType::Scalar, [&](const unsigned int i) -> Real* { return &m_curvatures[i]; }, false });
}

SurfaceTension_Blank2024::~SurfaceTension_Blank2024(void)
{
	m_normals.clear();
	m_betterNormals.clear();
	m_curvatures.clear();
	m_model->removeFieldByName("normal");
	m_model->removeFieldByName("~^normal");
	m_model->removeFieldByName("curvature");
}

void SurfaceTension_Blank2024::step()
{
	Simulation* sim = Simulation::getCurrent();
	const unsigned int numParticles = m_model->numActiveParticles();
	const Real k = m_surfaceTension;
	const Real kb = m_surfaceTensionBoundary;
	const unsigned int fluidModelIndex = m_model->getPointSetIndex();
	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int nBoundaries = sim->numberOfBoundaryModels();
	FluidModel* model = m_model;

	auto tsISPH = dynamic_cast<TimeStepISPH*>(sim->getTimeStep());

	if (tsISPH) // Only works with ISPH
	{
		// Compute pressures

		computeRawNormals();
		smoothAndCapNormals();

#pragma omp parallel default(shared)
		{
#pragma omp for schedule(static)  
			for (int i = 0; i < (int)numParticles; i++)
			{
				const Vector3r& xi = m_model->getPosition(i);
				const Real& mi = model->getMass(i);
				const Real& rho_i = model->getDensity(i);

				Real& kappa_i = m_curvatures[i];
				const Vector3r& n_i = m_betterNormals[i];
				Real& pZero_i = tsISPH->getSimulationData().getPZero(fluidModelIndex, i);

				//////////////////////////////////////////////////////////////////////////
				// Fluid
				//////////////////////////////////////////////////////////////////////////

				// Calculate Sh_i

				Real Sh_i = mi / rho_i * sim->W_zero();

				forall_fluid_neighbors(
					const Real & mj = fm_neighbor->getMass(neighborIndex);
					const Real & rho_j = fm_neighbor->getDensity(neighborIndex);

					Sh_i += mj / rho_j * sim->W(xi - xj);
				);

				// Compute gradient correction matrix

				Matrix3r L = computeGradCorrectionMatrix(i);

				// Curvature

				kappa_i = 0.0;

				forall_fluid_neighbors(
					const Real & mj = fm_neighbor->getMass(neighborIndex);
					const Real & rho_j = fm_neighbor->getDensity(neighborIndex);
					const Vector3r& n_j = m_betterNormals[neighborIndex];

					//Vector3r correctedGradW = sim->gradW(xi - xj);
					Vector3r correctedGradW = L * sim->gradW(xi - xj);

					kappa_i -= 0.5 * mj / rho_j * (n_i - n_j).dot(correctedGradW);
				);

				pZero_i += (1.0 + 1.0 / Sh_i) * k * kappa_i;

				//////////////////////////////////////////////////////////////////////////
				// Boundary (TODO)
				//////////////////////////////////////////////////////////////////////////

			}
		}
	}
}

void SurfaceTension_Blank2024::computeRawNormals()
{
	Simulation* sim = Simulation::getCurrent();
	const unsigned int numParticles = m_model->numActiveParticles();
	const unsigned int fluidModelIndex = m_model->getPointSetIndex();
	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int nBoundaries = sim->numberOfBoundaryModels();
	FluidModel* model = m_model;

#pragma omp parallel default(shared)
	{
#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			const Vector3r& xi = m_model->getPosition(i);
			Vector3r& n_i = m_normals[i];

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////

			n_i = Vector3r::Zero();

			forall_fluid_neighbors(
				const Real & mj = fm_neighbor->getMass(neighborIndex);
				const Real & rho_j = fm_neighbor->getDensity(neighborIndex);

				n_i -= mj / rho_j * sim->gradW(xi - xj);
			);
		}
	}
}

void SurfaceTension_Blank2024::smoothAndCapNormals()
{
	Simulation* sim = Simulation::getCurrent();
	const unsigned int numParticles = m_model->numActiveParticles();
	const unsigned int fluidModelIndex = m_model->getPointSetIndex();
	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int nBoundaries = sim->numberOfBoundaryModels();
	FluidModel* model = m_model;

	// Smooth, cap and normalize

#pragma omp parallel default(shared)
	{
#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			const Vector3r& xi = m_model->getPosition(i);
			const Real& mi = model->getMass(i);
			const Real& rho_i = model->getDensity(i);
			const Vector3r& n_i = m_normals[i];
			Vector3r& newN_i = m_betterNormals[i];

			const auto tsISPH = dynamic_cast<TimeStepISPH*>(sim->getTimeStep());

			// Smooth

			newN_i = mi / rho_i * n_i * sim->W_zero();

			forall_fluid_neighbors(
				const Real & mj = fm_neighbor->getMass(neighborIndex);
				const Real & rho_j = fm_neighbor->getDensity(neighborIndex);
				const Vector3r& n_j = m_normals[neighborIndex];

				newN_i += mj/rho_j * n_j * sim->W(xi-xj);
			);

			// Calculate Sh_i

			Real Sh_i = mi / rho_i * sim->W_zero();

			forall_fluid_neighbors(
				const Real & mj = fm_neighbor->getMass(neighborIndex);
				const Real & rho_j = fm_neighbor->getDensity(neighborIndex);

				Sh_i += mj / rho_j * sim->W(xi - xj);
			);

			newN_i /= Sh_i;

			// Cap

			const Real& h_smooth = sim->getSupportRadius() / 2.0;

			if (newN_i.norm() <= 0.2/h_smooth)
			{
				newN_i = Vector3r::Zero();
			}
			else	// Normalize
			{
				newN_i /= newN_i.norm();
			}
		}
	}
}

Matrix3r SurfaceTension_Blank2024::computeGradCorrectionMatrix(const unsigned int i)
{
	Simulation* sim = Simulation::getCurrent();
	const unsigned int fluidModelIndex = m_model->getPointSetIndex();
	const unsigned int nFluids = sim->numberOfFluidModels();

	Matrix3r A = Matrix3r::Zero();

	const Vector3r& xi = m_model->getPosition(i);

	forall_fluid_neighbors(
		const Real & mj = fm_neighbor->getMass(neighborIndex);
		const Real & rho_j = fm_neighbor->getDensity(neighborIndex);
		Vector3r x_ji = xj - xi;

		A += mj / rho_j * x_ji * sim->gradW(xi - xj).transpose();
	);

	return A.inverse();
}

void SurfaceTension_Blank2024::reset()
{
}