#include "SurfaceTension_Tartakovsky2016.h"
#include <iostream>
#include "../Simulation.h"
#include "SPlisHSPlasH/BoundaryModel_Akinci2012.h"
#include "SPlisHSPlasH/BoundaryModel_Koschier2017.h"
#include "SPlisHSPlasH/BoundaryModel_Bender2019.h"

using namespace SPH;

SurfaceTension_Tartakovsky2016::SurfaceTension_Tartakovsky2016(FluidModel *model) :
	SurfaceTensionBase(model)
{
}

SurfaceTension_Tartakovsky2016::~SurfaceTension_Tartakovsky2016(void)
{
}

void SurfaceTension_Tartakovsky2016::step()
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int numParticles = m_model->numActiveParticles();
	const Real k = m_surfaceTension;
	const Real kb = m_surfaceTensionBoundary;
	const Real h = sim->getSupportRadius() / 2.0;

	const unsigned int fluidModelIndex = m_model->getPointSetIndex();
	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int nBoundaries = sim->numberOfBoundaryModels();
	const Real density0 = m_model->getDensity0();
	FluidModel *model = m_model;

	// Compute forces
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			const Vector3r &xi = m_model->getPosition(i);
			Vector3r &ai = m_model->getAcceleration(i);
			const Real& mi = m_model->getMass(i);

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			forall_fluid_neighbors_in_same_phase(

				Real r_ij = (xi - xj).norm();
				ai -= k / mi * cos(3.0 * M_PI / (4.0 * h) * r_ij) / r_ij * (xj - xi);
			);

			//////////////////////////////////////////////////////////////////////////
			// Boundary TODO
			//////////////////////////////////////////////////////////////////////////
			if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
			{
				forall_boundary_neighbors(

					Real r_ij = (xi - xj).norm();
					
					ai -= kb / mi * cos(3.0 * M_PI / (4.0 * h) * r_ij) / r_ij * (xj - xi);

					bm_neighbor->addForce(xj, -model->getMass(i) * ai);
				);
			}
		}
	}
}


void SurfaceTension_Tartakovsky2016::reset()
{
}

