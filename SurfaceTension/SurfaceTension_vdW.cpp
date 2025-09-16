#include "SurfaceTension_vdW.h"
#include <iostream>
#include "../Simulation.h"
#include "SPlisHSPlasH/BoundaryModel_Akinci2012.h"
#include "SPlisHSPlasH/BoundaryModel_Koschier2017.h"
#include "SPlisHSPlasH/BoundaryModel_Bender2019.h"

using namespace SPH;

SurfaceTension_vdW::SurfaceTension_vdW(FluidModel* model) :
	SurfaceTensionBase(model)
{
	Simulation* sim = Simulation::getCurrent();
	Real supportRadius = static_cast<Real>(4.0) * sim->getSupportRadius();

	if (m_neighborhoodSearch == NULL)
#ifdef GPU_NEIGHBORHOOD_SEARCH
		m_neighborhoodSearch = new NeighborhoodSearch(supportRadius);
#else
		m_neighborhoodSearch = new NeighborhoodSearch(supportRadius, false);
#endif
	m_neighborhoodSearch->set_radius(supportRadius);

	for (int i = 0; i < (int)sim->numberOfFluidModels(); i++)
	{
		FluidModel* fm = sim->getFluidModel(i);
		m_neighborhoodSearch->add_point_set(&fm->getPosition(0)[0], fm->numActiveParticles(), true, true, true, fm);
	}

	/*model->addField({ "normal", FieldType::Vector3, [&](const unsigned int i) -> Real* { return &m_normals[i][0]; }, false });*/
}

SurfaceTension_vdW::~SurfaceTension_vdW(void)
{
	/*m_model->removeFieldByName("normal");*/
	delete m_neighborhoodSearch;
}

void SurfaceTension_vdW::step()
{
	Simulation* sim = Simulation::getCurrent();
	const Real density0 = m_model->getDensity0();
	const unsigned int numParticles = m_model->numActiveParticles();
	const Real k = m_surfaceTension;
	const Real kb = m_surfaceTensionBoundary;
	const unsigned int fluidModelIndex = m_model->getPointSetIndex();
	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int nBoundaries = sim->numberOfBoundaryModels();
	FluidModel* model = m_model;

	// Neighborhood Search with bigger supportRadius

	m_neighborhoodSearch->find_neighbors();

	// Compute forces
#pragma omp parallel default(shared)
	{
#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			const Vector3r& xi = m_model->getPosition(i);
			const Real& rhoi = m_model->getDensity(i);
			const Real& mi = m_model->getMass(i);
			const Real a_bar = k / (mi * mi);
			Vector3r& ai = m_model->getAcceleration(i);

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////

			// Literally copypasted from forall_fluid_neighbors_in_same_phase macro
			for (unsigned int j = 0; j < numberOfNeighbors(fluidModelIndex, fluidModelIndex, i); j++)
			{
				const unsigned int neighborIndex = getNeighbor(fluidModelIndex, fluidModelIndex, i, j);
				const Vector3r& xj = model->getPosition(neighborIndex);
				const Real& rhoj = m_model->getDensity(neighborIndex);
				const Real& mj = m_model->getMass(neighborIndex);

				Vector3r xixj = (xi - xj);

				ai += 2.0 * a_bar * mj * CubicKernelvdW::gradW(xixj);
			}
		}
	}
}


void SurfaceTension_vdW::reset()
{
}

unsigned int SurfaceTension_vdW::getNeighbor(const unsigned int pointSetIndex, const unsigned int neighborPointSetIndex, const unsigned int index, const unsigned int k) const
{
	return m_neighborhoodSearch->point_set(pointSetIndex).neighbor(neighborPointSetIndex, index, k);
}

unsigned int SurfaceTension_vdW::numberOfNeighbors(const unsigned int pointSetIndex, const unsigned int neighborPointSetIndex, const unsigned int index) const
{
	return static_cast<unsigned int>(m_neighborhoodSearch->point_set(pointSetIndex).n_neighbors(neighborPointSetIndex, index));
}