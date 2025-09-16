#include "SurfaceTension_Li2018.h"
#include <iostream>
#include "../Simulation.h"
#include "SPlisHSPlasH/BoundaryModel_Akinci2012.h"
#include "SPlisHSPlasH/BoundaryModel_Koschier2017.h"
#include "SPlisHSPlasH/BoundaryModel_Bender2019.h"

using namespace SPH;

int SurfaceTension_Li2018::A1 = -1;
int SurfaceTension_Li2018::A2 = -1;
int SurfaceTension_Li2018::H1 = -1;
int SurfaceTension_Li2018::H2 = -1;

SurfaceTension_Li2018::SurfaceTension_Li2018(FluidModel* model) :
	SurfaceTensionBase(model)
{
	// TODO: Check different parameters for boundary force
	m_a1 = 0.0;
	m_a2 = 0.0;
	m_h1 = 0.0;
	m_h2 = 0.0;

	// P_virial

	m_virialPG = 0.0;
	m_virialPL.resize(model->numParticles(), 0.0);
	m_tensors.resize(model->numParticles(), Matrix3r::Zero());

	model->addField({ "virialP_Global", FieldType::Scalar, [&](const unsigned int i) -> Real* { return &m_virialPG; }, false });
	model->addField({ "virialP_Local", FieldType::Scalar, [&](const unsigned int i) -> Real* { return &m_virialPL[i]; } });
}

SurfaceTension_Li2018::~SurfaceTension_Li2018(void)
{
	delete m_neighborhoodSearch;
	m_model->removeFieldByName("virialP_Global");
	m_model->removeFieldByName("virialP_Local");
}

void SurfaceTension_Li2018::initParameters()
{
	SurfaceTensionBase::initParameters();

	A1 = createNumericParameter("A1", "Strength of Repulsive Force", &m_a1);
	setGroup(A1, "Fluid Model|Surface tension");
	setDescription(A1, "Change the Strength of Repulsive Force.");

	A2 = createNumericParameter("A2", "Strength of Attractive Force", &m_a2);
	setGroup(A2, "Fluid Model|Surface tension");
	setDescription(A2, "Change the Strength of Attractive Force.");

	H1 = createNumericParameter("H1", "Support Multiplier for Repulsive Force", &m_h1);
	setGroup(H1, "Fluid Model|Surface tension");
	setDescription(H1, "Change the Support Multiplier for Repulsive Force.");

	H2 = createNumericParameter("H2", "Support Multiplier for Attractive Force", &m_h2);
	setGroup(H2, "Fluid Model|Surface tension");
	setDescription(H2, "Change the Support Multiplier for Attractive Force.");
}

void SurfaceTension_Li2018::deferredInit()
{
	// TODO: Add optimization when STSupportRadius == sim->getSupportRadius()

	Simulation* sim = Simulation::getCurrent();

	Real biggerMultiplier = (m_h1 >= m_h2) ? m_h1 : m_h2;
	Real STSupportRadius = 3.0 * biggerMultiplier * sim->getSupportRadius() / 2.0; // 3 from gaussian, support = 3/2 h

	m_kernel1.setRadius(m_h1 * sim->getSupportRadius());
	m_kernel2.setRadius(m_h2 * sim->getSupportRadius());
	m_boundaryKernel.setRadius(sim->getSupportRadius() / 2.0); // h_bar = h for Wendland

	if (m_neighborhoodSearch == NULL)
#ifdef GPU_NEIGHBORHOOD_SEARCH
		m_neighborhoodSearch = new NeighborhoodSearch(STSupportRadius);
#else
		m_neighborhoodSearch = new NeighborhoodSearch(STSupportRadius, false);
#endif

	m_neighborhoodSearch->set_radius(STSupportRadius);

	for (int i = 0; i < (int)sim->numberOfFluidModels(); i++)
	{
		FluidModel* fm = sim->getFluidModel(i);
		m_neighborhoodSearch->add_point_set(&fm->getPosition(0)[0], fm->numActiveParticles(), true, true, true, fm);
	}

	for (int i = 0; i < (int)sim->numberOfBoundaryModels(); i++)
	{
		if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
		{
			BoundaryModel_Akinci2012* bm = static_cast<BoundaryModel_Akinci2012*>(sim->getBoundaryModel(i));
			m_neighborhoodSearch->add_point_set(&bm->getPosition(0)[0], bm->numberOfParticles(), false, true, true, bm);
		}
	}
}

void SurfaceTension_Li2018::step()
{
	Simulation* sim = Simulation::getCurrent();
	const unsigned int numParticles = m_model->numActiveParticles();
	const Real k = m_surfaceTension;
	const Real kb = m_surfaceTensionBoundary;
	
	const unsigned int fluidModelIndex = m_model->getPointSetIndex();
	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int nBoundaries = sim->numberOfBoundaryModels();
	const Real density0 = m_model->getDensity0();
	FluidModel* model = m_model;

	Real biggerMultiplier = (m_h1 >= m_h2) ? m_h1 : m_h2;
	Real h_bar = 3.0 * biggerMultiplier * sim->getSupportRadius() / 2.0;

	// Neighborhood Search with bigger supportRadius

	m_neighborhoodSearch->find_neighbors();

	// Reset global virial pressure

	m_virialPG = 0.0;

	// Compute forces
#pragma omp parallel default(shared)
	{
#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			const Vector3r& xi = m_model->getPosition(i);
			const Real& mi = m_model->getMass(i);
			Vector3r& ai = m_model->getAcceleration(i);
			Real& virialPi = getVirialP(i);
			Matrix3r& tensor = getTensor(i);

			// Reset local virial pressure
			virialPi = 0.0;
			tensor = Matrix3r::Zero();

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////

			const unsigned int numberOfFluidNeighbors = static_cast<unsigned int>(m_neighborhoodSearch->point_set(fluidModelIndex).n_neighbors(fluidModelIndex, i));

			for (unsigned int j = 0; j < numberOfFluidNeighbors; j++)
			{
				const unsigned int neighborIndex = m_neighborhoodSearch->point_set(fluidModelIndex).neighbor(fluidModelIndex, i, j);
				const Vector3r& xj = model->getPosition(neighborIndex);

				const Real& mj = m_model->getMass(neighborIndex);
				Vector3r xixj = (xi - xj);
				Real r_ij = (xi - xj).norm();
				// Maybe I should make sure xixj!=0

				Vector3r a_ij = -1.0 * k / mi * (m_a1 * m_kernel1.gradW(xixj) - m_a2 * m_kernel2.gradW(xixj));

				ai += a_ij;

				Real h_particle = sim->getSupportRadius() / 2.0;
				m_virialPG += mi * (xi - xj).dot(a_ij) / model->numParticles() / (h_particle * h_particle * h_particle) / 6.0; // Nh^3 \approx Vr, 2d = 6. 
				//virialPi += 0.5 * mi * (xi - xj).dot(accel) / (h_particle * h_particle * h_particle);

				// Virial Pressure geometric model

				Real volume = 4.0 / 3.0 * M_PI * h_bar * h_bar * h_bar;

				tensor -= 0.5 / volume * mi * a_ij * (xi - xj).transpose(); // Self-term

				const unsigned int numberOfFluidNeighbors_j = static_cast<unsigned int>(m_neighborhoodSearch->point_set(fluidModelIndex).n_neighbors(fluidModelIndex, j));

				for (unsigned int k = 0; k < numberOfFluidNeighbors_j; k++)
				{
					const unsigned int secondNeighborIndex = m_neighborhoodSearch->point_set(fluidModelIndex).neighbor(fluidModelIndex, j, k);
					const Vector3r& xk = model->getPosition(secondNeighborIndex);
					Real r_jk = (xj - xk).norm();
					Real r_ik = (xi - xk).norm();

					// Integral calculation

					Real d = r_jk;
					if (r_ik > h_bar)
					{
						Real cos_j = (xk - xj).dot(xi - xj) / r_jk / r_ij;

						d = r_ij * (cos_j + sqrtf(cos_j * cos_j + h_bar * h_bar / r_ij / r_ij - 1));

					}

					Vector3r f_jk = -1.0 * k * (m_a1 * m_kernel1.gradW(xj-xk) - m_a2 * m_kernel2.gradW(xj-xk));;

					tensor -= 0.5 / volume * d / r_jk * f_jk * (xj - xk).transpose();
				}
			}

			virialPi = -1.0 / 3.0 * tensor.trace();

			//////////////////////////////////////////////////////////////////////////
			// Boundary TODO (other boundary handling)
			//////////////////////////////////////////////////////////////////////////

			if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
			{
				for (unsigned int pid = nFluids; pid < static_cast<unsigned int>(m_neighborhoodSearch->n_point_sets()); pid++)
				{
					BoundaryModel_Akinci2012* bm_neighbor = static_cast<BoundaryModel_Akinci2012*>(m_neighborhoodSearch->point_set(pid).get_user_data());
					const unsigned int numberOfBoundaryNeighbors = static_cast<unsigned int>(m_neighborhoodSearch->point_set(fluidModelIndex).n_neighbors(pid, i));

					for (unsigned int j = 0; j < numberOfBoundaryNeighbors; j++)
					{
						const unsigned int neighborIndex = m_neighborhoodSearch->point_set(fluidModelIndex).neighbor(pid, i, j);

						const Vector3r& xj = bm_neighbor->getPosition(neighborIndex);

						Real r_ij = (xi - xj).norm();

						ai += -1.0 * kb / mi * (m_a1 * m_kernel1.gradW(xi-xj) - m_a2 * m_kernel2.gradW(xi-xj));

						// Monaghan & Kajtar boundary force

						if (r_ij <= m_boundaryKernel.getRadius() && r_ij > 1.0e-9)
						{
							// Units?
							ai += kb / mi * m_boundaryKernel.W(r_ij) * (xi - xj) / (r_ij * r_ij);
						}

						bm_neighbor->addForce(xj, -model->getMass(i) * ai);
					}
				}
			}
		}
	}
}


void SurfaceTension_Li2018::reset()
{
}
