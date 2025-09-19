#include "FispactCompiledSource.hpp"

#include "mpi.h"
#include "openmc/constants.h"
#include "openmc/distribution.h"
#include "openmc/distribution_multi.h"
#include "openmc/mesh.h"
#include "openmc/message_passing.h"
#include "openmc/position.h"
#include "utils/PhotonSharingData.h"
#include <cstddef>
#include <cstdint>
#include <fcntl.h>
#include <iostream>
#include <string>
#include <sys/types.h>

FizzyCompiledSource::FizzyCompiledSource(int32_t mesh_id)
    : _mesh_id(mesh_id), angle_(openmc::UPtrAngle(new openmc::Isotropic)) {
  // Set up time_
  double T[]{0.0};
  double p[]{1.0};
  time_ = openmc::UPtrDist{new openmc::Discrete{T, p, 1}};

  // Get shared interprocess data
  const std::string shared_data_name = generateInterprocessName();
  try {
    segment_ = boost::interprocess::managed_shared_memory(
        boost::interprocess::open_only, shared_data_name.c_str());
  } catch (bi::interprocess_exception) {
    std::cerr << "Could not find interprocess segment with name "
              << shared_data_name << std::endl;

    exit(-1);
  }

  // Sets strength for entire source term, not just this local contribution
  // Mostly used for tally normalisation later on!
  strength_ = 1;
}

const std::string FizzyCompiledSource::generateInterprocessName() {
  char mpi_proc_name[MPI_MAX_PROCESSOR_NAME];
  int len = 0;
  int err = MPI_Get_processor_name(mpi_proc_name, &len);

  std::string ipc_name = std::string(mpi_proc_name);
  ipc_name += "_" + std::to_string(openmc::mpi::rank);

  return ipc_name;
}

FizzyCompiledSource::~FizzyCompiledSource() {
  const std::string shared_data_name = generateInterprocessName();
  bi::shared_memory_object::remove(shared_data_name.c_str());
}

double FizzyCompiledSource::calculateParticleWeight(
    const double &local_domain_strength,
    const double &total_domain_strength) const {
  return local_domain_strength / total_domain_strength;
}

double FizzyCompiledSource::calculateParticleWeight(
    const PhotonSharingData *shared_data) const {
  return calculateParticleWeight(shared_data->_local_domain_strength,
                                 shared_data->_total_domain_strength);
}

void FizzyCompiledSource::setupLocalElementsDiscreteIndex(
    const PhotonSharingData *shared_data) {

  const BoostIpIntDoubMap &elem_id_to_strength_map =
      shared_data->_elem_strength;

  std::vector<double> element_strengths(elem_id_to_strength_map.size(), 0);

  element_ids_.assign(elem_id_to_strength_map.size(), 0);

  int i = 0;
  for (auto &pair : elem_id_to_strength_map) {
    element_ids_.at(i) = pair.first;
    element_strengths.at(i) = pair.second;
    i++;
  }

  // Set up discrete index to sample element id's from, copying behavoir from
  // openmc::MeshSource
  di_.assign(element_strengths);
}

int32_t FizzyCompiledSource::sampleLocalElementsIndex(uint64_t *seed) const {

  int index = di_.sample(seed);
  return element_ids_.at(index);
}

openmc::Position
FizzyCompiledSource::sampleElementVolume(uint64_t *seed, int32_t mesh_id,
                                         int element_id) const {
  // Get openmc mesh
  openmc::Position r;

  int32_t mesh_idx = openmc::model::mesh_map.at(mesh_id);
  std::unique_ptr<openmc::Mesh> &mesh = openmc::model::meshes[mesh_idx];

  openmc::LibMesh *derived_libmesh_ptr =
      dynamic_cast<openmc::LibMesh *>(mesh.get());

  do {
    r = derived_libmesh_ptr->sample_element(element_id, seed);
  } while (!this->satisfies_spatial_constraints(r));
  return r;
}

double
FizzyCompiledSource::sampleElementEnergy(uint64_t *seed,
                                         const PhotonSharingData *shared_data,
                                         const int32_t &element_id) const {

  // Set up ptrs to element energy distribution and photon bins in shared data,
  // this is a bit messy but avoids doing a copy!
  const double *element_energy =
      &shared_data->_photon_fluxes.at(element_id).at(0);
  const double *photon_bins = &shared_data->_photon_bins.at(0);
  size_t n_bins = shared_data->_photon_bins.size();
  openmc::Tabular distribution(photon_bins, element_energy, n_bins,
                               openmc::Interpolation::histogram, nullptr);

  while (true) {
    double energy = distribution.sample(seed);

    if (satisfies_energy_constraints(energy)) {
      return energy;
    }
  }
}

openmc::SourceSite FizzyCompiledSource::sample(uint64_t *seed) const {

  // If this is the first sample, do setup.
  // This setup should really be in the constructor, but doing it here makes it
  // compatible with MOOSE-multiapp runs. When using as a MOOSE-multiapp, the
  // constructor gets called before the shared data is ready to be read
  if (!setup_) {
    auto *p_this = const_cast<FizzyCompiledSource *>(this);
    p_this->instance_ = p_this->segment_.find<PhotonSharingData>(
        "PhotonSharingData photon_sharing_instance");

    // Set class member ptr to shared data
    p_this->shared_data_ = instance_.first;

    p_this->setupLocalElementsDiscreteIndex(shared_data_);
    p_this->setup_ = true;
  }

  // init particle
  openmc::SourceSite particle;

  // weight
  particle.particle = openmc::ParticleType::photon;

  // Currently multiplying by number of ranks, mimicing behavoir in
  // openmc/src/source.cpp:sample_external_source (line 696)
  particle.wgt = calculateParticleWeight(shared_data_) * openmc::mpi::n_procs;

  // Calculate position
  // Get element index of sampled element
  int32_t element_id = sampleLocalElementsIndex(seed);
  // Sample a position within the volume of our chosen chosen element
  particle.r = sampleElementVolume(seed, _mesh_id, element_id);
  // Sample an isotropic angle
  particle.u = angle_->sample(seed);
  // Sample isotropic angle
  particle.E = sampleElementEnergy(seed, shared_data_, element_id);
  // Sample time
  particle.time = time_->sample(seed);

  return particle;
}

extern "C" std::unique_ptr<FizzyCompiledSource>
openmc_create_source(std::string parameters) {
  int32_t mesh_id = std::stoi(parameters);
  return std::make_unique<FizzyCompiledSource>(mesh_id);
}
