#include "FispactCompiledSource.hpp"

#include "mpi.h"
#include "openmc/capi.h"
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
#include <memory>
#include <stdexcept>
#include <string>
#include <sys/types.h>

FizzyCompiledSource::FizzyCompiledSource(int32_t mesh_id)
    : mesh_id_(mesh_id), angle_(openmc::UPtrAngle(new openmc::Isotropic)) {
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
}

const std::string FizzyCompiledSource::generateInterprocessName() {
  char mpi_proc_name[MPI_MAX_PROCESSOR_NAME];
  int len = 0;
  int err = MPI_Get_processor_name(mpi_proc_name, &len);

  MPI_Comm node_comm;
  int local_rank;
  MPI_Comm_split_type(openmc::mpi::intracomm, MPI_COMM_TYPE_SHARED, 0,
                      MPI_INFO_NULL, &node_comm);
  MPI_Comm_rank(node_comm, &local_rank);

  std::string ipc_name = std::string(mpi_proc_name);

  ipc_name += "_" + std::to_string(local_rank);

  return ipc_name;
}

FizzyCompiledSource::~FizzyCompiledSource() {
  // const std::string shared_data_name = generateInterprocessName();
  //  bi::shared_memory_object::remove(shared_data_name.c_str());
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

  const BoostIpVector &boost_element_strength = shared_data->_elem_strength;

  std::vector<double> element_strengths(boost_element_strength.begin(),
                                        boost_element_strength.end());

  element_ids_.clear();
  element_ids_.reserve(element_strengths.size());

  for (const auto &[global_elem_id, local_elem_id] :
       shared_data->_local_elem_idx_map) {
    element_ids_.push_back(global_elem_id);
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

  int32_t mesh_idx;

  try {
    mesh_idx = openmc::model::mesh_map.at(mesh_id);
  } catch (std::out_of_range) {
    std::cerr << "No mesh with mesh_idx " + std::to_string(mesh_id) + " found"
              << std::endl;
  }
  std::unique_ptr<openmc::Mesh> &mesh = openmc::model::meshes[mesh_idx];

  openmc::LibMesh *derived_libmesh_ptr =
      dynamic_cast<openmc::LibMesh *>(mesh.get());

  do {
    r = derived_libmesh_ptr->sample_element(element_id, seed);
  } while (!this->satisfies_spatial_constraints(r));
  return r;
}

size_t FizzyCompiledSource::getSpectraIdx(const PhotonSharingData *shared_data,
                                          const int32_t &element_id) const {
  size_t n_bins = shared_data->_photon_bins.size() - 1;
  return n_bins * shared_data->_local_elem_idx_map.at(element_id);
}

double
FizzyCompiledSource::sampleElementEnergy(uint64_t *seed,
                                         const PhotonSharingData *shared_data,
                                         const int32_t &element_id) const {

  // Set up ptrs to element energy distribution and photon bins in shared data,
  // this is a bit messy but avoids doing a copy!
  // const double *element_energy =
  //     &shared_data->_photon_fluxes.at(getSpectraIdx(shared_data,
  //     element_id));

  /// _photon_bins represents the boundaries of each bin, so 25 values represent
  /// 24 bins, hence the need for -1

  auto idx = shared_data->_local_elem_idx_map.at(element_id);
  const openmc::Tabular &energy_distribution = *energy_distributions_[idx];

  while (true) {
    double energy = energy_distribution.sample(seed).first;

    if (satisfies_energy_constraints(energy)) {
      return energy;
    }
  }
}

void FizzyCompiledSource::constructEnergyDistributions(
    const PhotonSharingData *shared_data) {

  energy_distributions_.reserve(shared_data->_local_elem_idx_map.size());
  for (auto &element_id : element_ids_) {
    size_t n_bin_boundaries = shared_data->_photon_bins.size();
    size_t n_bins = n_bin_boundaries - 1;

    std::vector<double> element_energy(
        shared_data->_photon_fluxes.begin() +
            getSpectraIdx(shared_data, element_id),
        shared_data->_photon_fluxes.begin() +
            getSpectraIdx(shared_data, element_id) + n_bins);

    const double *photon_bins = &shared_data->_photon_bins.at(0);

    /// Must divide element_energy by bin width to get proper histogram sampling
    for (int i = 1; i < n_bin_boundaries; i++) {
      double bin_width = *(photon_bins + i) - *(photon_bins + (i - 1));
      element_energy.at(i - 1) /= bin_width;
    }

    energy_distributions_.push_back(std::make_unique<openmc::Tabular>(
        photon_bins, element_energy.data(), n_bins,
        openmc::Interpolation::histogram, nullptr));
  }
}

void FizzyCompiledSource::sharedDataInit() const {

  if (initialised_) {
    return;
  }
  auto *p_this = const_cast<FizzyCompiledSource *>(this);
  p_this->instance_ =
      p_this->segment_.find<PhotonSharingData>("photon_sharing_instance");
  // Set class member ptr to shared data
  p_this->shared_data_ = instance_.first;
  p_this->initialised_ = true;
}

void FizzyCompiledSource::timestepInit() const {
  if (!shared_data_->_is_setup) {
    auto *p_this = const_cast<FizzyCompiledSource *>(this);
    p_this->setupLocalElementsDiscreteIndex(shared_data_);
    p_this->constructEnergyDistributions(shared_data_);
    p_this->strength_ = shared_data_->_total_domain_strength;

    std::vector<double> source_strengths = {strength_};
    openmc::model::external_sources_probability.assign(source_strengths);
    shared_data_->_is_setup = true;
  }
}

openmc::SourceSite FizzyCompiledSource::sample(uint64_t *seed) const {

  // If this is the first sample, do setup.
  // This setup should really be in the constructor, but doing it here makes it
  // compatible with MOOSE-multiapp runs. When using as a MOOSE-multiapp, the
  // constructor gets called before the shared data is ready to be read

  openmc::SourceSite particle;
// init particle
#pragma omp single
  {
    // Perform one time setup
    sharedDataInit();
    // Perform once per timestep setup
    timestepInit();
  }
#pragma omp barrier

  int32_t element_id = sampleLocalElementsIndex(seed);

  particle.particle = openmc::ParticleType::photon();

  // Currently multiplying by number of ranks, mimicing behavoir in
  // openmc/src/source.cpp:sample_external_source (line 696)
  particle.wgt = calculateParticleWeight(shared_data_) * openmc::mpi::n_procs;

  // Calculate position
  // Get element index of sampled element
  // Sample a position within the volume of our chosen chosen element
  particle.r = sampleElementVolume(seed, mesh_id_, element_id);
  // Sample an isotropic angle
  particle.u = angle_->sample(seed).first;
  // Sample isotropic angle
  particle.E = sampleElementEnergy(seed, shared_data_, element_id);
  // Sample time
  particle.time = time_->sample(seed).first;

  return particle;
}

extern "C" std::unique_ptr<FizzyCompiledSource>
openmc_create_source(std::string parameters) {
  int32_t mesh_id = std::stoi(parameters);
  return std::make_unique<FizzyCompiledSource>(mesh_id);
}
