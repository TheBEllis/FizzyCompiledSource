#include "FispactCompiledSource.hpp"
#include "mpi.h"
#include "openmc/constants.h"
#include "openmc/distribution.h"
#include "openmc/distribution_multi.h"
#include "openmc/mesh.h"
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

  // Get MPI rank
  auto return_err = MPI_Comm_rank(MPI_COMM_WORLD, &_my_rank);
  if (return_err != MPI_SUCCESS) {
    std::cerr << "MPI_FAILURE" << std::endl;
  }
  return_err = MPI_Comm_size(MPI_COMM_WORLD, &_num_ranks);
  if (return_err != MPI_SUCCESS) {
    std::cerr << "MPI_FAILURE" << std::endl;
  }

  // Get shared interprocess data
  std::string shared_data_name = "SHARING_DATA_" + std::to_string(_my_rank);
  try {
    segment_ = boost::interprocess::managed_shared_memory(
        boost::interprocess::open_only, shared_data_name.c_str());
  } catch (bi::interprocess_exception) {
    std::cerr << "Could not find interprocess segment with name "
              << shared_data_name << std::endl;

    exit(-1);
  }
  instance_ = segment_.find<PhotonSharingData>(
      "PhotonSharingData photon_sharing_instance");

  // Set class member ptr
  shared_data_ = instance_.first;

  // strength_ = calculateParticleWeight(shared_data);
  // strength_ = shared_data->_total_domain_strength;
  strength_ = 1;

  setupLocalElementsDiscreteIndex(shared_data_);
}

FizzyCompiledSource::~FizzyCompiledSource() {
  std::string shared_data_name = "SHARING_DATA_" + std::to_string(_my_rank);
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

  openmc::span probs(element_strengths);
  di_ = openmc::DiscreteIndex(probs);
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
  const auto &mesh = openmc::model::meshes[mesh_idx];
  do {
    r = mesh->sample_element(element_id, seed);
  } while (!this->satisfies_spatial_constraints(r));
  return r;
}

double
FizzyCompiledSource::sampleElementEnergy(uint64_t *seed,
                                         const PhotonSharingData *shared_data,
                                         const int32_t &element_id) const {
  const double *p = &shared_data->_photon_fluxes.at(element_id).at(0);

  const double *x = &shared_data->_photon_bins.at(0);
  size_t n_bins = shared_data->_photon_bins.size();
  openmc::Tabular distribution(x, p, n_bins, openmc::Interpolation::histogram,
                               nullptr);

  while (true) {
    double energy = distribution.sample(seed);

    if (satisfies_energy_constraints(energy)) {
      return energy;
    }
  }
}

openmc::SourceSite FizzyCompiledSource::sample(uint64_t *seed) const {
  // init particle
  openmc::SourceSite particle;

  // A bit messy to do setup here, but const can't do it in constructor, as
  // the constructor is called before the shared memory is ready to be read
  // Also trying to put this in a function that returns a PhotonSharingData*
  // or something like that doesn't work and segfaults

  // Get boost::interprocess shared memory
  // std::string shared_data_name = "SHARING_DATA_" + std::to_string(_my_rank);
  // boost::interprocess::managed_shared_memory segment;
  // try {
  //   segment = boost::interprocess::managed_shared_memory(
  //       boost::interprocess::open_only, shared_data_name.c_str());
  // } catch (bi::interprocess_exception) {
  //   std::cerr << "Could not find interprocess segment with name "
  //             << shared_data_name << std::endl;
  //
  //   exit(-1);
  // }

  // std::pair<PhotonSharingData *, std::size_t> instance;
  // instance = segment.find<PhotonSharingData>(
  //     "PhotonSharingData photon_sharing_instance");
  // const PhotonSharingData *shared_data = instance.first;

  // weight
  particle.particle = openmc::ParticleType::photon;
  particle.wgt = calculateParticleWeight(shared_data_) * _num_ranks;
  // position
  // Get element index of sampled element
  int32_t element_id = sampleLocalElementsIndex(seed);
  // Sample within chosen elements volume
  particle.r = sampleElementVolume(seed, _mesh_id, element_id);
  // angle
  particle.u = angle_->sample(seed);
  // energy
  particle.E = sampleElementEnergy(seed, shared_data_, element_id);
  particle.time = time_->sample(seed);

  return particle;
}

extern "C" std::unique_ptr<FizzyCompiledSource>
openmc_create_source(std::string parameters) {
  int32_t mesh_id = std::stoi(parameters);
  return std::make_unique<FizzyCompiledSource>(mesh_id);
}
