#pragma once

//
#include <boost/interprocess/interprocess_fwd.hpp>
#include <cstdint>
#include <memory> // for unique_ptr
#include <sys/types.h>

#include "openmc/distribution_multi.h"
#include "openmc/position.h"
#include "utils/PhotonSharingData.h"

#include "openmc/particle.h"
#include "openmc/particle_data.h"
#include "openmc/random_lcg.h"
#include "openmc/source.h"

class FizzyCompiledSource : public openmc::Source {
public:
  FizzyCompiledSource() {};

  FizzyCompiledSource(int32_t mesh_id);

  ~FizzyCompiledSource();

  void sampleElementStrengths();
  //
  double calculateParticleWeight(const double &local_domain_strength,
                                 const double &total_domain_strength) const;

  double calculateParticleWeight(const PhotonSharingData *shared_data,
                                 bool uniform, int32_t element_id,
                                 int32_t mesh_id) const;

  double
  calculateParticleWeightUniform(const double &element_strength,
                                 const double &element_volume,
                                 const double &total_domain_strength) const;

  int32_t sampleLocalElementsIndex(uint64_t *seed) const;

  void setupLocalElementsDiscreteIndex(const PhotonSharingData *shared_data,
                                       bool uniform, int32_t mesh_id);

  void constructEnergyDistributions(const PhotonSharingData *shared_data);

  openmc::Position sampleElementVolume(u_int64_t *seed, int mesh_id,
                                       int32_t element_id) const;

  double sampleElementEnergy(uint64_t *seed,
                             const PhotonSharingData *shared_data,
                             const int32_t &element_id) const;

  size_t getSpectraIdx(const PhotonSharingData *shared_data,
                       const int32_t &element_id) const;

  void sharedDataInit() const;

  void timestepInit() const;

  bool constraints_applied() const override { return true; }

  const std::string generateInterprocessName();

  openmc::SourceSite sample(uint64_t *seed) const;

  // Data members
  openmc::UPtrAngle angle_;
  openmc::UPtrDist time_;

  //
  int32_t mesh_id_;

  // Interprocess data structures
  boost::interprocess::managed_shared_memory segment_;
  std::pair<PhotonSharingData *, std::size_t> instance_;
  PhotonSharingData *shared_data_;

  std::vector<int> element_ids_;
  std::vector<std::unique_ptr<openmc::Tabular>> energy_distributions_;
  openmc::DiscreteIndex di_;

  bool uniform_ = false;
  bool initialised_ = false;

  double rank_elements_volume_ = 0;
};
