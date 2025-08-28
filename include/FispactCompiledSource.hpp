#pragma once

//
#include <boost/interprocess/interprocess_fwd.hpp>
#include <cstdint>
#include <memory> // for unique_ptr
#include <sys/types.h>

// Define LIBMESH so that openmc libmesh gets included
#define LIBMESH

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

  double calculateParticleWeight(const PhotonSharingData *shared_data) const;

  int32_t sampleLocalElementsIndex(uint64_t *seed) const;

  void setupLocalElementsDiscreteIndex(const PhotonSharingData *shared_data);

  openmc::Position sampleElementVolume(u_int64_t *seed, int mesh_id,
                                       int element_id) const;

  double sampleElementEnergy(uint64_t *seed,
                             const PhotonSharingData *shared_data,
                             const int32_t &element_id) const;

  bool constraints_applied() const override { return true; }

  openmc::SourceSite sample(u_int64_t *seed) const;

  // Data members
  openmc::UPtrAngle angle_;
  openmc::UPtrDist time_;
  //
  int _num_ranks = -1;
  int _my_rank = -1;
  int32_t _mesh_id;

  // Interprocess data structures
  boost::interprocess::managed_shared_memory segment_;
  std::pair<PhotonSharingData *, std::size_t> instance_;
  const PhotonSharingData *shared_data_;

  std::vector<int> element_ids_;
  openmc::DiscreteIndex di_;

  bool setup_ = false;
};
