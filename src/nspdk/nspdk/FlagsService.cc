#include "nspdk/FlagsService.h"


namespace nspdk {

	FlagsService * FlagsService::The_FlagsService = NULL;

FlagsServiceClient::FlagsServiceClient(const std::string& id_for_flags_service) {
  this->id_for_flags_service = id_for_flags_service;
  FlagsService::get_instance().register_flags_service_client(this);
}

FlagsServiceClient::~FlagsServiceClient() {
  for (FMap::const_iterator i=flags_traits.begin(); i!=flags_traits.end(); ++i)
    delete i->second;
  FlagsService::get_instance().unregister_flags_service_client(this);
}

} // namespace
