#include "api_v1_pServer.h"

using namespace api::v1;
using namespace std;

PServer::PServer(): plain(){
    api_keys = std::unordered_set<std::string>{MasterApiKey};
}
