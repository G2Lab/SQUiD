#pragma once

#include "../../src/databases/FHE_disk_database.hpp"

#include <drogon/HttpController.h>
#include <vector>
#include <helib/helib.h>
#include <string>
#include <sstream>

using namespace drogon;
using namespace std;


namespace api
{
namespace v1
{
class Server : public drogon::HttpController<Server>
{
  public:
    METHOD_LIST_BEGIN
    METHOD_ADD(Server::getContext,"/getContext?key={1}", Get);
    METHOD_ADD(Server::printDB,"/printDB?key={1}", Get);
    METHOD_ADD(Server::authorizeAPI,"/au?key={1}", Post);
    METHOD_ADD(Server::countingQueryAPI,"/countingQuery?query={1}&conj={2}&key={3}", Get);
    METHOD_ADD(Server::mafQueryAPI,"/mafQuery?query={1}&conj={2}&target={3}&key={4}", Get);
    METHOD_ADD(Server::PRSQueryAPI,"/PRSQuery?params={1}&key={2}", Get);
    METHOD_ADD(Server::similarityQueryAPI,"/simQuery?threshold={1}&binaryPheno={2}&key={3}", Post);
    METHOD_ADD(Server::getHeadersAPI,"/headers?key={1}", Get);
    METHOD_ADD(Server::getOwnerPublicKey,"/getOwnerPublicKey?key={1}", Get);
    METHOD_LIST_END

    Server();

    void getContext(const HttpRequestPtr &req,
                 std::function<void (const HttpResponsePtr &)> &&callback,
                 const std::string &apikey) const;
    
    void printDB(const HttpRequestPtr &req,
                 std::function<void (const HttpResponsePtr &)> &&callback,
                 const std::string &apikey) const;

    void authorizeAPI(const HttpRequestPtr &req,
                 std::function<void (const HttpResponsePtr &)> &&callback,
                 const std::string &apikey);
      
    void countingQueryAPI(const HttpRequestPtr &req,
                 std::function<void (const HttpResponsePtr &)> &&callback,
                 std::string query,
                 std::string conj,
                 const std::string &apikey) const;
    void mafQueryAPI(const HttpRequestPtr &req,
                 std::function<void (const HttpResponsePtr &)> &&callback,
                 std::string query,
                 std::string conj,
                 std::string target,
                 const std::string &apikey) const;
    void PRSQueryAPI(const HttpRequestPtr &req,
                 std::function<void (const HttpResponsePtr &)> &&callback,
                 std::string params,
                 const std::string &apikey) const;
    void similarityQueryAPI(const HttpRequestPtr &req,
                  std::function<void (const HttpResponsePtr &)> &&callback,
                  std::string threshold,
                  std::string binaryPheno,
                  const std::string &apikey) const;
    void getHeadersAPI(const HttpRequestPtr &req,
                 std::function<void (const HttpResponsePtr &)> &&callback,
                 const std::string &apikey) const;
    void getOwnerPublicKey(const HttpRequestPtr &req,
                  std::function<void (const HttpResponsePtr &)> &&callback,
                  const std::string &apikey) const;
    
    void addKSK(helib::PubKey& client_pk, string id);

  private:
    const std::string MasterApiKey = "nNCHuSdBWZsDJNFOJqUWDAUibEvVcVniRqbiIoM";
    std::unordered_set<std::string> api_keys;
    std::map<std::string, std::pair<std::vector<helib::DoubleCRT>, std::vector<helib::DoubleCRT>>> key_switch_store;
    FHEDiskDatabase squid;
};
}
}