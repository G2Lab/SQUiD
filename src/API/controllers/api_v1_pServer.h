#pragma once

#include "src/plain.h"

#include <drogon/HttpController.h>
#include <vector>
#include <string>
#include <sstream>

using namespace drogon;
using namespace std;


namespace api
{
namespace v1
{
class PServer : public drogon::HttpController<PServer>
{
  public:
  
    METHOD_LIST_BEGIN
    //METHOD_ADD(Server::countingQueryAPI,"/countingQuery?query={1}&conj={2}&key={3}", Get);
    //METHOD_ADD(Server::mafQueryAPI,"/mafQuery?query={1}&conj={2}&target={3}&key={4}", Get);
    //METHOD_ADD(Server::distributionQueryAPI,"/distQuery?params={1}&key={2}", Get);
    METHOD_LIST_END
  

    PServer();
      
    /*
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
    void distributionQueryAPI(const HttpRequestPtr &req,
                 std::function<void (const HttpResponsePtr &)> &&callback,
                 std::string params,
                 const std::string &apikey) const;
    void getHeadersAPI(const HttpRequestPtr &req,
                 std::function<void (const HttpResponsePtr &)> &&callback,
                 const std::string &apikey) const;
    */

  private:
    const std::string MasterApiKey = "nNCHuSdBWZsDJNFOJqUWDAUibEvVcVniRqbiIoM";
    std::unordered_set<std::string> api_keys;

    Plain plain;
};
}
}
